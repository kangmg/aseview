"""
Command-line interface for aseview molecular viewer.

Usage:
    aseview2 molecule.xyz
    aseview2 trajectory.xyz -i 0:10
    aseview2 structure.cif -p 8888
"""
import http.server
import socketserver
import threading
import webbrowser
import sys
import os
import signal
from typing import Optional, List

import typer
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich import print as rprint

app = typer.Typer(
    name="aseview2",
    help="Molecular structure viewer for ASE-supported file formats",
    add_completion=False,
    rich_markup_mode="rich",
)
console = Console()


def kill_port(port: int) -> bool:
    """
    Kill process using the specified port.

    Returns True if a process was killed, False otherwise.
    """
    import socket

    # Check if port is in use first
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        sock.bind(('', port))
        sock.close()
        return False  # Port is free
    except OSError:
        pass  # Port is in use

    # Find PID using the port by parsing /proc/net/tcp
    pids_to_kill = set()

    try:
        # Read /proc/net/tcp to find socket inode
        with open('/proc/net/tcp', 'r') as f:
            lines = f.readlines()[1:]  # Skip header

        target_inodes = set()
        for line in lines:
            parts = line.split()
            if len(parts) >= 10:
                local_addr = parts[1]
                # local_addr format: IP:PORT (hex)
                local_port = int(local_addr.split(':')[1], 16)
                if local_port == port:
                    inode = parts[9]
                    target_inodes.add(inode)

        if not target_inodes:
            return False

        # Find PIDs that have these inodes
        for pid in os.listdir('/proc'):
            if not pid.isdigit():
                continue

            fd_dir = f'/proc/{pid}/fd'
            try:
                for fd in os.listdir(fd_dir):
                    try:
                        link = os.readlink(f'{fd_dir}/{fd}')
                        if link.startswith('socket:['):
                            inode = link[8:-1]
                            if inode in target_inodes:
                                pids_to_kill.add(int(pid))
                    except (OSError, PermissionError):
                        continue
            except (OSError, PermissionError):
                continue

        # Kill the processes
        for pid in pids_to_kill:
            try:
                os.kill(pid, signal.SIGKILL)
                console.print(f"  [yellow]Killed process {pid} on port {port}[/yellow]")
            except (OSError, PermissionError) as e:
                console.print(f"  [red]Failed to kill process {pid}: {e}[/red]")

        return len(pids_to_kill) > 0

    except Exception as e:
        console.print(f"  [red]Error finding process: {e}[/red]")
        return False


def parse_index(index_str: Optional[str]):
    """
    Parse ASE-style index string.

    Examples:
        "0" -> 0
        ":" -> ":"
        "0:10" -> slice(0, 10)
        "-1" -> -1
        "::2" -> slice(None, None, 2)
    """
    if index_str is None:
        return ":"

    index_str = index_str.strip()

    if index_str == ":":
        return ":"

    if ":" in index_str:
        parts = index_str.split(":")
        if len(parts) == 2:
            start = int(parts[0]) if parts[0] else None
            stop = int(parts[1]) if parts[1] else None
            return slice(start, stop)
        elif len(parts) == 3:
            start = int(parts[0]) if parts[0] else None
            stop = int(parts[1]) if parts[1] else None
            step = int(parts[2]) if parts[2] else None
            return slice(start, stop, step)

    return int(index_str)


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    """HTTP handler that suppresses logging."""

    def __init__(self, *args, directory=None, **kwargs):
        self.directory = directory
        super().__init__(*args, directory=directory, **kwargs)

    def log_message(self, format, *args):
        pass  # Suppress logging


def serve_html(html_content: str, port: int = 8080, open_browser: bool = True) -> None:
    """
    Serve HTML content via a simple HTTP server.
    """
    import tempfile

    # Create temporary directory with HTML file
    temp_dir = tempfile.mkdtemp()
    html_path = os.path.join(temp_dir, "index.html")

    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    # Find available port
    original_port = port
    max_attempts = 100

    for attempt in range(max_attempts):
        try:
            handler = lambda *args, **kwargs: QuietHandler(*args, directory=temp_dir, **kwargs)
            socketserver.TCPServer.allow_reuse_address = True
            httpd = socketserver.TCPServer(("", port), handler)
            break
        except OSError:
            port += 1
            if attempt == max_attempts - 1:
                console.print(f"[red]Error: Could not find available port (tried {original_port}-{port})[/red]")
                sys.exit(1)

    url = f"http://localhost:{port}"

    # Print beautiful banner
    console.print()

    link_text = Text()
    link_text.append("➜ ", style="green bold")
    link_text.append("Open in browser: ", style="white")
    link_text.append(url, style="cyan bold underline link " + url)

    panel = Panel(
        link_text,
        title="[bold green]aseview2 server running[/bold green]",
        subtitle="[dim]Press Ctrl+C to stop[/dim]",
        border_style="green",
        padding=(1, 2),
    )
    console.print(panel)

    console.print(f"  [dim]SSH forwarding: ssh -L {port}:localhost:{port} ...[/dim]")
    console.print()

    # Open browser if requested and not in SSH session
    if open_browser and not os.environ.get('SSH_CONNECTION'):
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()

    # Run server
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        console.print("\n  [yellow]Shutting down...[/yellow]")
        httpd.server_close()
        # Cleanup temp files
        try:
            os.remove(html_path)
            os.rmdir(temp_dir)
        except:
            pass


@app.command()
def main(
    files: List[str] = typer.Argument(
        ...,
        help="Input file(s) in any ASE-supported format (xyz, cif, pdb, etc.)",
    ),
    index: Optional[str] = typer.Option(
        None, "-i", "--index",
        help="Index or slice for trajectory (e.g., 0, -1, :, 0:10, ::2)",
    ),
    format: Optional[str] = typer.Option(
        None, "-f", "--format",
        help="File format (auto-detected if not specified)",
    ),
    port: int = typer.Option(
        8080, "-p", "--port",
        help="Port for HTTP server",
    ),
    no_browser: bool = typer.Option(
        False, "--no-browser",
        help="Don't open browser automatically",
    ),
    output: Optional[str] = typer.Option(
        None, "-o", "--output",
        help="Save HTML to file instead of serving",
    ),
    viewer: Optional[str] = typer.Option(
        None, "-v", "--viewer",
        help="Viewer type: molecular, normal, overlay (auto-selected)",
    ),
    style: str = typer.Option(
        "cartoon", "--style",
        help="Visual style: default, cartoon, neon, glossy, metallic, rowan, grey",
    ),
    kill: bool = typer.Option(
        False, "-k", "--kill",
        help="Kill existing process on the port before starting",
    ),
):
    """
    Molecular structure viewer for ASE-supported file formats.

    \b
    Examples:
      aseview2 molecule.xyz              # View single structure
      aseview2 trajectory.xyz -i :       # View all frames
      aseview2 trajectory.xyz -i 0:10    # View frames 0-9
      aseview2 molecule.xyz -k           # Kill existing server, then start

    \b
    SSH port forwarding:
      Server:  aseview2 molecule.xyz -p 8080
      Local:   ssh -L 8080:localhost:8080 user@remote
      Browser: http://localhost:8080
    """
    # Kill existing process on port if requested
    if kill:
        if kill_port(port):
            import time
            time.sleep(0.5)  # Wait for port to be released

    # Import ASE for reading files
    try:
        from ase.io import read
    except ImportError:
        console.print("[red]Error: ASE is required. Install with: pip install ase[/red]")
        raise typer.Exit(1)

    # Import viewer classes
    from .wrapper import MolecularViewer, NormalViewer, OverlayViewer

    # Parse index
    idx = parse_index(index)

    # Read structures
    all_atoms = []
    for filepath in files:
        if not os.path.exists(filepath):
            console.print(f"[red]Error: File not found: {filepath}[/red]")
            raise typer.Exit(1)

        try:
            atoms = read(filepath, index=idx, format=format)
            if isinstance(atoms, list):
                all_atoms.extend(atoms)
            else:
                all_atoms.append(atoms)
        except Exception as e:
            console.print(f"[red]Error reading {filepath}: {e}[/red]")
            raise typer.Exit(1)

    if not all_atoms:
        console.print("[red]Error: No structures loaded[/red]")
        raise typer.Exit(1)

    # Show loading info
    n_atoms = len(all_atoms[0]) if all_atoms else 0
    n_frames = len(all_atoms)
    console.print(f"  [green]✓[/green] Loaded [bold]{n_frames}[/bold] frame(s), [bold]{n_atoms}[/bold] atoms")

    # Determine viewer type
    viewer_type = viewer
    if viewer_type is None:
        if len(files) > 1:
            viewer_type = "overlay"
        elif len(all_atoms) > 1:
            viewer_type = "molecular"
        else:
            viewer_type = "molecular"

    # Create viewer
    viewer_kwargs = {"style": style}

    if viewer_type == "molecular":
        viewer_obj = MolecularViewer(all_atoms, **viewer_kwargs)
    elif viewer_type == "normal":
        viewer_obj = NormalViewer(all_atoms, **viewer_kwargs)
    elif viewer_type == "overlay":
        viewer_obj = OverlayViewer(all_atoms, **viewer_kwargs)
    else:
        console.print(f"[red]Error: Unknown viewer type: {viewer_type}[/red]")
        raise typer.Exit(1)

    # Get HTML content
    html_content = viewer_obj.get_html()

    # Output
    if output:
        with open(output, 'w', encoding='utf-8') as f:
            f.write(html_content)
        console.print(f"  [green]✓[/green] Saved to [bold]{output}[/bold]")
    else:
        serve_html(html_content, port=port, open_browser=not no_browser)


def cli():
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    cli()
