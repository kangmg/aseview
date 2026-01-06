"""
Command-line interface for aseview molecular viewer.

Usage:
    aseview2 molecule.xyz
    aseview2 trajectory.xyz -i 0:10
    aseview2 structure.cif -p 8888
"""
import argparse
import http.server
import socketserver
import threading
import webbrowser
import sys
import os


def parse_index(index_str: str):
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

    Args:
        html_content: HTML string to serve
        port: Port number to bind to
        open_browser: Whether to open browser automatically
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
                print(f"Error: Could not find available port (tried {original_port}-{port})")
                sys.exit(1)

    url = f"http://localhost:{port}"

    print(f"\n  aseview2 server running at: {url}")
    print(f"  (for SSH: ssh -L {port}:localhost:{port} ...)")
    print(f"\n  Press Ctrl+C to stop\n")

    # Open browser if requested and not in SSH session
    if open_browser and not os.environ.get('SSH_CONNECTION'):
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()

    # Run server
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        print("\n  Shutting down...")
        httpd.server_close()
        # Cleanup temp files
        try:
            os.remove(html_path)
            os.rmdir(temp_dir)
        except:
            pass


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="aseview2",
        description="Molecular structure viewer for ASE-supported file formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  aseview2 molecule.xyz              # View single structure
  aseview2 trajectory.xyz -i :       # View all frames as animation
  aseview2 trajectory.xyz -i 0:10    # View frames 0-9
  aseview2 trajectory.xyz -i -1      # View last frame
  aseview2 structure.cif -p 9000     # Use custom port
  aseview2 file1.xyz file2.xyz       # Overlay multiple structures

SSH port forwarding:
  On remote server: aseview2 molecule.xyz -p 8080
  On local machine: ssh -L 8080:localhost:8080 user@remote
  Then open http://localhost:8080 in your browser
        """
    )

    parser.add_argument(
        "files",
        nargs="+",
        help="Input file(s) in any ASE-supported format (xyz, cif, pdb, etc.)"
    )

    parser.add_argument(
        "-i", "--index",
        type=str,
        default=None,
        metavar="INDEX",
        help="Index or slice for trajectory (e.g., 0, -1, :, 0:10, ::2)"
    )

    parser.add_argument(
        "-f", "--format",
        type=str,
        default=None,
        metavar="FORMAT",
        help="File format (auto-detected if not specified)"
    )

    parser.add_argument(
        "-p", "--port",
        type=int,
        default=8080,
        metavar="PORT",
        help="Port for HTTP server (default: 8080)"
    )

    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Don't open browser automatically"
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        metavar="FILE",
        help="Save HTML to file instead of serving"
    )

    parser.add_argument(
        "-v", "--viewer",
        type=str,
        choices=["molecular", "normal", "overlay"],
        default=None,
        help="Viewer type (auto-selected based on input)"
    )

    parser.add_argument(
        "--style",
        type=str,
        choices=["default", "cartoon", "neon", "glossy", "metallic", "rowan", "grey"],
        default="cartoon",
        help="Visual style (default: cartoon)"
    )

    args = parser.parse_args()

    # Import ASE for reading files
    try:
        from ase.io import read
    except ImportError:
        print("Error: ASE is required. Install with: pip install ase")
        sys.exit(1)

    # Import viewer classes
    from .wrapper import MolecularViewer, NormalViewer, OverlayViewer

    # Parse index
    index = parse_index(args.index)

    # Read structures
    all_atoms = []
    for filepath in args.files:
        if not os.path.exists(filepath):
            print(f"Error: File not found: {filepath}")
            sys.exit(1)

        try:
            atoms = read(filepath, index=index, format=args.format)
            if isinstance(atoms, list):
                all_atoms.extend(atoms)
            else:
                all_atoms.append(atoms)
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            sys.exit(1)

    if not all_atoms:
        print("Error: No structures loaded")
        sys.exit(1)

    # Determine viewer type
    viewer_type = args.viewer
    if viewer_type is None:
        if len(args.files) > 1:
            viewer_type = "overlay"
        elif len(all_atoms) > 1:
            viewer_type = "molecular"  # Use molecular viewer for trajectories
        else:
            viewer_type = "molecular"

    # Create viewer
    viewer_kwargs = {"style": args.style}

    if viewer_type == "molecular":
        viewer = MolecularViewer(all_atoms, **viewer_kwargs)
    elif viewer_type == "normal":
        viewer = NormalViewer(all_atoms, **viewer_kwargs)
    elif viewer_type == "overlay":
        viewer = OverlayViewer(all_atoms, **viewer_kwargs)

    # Get HTML content
    html_content = viewer.get_html()

    # Output
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Saved to {args.output}")
    else:
        serve_html(html_content, port=args.port, open_browser=not args.no_browser)


if __name__ == "__main__":
    main()
