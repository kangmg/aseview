#!/usr/bin/env python3
"""
Render ASE-supported structure file to ASEView HTML for VS Code custom editor.
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional


def parse_index(index_str: Optional[str]):
    """Parse ASE-style index string."""
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
        if len(parts) == 3:
            start = int(parts[0]) if parts[0] else None
            stop = int(parts[1]) if parts[1] else None
            step = int(parts[2]) if parts[2] else None
            return slice(start, stop, step)

    return int(index_str)


def build_viewer(viewer_name: str, atoms_list, style: str):
    from aseview.wrapper import FragSelector, MolecularViewer, NormalViewer, OverlayViewer

    if viewer_name == "frag":
        return FragSelector(atoms_list[0], style=style)
    if viewer_name == "overlay":
        return OverlayViewer(atoms_list, style=style)
    if viewer_name == "normal":
        return NormalViewer(atoms_list[0], style=style)

    # molecular and auto fallback
    return MolecularViewer(atoms_list, style=style)


def main() -> int:
    parser = argparse.ArgumentParser(description="Render file with ASEView and print HTML.")
    parser.add_argument("--file", required=True, help="Path to structure file.")
    parser.add_argument("--format", default=None, help="Optional ASE format override.")
    parser.add_argument("--index", default=":", help="ASE index/slice expression.")
    parser.add_argument(
        "--viewer",
        default="auto",
        choices=["auto", "molecular", "overlay", "frag", "normal"],
        help="Viewer mode.",
    )
    parser.add_argument("--style", default="cartoon", help="ASEView style.")
    args = parser.parse_args()

    try:
        from ase.io import read

        idx = parse_index(args.index)
        atoms = read(args.file, index=idx, format=args.format)
        atoms_list = atoms if isinstance(atoms, list) else [atoms]
        if not atoms_list:
            raise ValueError("No structure frames were loaded.")

        viewer_name = args.viewer
        if viewer_name == "auto":
            viewer_name = "molecular"

        viewer = build_viewer(viewer_name, atoms_list, args.style)
        html = viewer.get_html()
        sys.stdout.write(html)
        return 0
    except Exception as exc:
        # Keep error text in stderr so extension can display it cleanly.
        sys.stderr.write(f"{type(exc).__name__}: {exc}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
