"""Pretvori Julia percent ćelije u Jupyter notebook, bez dodatnih paketa."""

import json
import sys
from pathlib import Path


root = Path(__file__).resolve().parent.parent
cells = []
kind = None
lines = []


def flush_cell():
    if kind is None:
        return
    source = "".join(lines).strip("\n") + "\n"
    cell = {"cell_type": kind, "metadata": {}, "source": source.splitlines(True)}
    if kind == "code":
        cell.update(execution_count=None, outputs=[])
    cells.append(cell)


source_path = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "primjer.jl"
output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else root / "Sorensen-Zhou-Lyap.ipynb"

for line in source_path.read_text().splitlines(True):
    if line.startswith("# %%"):
        flush_cell()
        kind = "markdown" if "[markdown]" in line else "code"
        lines = []
    elif kind == "markdown":
        lines.append(line[2:] if line.startswith("# ") else line.removeprefix("#"))
    else:
        lines.append(line)
flush_cell()

notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Julia Hammarling", "language": "julia", "name": "julia-hammarling"},
        "language_info": {"name": "julia", "file_extension": ".jl"},
    },
    "nbformat": 4,
    "nbformat_minor": 4,
}
output_path.write_text(
    json.dumps(notebook, ensure_ascii=False, indent=1) + "\n"
)
