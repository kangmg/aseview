# Installation

## Requirements

- Python 3.9+
- ASE (Atomic Simulation Environment)

## Install from PyPI

```bash
pip install aseview
```

## Install from Source

```bash
git clone https://github.com/kangmg/aseview.git
cd aseview
pip install -e .
```

## Dependencies

aseview will automatically install the following dependencies:

| Package | Purpose |
|---------|---------|
| `ase` | Atomic structure handling |
| `numpy` | Numerical operations |
| `typer` | CLI framework |
| `rich` | Terminal formatting |

## Verify Installation

```bash
# Check CLI
aseview --help

# Check Python import
python -c "from aseview import MolecularViewer; print('OK')"
```

## Note on aseview2

`aseview2` is kept as a backward-compatible alias for `aseview`. Both commands are identical:

```bash
aseview molecule.xyz   # recommended
aseview2 molecule.xyz  # alias (backward compatibility)
```
