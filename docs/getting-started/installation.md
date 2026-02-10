# Installation

## Requirements

- Python 3.8+
- ASE (Atomic Simulation Environment)

## Install from PyPI

```bash
pip install aseview2
```

## Install from Source

```bash
git clone https://github.com/kangmg/aseview_v2_dev.git
cd aseview_v2_dev
pip install -e .
```

## Dependencies

aseview2 will automatically install the following dependencies:

| Package | Purpose |
|---------|---------|
| `ase` | Atomic structure handling |
| `numpy` | Numerical operations |
| `typer` | CLI framework |
| `rich` | Terminal formatting |

## Verify Installation

```bash
# Check CLI
aseview2 --help

# Check Python import
python -c "from aseview import MolecularViewer; print('OK')"
```
