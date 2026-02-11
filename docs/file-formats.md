# Supported File Formats

ASE automatically detects file formats using three methods:

1. **Extension** - File extension (e.g., `.xyz`, `.cif`)
2. **Filename pattern** - Glob matching (e.g., `*POSCAR*`, `GEOMETRY.OUT`)
3. **Magic bytes** - Content signature detection

## Format Detection Priority

```
1. Filename pattern match (glob)
2. File extension match
3. Magic bytes detection (reads first ~50KB)
4. Fallback to extension as format name
```

## Quick Reference

| Format | Extension | Filename Pattern | Magic Detection |
|--------|-----------|------------------|-----------------|
| `extxyz` | `.xyz` | — | — |
| `cif` | `.cif` | — | — |
| `vasp` | `.poscar` | `*POSCAR*`, `*CONTCAR*` | — |
| `traj` | `.traj` | — | ✓ |
| `json` | `.json` | — | — |
| `pdb` | `.pdb` | — | — |

---

## All Formats by Category

### Molecular Formats

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `extxyz` | xyz | — | — | Extended XYZ file |
| `xyz` | — | — | — | XYZ file (plain) |
| `cif` | cif | — | — | Crystallographic Information File |
| `proteindatabank` | pdb | — | — | Protein Data Bank |
| `mol` | — | — | — | MDL Molfile |
| `sdf` | — | — | — | SDF format |
| `cjson` | cjson | — | — | Chemical JSON |
| `gen` | — | — | — | DFTB+ GEN format |

### ASE Native

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `traj` | traj | — | `- of UlmASE-Trajectory` | ASE trajectory |
| `json` | json | — | — | ASE JSON database |
| `db` | — | — | — | ASE SQLite database |
| `bundletrajectory` | — | — | — | ASE bundle trajectory |

### VASP

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `vasp` | poscar | `*POSCAR*`, `*CONTCAR*`, `*CENTCAR*` | — | POSCAR/CONTCAR |
| `vasp-out` | — | `*OUTCAR*` | — | OUTCAR file |
| `vasp-xdatcar` | — | `*XDATCAR*` | — | XDATCAR file |
| `vasp-xml` | — | `*vasp*.xml` | — | vasprun.xml |

### Quantum ESPRESSO

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `espresso-in` | pwi | — | `*\n&system`, `*\n&SYSTEM` | QE input |
| `espresso-out` | pwo, out | — | `*Program PWSCF` | QE output |

### Gaussian

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `gaussian-in` | com, gjf | — | — | Gaussian input |
| `gaussian-out` | log | — | `*Entering Gaussian System` | Gaussian output |

### ORCA

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `orca-output` | — | — | `* O   R   C   A *` | ORCA output |

### FHI-aims

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `aims` | in | — | — | FHI-aims geometry |
| `aims-output` | — | — | `*Invoking FHI-aims ...` | FHI-aims output |

### CASTEP

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `castep-castep` | castep | — | — | CASTEP output |
| `castep-cell` | cell | — | — | CASTEP geometry |
| `castep-geom` | geom | — | — | CASTEP trajectory |
| `castep-md` | md | — | — | CASTEP MD file |
| `castep-phonon` | phonon | — | — | CASTEP phonon |

### CP2K

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `cp2k-dcd` | dcd | — | — | CP2K DCD file |
| `cp2k-restart` | restart | — | — | CP2K restart |

### GPAW

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `gpaw-out` | — | — | `*  ___ ___ ___ _ _ _` | GPAW text output |
| `gpw` | — | — | `- of UlmGPAW` | GPAW restart file |

### ABINIT

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `abinit-gsr` | — | `*o_GSR.nc` | — | ABINIT GSR file |
| `abinit-in` | — | — | `*znucl *` | ABINIT input |
| `abinit-out` | — | — | `*.Version * of ABINIT` | ABINIT output |

### LAMMPS

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `lammps-data` | — | — | — | LAMMPS data file |
| `lammps-dump-text` | — | — | `ITEM: TIMESTEP` (regex) | LAMMPS text dump |
| `lammps-dump-binary` | — | — | — | LAMMPS binary dump |

### NWChem

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `nwchem-in` | nwi | — | — | NWChem input |
| `nwchem-out` | nwo | — | `*Northwest Computational...` | NWChem output |

### TURBOMOLE

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `turbomole` | — | `coord` | `$coord` | TURBOMOLE coord |
| `turbomole-gradient` | — | `gradient` | `$grad` | TURBOMOLE gradient |

### Other DFT Codes

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `dftb` | — | — | `Geometry` | DFTB+ input |
| `elk` | — | `GEOMETRY.OUT` | — | ELK output |
| `exciting` | — | `input.xml`, `INFO.out` | — | exciting |
| `siesta-xv` | — | `*.XV` | — | Siesta XV file |
| `octopus-in` | — | `inp` | — | Octopus input |
| `onetep-out` | — | — | `*Linear-Scaling Ab Initio*` | ONETEP output |
| `qbox` | — | — | `*:simulation xmlns:` | QBOX output |

### Crystal Formats

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `crystal` | f34, 34 | `f34`, `34` | — | Crystal fort.34 |
| `xsf` | — | — | `*\nCRYSTAL`, `*\nATOMS`, etc. | XCrySDen |
| `cube` | cube | — | — | Gaussian CUBE |
| `struct` | — | — | — | WIEN2k structure |
| `res` | shelx | — | — | SHELX format |

### Visualization/Export

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `png` | — | — | — | PNG image |
| `eps` | — | — | — | Encapsulated PostScript |
| `pov` | — | — | — | POV-Ray |
| `gif` | — | — | — | GIF animation |
| `mp4` | — | — | — | MP4 video |
| `html` | — | — | — | X3DOM HTML |
| `x3d` | — | — | — | X3D format |

### Other

| Format | Ext | Glob | Magic | Description |
|--------|-----|------|-------|-------------|
| `cfg` | — | — | — | AtomEye configuration |
| `gromacs` | gro | — | — | Gromacs coordinates |
| `gromos` | g96 | — | — | Gromos96 geometry |
| `dlp4` | config | `*CONFIG*` | — | DL_POLY_4 CONFIG |
| `dlp-history` | — | `HISTORY` | — | DL_POLY HISTORY |
| `eon` | con | — | — | EON CON file |
| `gpumd` | — | `xyz.in` | — | GPUMD input |
| `magres` | — | — | — | MAGRES NMR data |
| `netcdftrajectory` | — | — | `CDF` | AMBER NetCDF |
| `v-sim` | ascii | — | — | V_Sim ascii |
| `rmc6f` | rmc6f | — | — | RMCProfile |

---

## Compression Support

ASE automatically handles compressed files:

| Extension | Compression |
|-----------|-------------|
| `.gz` | gzip |
| `.bz2` | bzip2 |
| `.xz` | lzma |

```python
# These all work automatically:
atoms = read("structure.xyz.gz")
atoms = read("trajectory.traj.bz2")
write("output.cif.xz", atoms)
```

---

## Format Codes

Each format has a 2-character code:

| 1st char | Meaning |
|----------|---------|
| `1` | Single Atoms object |
| `+` | Multiple Atoms objects (trajectory) |

| 2nd char | Meaning |
|----------|---------|
| `F` | Accepts file descriptor |
| `B` | Binary mode |
| `S` | Requires filename string |

---

## Specifying Format Manually

```python
from ase.io import read, write

# Auto-detect (default)
atoms = read("structure.xyz")

# Explicit format
atoms = read("myfile.dat", format="extxyz")
atoms = read("POSCAR", format="vasp")

# Write with format
write("output", atoms, format="vasp")
```

---

## Programmatic Access

```python
from ase.io.formats import ioformats

# List all formats
for name, fmt in ioformats.items():
    print(f"{name}: ext={fmt.extensions}, glob={fmt.globs}")

# Get specific format info
fmt = ioformats['vasp']
print(fmt.full_description())
```
