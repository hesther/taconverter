taconverter
==============================

Converter for input files for the machine learning architectures [VASP ML](https://www.vasp.at/wiki/index.php/ML_AB), [N2P2](https://github.com/CompPhysVienna/n2p2) and models reading ase XYZ files (e.g. [MACE](https://github.com/ACEsuit/mace), developed within the [special research programme TACO](https://sfb-taco.at)

### Use

```
from taconverter import *
structures = read_N2P2("example_N2P2.txt") # Read N2P2 file
structures = reorder(structures) #Reorder atoms and structures for VASP (optional)
write_MLAB("example_MLAB",structures) # Write VASP ML_AB file
write_XYZ("example_XYZ",structures) # Write ase XYZ file
```
### Available functions:
- `read_N2P2()`: reads N2P2 input files
- `read_MLAB()`: reads VASP ML_AB input files
- `read_XYZ()`: reads ase extended XYZ files
- `write_N2P2()`: writes N2P2 input files
- `write_MLAB()`: writes VASP ML_AB input files
- `write_XYZ()`: writes ase extended XYZ files
- `reorder()`: reorders atoms within a conformation as well as the sequence of conformations for VASP. This is optional, and reordering can also be done by VASP later
