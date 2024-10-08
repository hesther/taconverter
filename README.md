taconverter
==============================

Converter for VASP ML, N2P2 and ase XYZ files

### Use

```
from taconverter import *
structures = read_N2P2("example_N2P2.txt") # Read N2P2 file
structures = reorder(structures) #Reorder atoms and structures for VASP (optional)
write_MLAB("example_MLAB",structures) # Write VASP ML_AB file
write_XYZ("example_XYZ",structures) # Write ase XYZ file
```

