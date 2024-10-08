from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase import io
from ase.build.tools import sort

def read_MLAB(filename):
    structures = []
    with open(filename) as f:
        lines = [line.rstrip().lstrip() for line in f]

    conf_i = [i for i, line in enumerate(lines) if line == "The number of atoms"]
    conf_i.append(len(lines)+1)

    for i, ii in enumerate(conf_i[:-1]):
        sublines = lines[ii+3:conf_i[i+1]]

        #Num atoms
        num = int(lines[ii+2])

        #Energy
        si = sublines.index("Total energy (eV)") + 2
        energy = float(sublines[si])

        #Symbols
        symbols = []
        si = sublines.index("Atom types and atom numbers") + 2
        new_line = sublines[si]
        while not "===" in new_line:
            symbol, times = new_line.split()
            symbols.extend([symbol] * int(times))
            si += 1
            new_line = sublines[si]

        #Positions
        positions = []
        si = sublines.index("Atomic positions (ang.)") + 2
        new_line = sublines[si]
        while not "===" in new_line:
            positions.append([float(x) for x in new_line.split()])
            si += 1
            new_line = sublines[si]

        #Forces
        forces = []
        si = sublines.index("Forces (eV ang.^-1)") + 2
        new_line = sublines[si]
        while not "===" in new_line:
            forces.append([float(x) for x in new_line.split()])
            si += 1
            new_line = sublines[si]

        #Cell
        cell = []
        si = sublines.index("Primitive lattice vectors (ang.)") + 2
        new_line = sublines[si]
        while not "===" in new_line:
            cell.append([float(x) for x in new_line.split()])
            si += 1
            new_line = sublines[si]

        #Stress
        stress = []
        si = sublines.index("Stress (kbar)") + 2
        stress.extend([float(x) for x in sublines[si+2].split()])
        stress.extend([float(x) for x in sublines[si+6].split()])

        #Make ase object
        atoms = Atoms(symbols, positions = positions)
        if cell != []:
            atoms.set_cell(cell)
            atoms.set_pbc(True)

        atoms.calc = SinglePointCalculator(
            atoms,
            forces = forces,
            energy = energy,
            stress = stress
        )
        
        structures.append(atoms)

    return structures

def read_XYZ(filename):
    structures = [structure for structure in io.iread(filename)]
    return structures

def read_N2P2(filename):
    structures = []
    with open(filename) as f:
        lines = [line.rstrip().lstrip() for line in f]
    conf_i = [i for i, line in enumerate(lines) if line == "begin"]
    conf_i.append(len(lines)+1)

    for i, ii in enumerate(conf_i[:-1]):
        sublines = lines[ii+2:conf_i[i+1]]
        

        #Cell
        cell = []
        new_lines = [l for l in sublines if 'lattice' in l]
        for new_line in new_lines:
            cell.append([float(x) for x in new_line.split()[1:]])

        #Symbols, positions and forces
        symbols = []
        positions = []
        forces = []
        new_lines = [l for l in sublines if 'atom' in l]
        for new_line in new_lines:
            values = new_line.split()
            symbols.append(values[4])
            positions.append([float(x) for x in values[1:4]])
            forces.append([float(x) for x in values[7:10]])

        #Energy
        energy = float([l for l in sublines if 'energy' in l][0].split()[1])

        #Make ase object
        atoms = Atoms(symbols, positions = positions)
        if cell != []:
            atoms.set_cell(cell)
            atoms.set_pbc(True)

        atoms.calc = SinglePointCalculator(
            atoms,
            forces = forces,
            energy = energy)
        
        structures.append(atoms)

    return structures
            
def write_MLAB(filename, structures):
    text = ""
    text += " 1.0 Version\n"
    text += "**************************************************\n"
    text += "     The number of configurations\n"
    text += "--------------------------------------------------\n"
    num_confs = len(structures)
    text += "        " + str(num_confs) + "\n"
    text += "**************************************************\n"
    text += "     The maximum number of atom type\n"
    text += "--------------------------------------------------\n"
    max_num_atom_type = max([len(structure.symbols.species()) for structure in structures])
    text += "       " + str(max_num_atom_type) + "\n"
    text += "**************************************************\n"
    text += "     The atom types in the data file\n"
    text += "--------------------------------------------------\n"
    text += "     "
    #TODO: sort types according to actual order in conformations
    atom_types = sorted(list(set([item for structure in structures for item in list(structure.symbols.species())])), key = lambda x: Atoms(x).symbols.numbers[0], reverse=True)
    ctr = 0
    for t in atom_types:
        ctr += 1
        text += str(t)
        if ctr == len(atom_types):
            text += "\n"
        elif ctr == 3:
            text += "\n     "
        else:
            text += " "
    text += "**************************************************\n"
    text += "     The maximum number of atoms per system\n"
    text += "--------------------------------------------------\n"
    max_num_atoms_system = max([len(structure) for structure in structures])
    text += "             " + str(max_num_atoms_system) + "\n"
    text += "**************************************************\n"
    text += "     The maximum number of atoms per atom type\n"
    text += "--------------------------------------------------\n"
    max_num_atoms_type = max([max([structure.symbols.count(s) for s in atom_types]) for structure in structures]) 
    text += "             " + str(max_num_atoms_type) + "\n"
    text += "**************************************************\n"
    text += "     Reference atomic energy (eV)\n"
    text += "--------------------------------------------------\n"
    text += "     "
    ctr = 0
    for t in atom_types:
        ctr += 1
        text += str("0.0")
        if ctr == len(atom_types):
            text += "\n"
        elif ctr == 3:
            text += "\n     "
        else:
            text += " "
    text += "**************************************************\n"
    text += "     Atomic mass\n"
    text += "--------------------------------------------------\n"
    text += "     "
    ctr = 0
    for t in atom_types:
        ctr += 1
        text += str(Atoms(t).get_masses()[0])
        if ctr == len(atom_types):
            text += "\n"
        elif ctr == 3:
            text += "\n     "
        else:
            text += " "
    text += "**************************************************\n"
    text += "     The numbers of basis sets per atom type\n"
    text += "--------------------------------------------------\n"
    text += "     "
    ctr = 0
    for t in atom_types:
        ctr += 1
        text += str("1")
        if ctr == len(atom_types):
            text += "\n"
        elif ctr == 3:
            text += "\n     "
        else:
            text += " "
    for t in atom_types:
        text += "**************************************************\n"
        text += "     Basis set for " + t + "\n"
        text += "--------------------------------------------------\n"
        text += "          1      1\n"

    ctr = 1
    for structure in structures:
        text += "**************************************************\n"
        text += "     Configuration num.      " + str(ctr) + "\n"
        text += "==================================================\n"
        text += "     System name\n"
        text += "--------------------------------------------------\n"
        text += "     Name\n"
        text += "==================================================\n"
        text += "     The number of atom types\n"
        text += "--------------------------------------------------\n"
        num_atom_type = len(structure.symbols.species())
        text += "       " + str(num_atom_type) + "\n"
        text += "==================================================\n"
        text += "     The number of atoms\n"
        text += "--------------------------------------------------\n"
        num_atoms = len(structure)
        text += "         " + str(num_atoms) + "\n"
        text += "**************************************************\n"
        text += "     Atom types and atom numbers\n"
        text += "--------------------------------------------------\n"
        type_and_num = {}
        for s in list(structure.symbols):
            if s not in type_and_num.keys():
                type_and_num[s] = 1
            else:
                type_and_num[s] += 1
        for t in type_and_num.keys():
            text += "     " + t +"      " + str(type_and_num[t]) + "\n"
        text += "==================================================\n"
        text += "     CTIFOR\n"
        text += "--------------------------------------------------\n"
        text += "   0.0\n"
        text += "=================================================\n"
        text += "     Primitive lattice vectors (ang.)\n"
        text += "--------------------------------------------------\n"
        for i in range(3):
            text += " " + " ".join(["{0:0.15E}".format(i) for i in structure.cell[:][i]]) + "\n"
        text += "==================================================\n"
        text += "     Atomic positions (ang.)\n"
        text += "--------------------------------------------------\n"
        for i in range(num_atoms):
            text += " " + " ".join(["{0:0.15E}".format(i) for i in structure.positions[i]]) + "\n"
        text += "==================================================\n"
        text += "     Total energy (eV)\n"
        text += "--------------------------------------------------\n"
        text += "  " + str(structure.get_potential_energy()) + "\n"
        text += "==================================================\n"
        text += "     Forces (eV ang.^-1)\n"
        text += "--------------------------------------------------\n"
        for i in range(num_atoms):
            text += " " + " ".join(["{0:0.15E}".format(i) for i in structure.get_forces()[i]]) + "\n"
        text += "==================================================\n"
        text += "     Stress (kbar)\n"
        text += "--------------------------------------------------\n"
        text += "     XX YY ZZ\n"
        text += "--------------------------------------------------\n"
        try:
            stress = structure.get_stress()
        except:
            stress = [0.0]*6
        text += " " + " ".join(["{0:0.15E}".format(i) for i in stress[:3]]) + "\n"
        text += "--------------------------------------------------\n"
        text += "     XY YZ ZX\n"
        text += "--------------------------------------------------\n"
        text += " " + " ".join(["{0:0.15E}".format(i) for i in stress[3:]]) + "\n"
        ctr += 1

    with open(filename, 'w') as f:
        print(text, file=f)

def write_XYZ(filename, structures):
    io.write(filename, structures)

def write_N2P2(filename, structures):
    text = ""
    for structure in structures:
        text += "begin\n"
        text += "comment\n"
        num_atoms = len(structure)
        if True in structure.get_pbc():
            for i in range(3):
                text += "lattice " + " ".join(["{0:0.15E}".format(i) for i in structure.cell[:][i]]) + "\n"
        for i in range(num_atoms):
            text += "atom " + " ".join(["{0:0.15E}".format(i) for i in structure.positions[i]]) + " " + list(structure.symbols)[i] + " 0.0 0.0 " + " ".join(["{0:0.15E}".format(i) for i in structure.get_forces()[i]]) + "\n"
        text += "energy " + str(structure.get_potential_energy()) + "\n"
        text += "charge 0.0\n"
        text += "end\n"

    with open(filename, 'w') as f:
        print(text, file=f)

def reorder(structures):
    reordered = []
    
    #Reorder atoms within each structure
    atom_types = sorted(list(set([item for structure in structures for item in list(structure.symbols.species())])), key = lambda x: Atoms(x).symbols.numbers[0], reverse=True)
    for structure in structures:
        sort_key = 1/structure.symbols.numbers
        symbols = list(structure.symbols)
        positions = structure.positions
        cell = structure.cell[:]
        energy = structure.get_potential_energy()
        forces = structure.get_forces()

        symbols = [x for _, x in sorted(zip(sort_key, symbols), key=lambda pair: pair[0])]
        positions = [x for _, x in sorted(zip(sort_key, positions), key=lambda pair: pair[0])]
        forces = [x for _, x in sorted(zip(sort_key, forces), key=lambda pair: pair[0])]

        atoms = Atoms(symbols, positions = positions)
        if True in structure.pbc:
            atoms.set_cell(cell)
            atoms.set_pbc(True)
            
        atoms.calc = SinglePointCalculator(
            atoms,
            forces = forces,
            energy = energy)
        reordered.append(atoms)
        

    #Reorder structures
    reordered = sorted(reordered, key = lambda x: str(x.symbols))
    return reordered
