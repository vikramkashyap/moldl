#!/usr/bin/env python3

def convert_mol_to_xyz(infile):
    name = infile.split('.mol')[0]
    with open(infile, 'r') as file:
        lines = file.readlines()

    istart = 0
    for l in lines:
        if 'V2000' in l.split() or 'V3000' in l.split():
            break
        else:
            istart = istart + 1

    atomre = re.compile(r'\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\w+)\s+')
    atoms = []
    for l in lines[istart + 1:]:
        m = atomre.match(l)
        if m is not None:
            atoms.append(m.groups())
        else:
            break

    with open(name + '.xyz', 'w') as outfile:
        outfile.write(str(len(atoms)) + '\n')
        outfile.write('File made automatically with Python script using input mol file. (v1)\n')
        for a in atoms:
            outfile.write('{}{:>15}{:>15}{:>15}\n'.format(a[-1], a[0], a[1], a[2]))
