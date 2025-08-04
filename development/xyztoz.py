import numpy as np
from ase.data import covalent_radii
from ase.io import read

def get_bonded_atoms(atoms, tol=0.45):
    nat = len(atoms)
    bonded = [[] for _ in range(nat)]
    positions = atoms.get_positions()
    numbers = atoms.get_atomic_numbers()

    for i in range(nat):
        for j in range(i+1, nat):
            dist = np.linalg.norm(positions[i] - positions[j])
            r_cov = covalent_radii[numbers[i]] + covalent_radii[numbers[j]] + tol
            if dist <= r_cov:
                bonded[i].append(j)
                bonded[j].append(i)
    return bonded

def angle_between(a, b, c):
    ba = a - b
    bc = c - b
    cos_angle = np.dot(ba, bc) / (np.linalg.norm(ba)*np.linalg.norm(bc))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    return np.arccos(cos_angle)

def dihedral(p0, p1, p2, p3):
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)

def build_zmatrix(atoms):
    nat = len(atoms)
    bonded = get_bonded_atoms(atoms)
    positions = atoms.get_positions()

    na = [None]*nat
    nb = [None]*nat
    nc = [None]*nat

    na[0] = None
    nb[0] = None
    nc[0] = None

    for i in range(1, nat):
        candidates = [x for x in bonded[i] if x < i]
        if not candidates:
            candidates = [j for j in range(i) if j != i]
        na[i] = min(candidates)

        if i == 1:
            nb[i] = None
            nc[i] = None
            continue

        candidates_nb = [x for x in bonded[na[i]] if x < i and x != i]
        if not candidates_nb:
            candidates_nb = [j for j in range(i) if j != i and j != na[i]]
        nb[i] = min(candidates_nb)

        if i == 2:
            nc[i] = None
            continue

        candidates_nc = [x for x in bonded[nb[i]] if x < i and x != na[i] and x != i]
        if not candidates_nc:
            candidates_nc = [j for j in range(i) if j not in (i, na[i], nb[i])]
        nc[i] = min(candidates_nc) if candidates_nc else None

    geo = np.zeros((nat, 3))
    for i in range(nat):
        if na[i] is None:
            geo[i, :] = 0.0
        else:
            geo[i, 0] = np.linalg.norm(positions[i] - positions[na[i]])
            if nb[i] is None:
                geo[i, 1] = 0.0
                geo[i, 2] = 0.0
            else:
                geo[i, 1] = angle_between(positions[i], positions[na[i]], positions[nb[i]])
                if nc[i] is None:
                    geo[i, 2] = 0.0
                else:
                    geo[i, 2] = dihedral(positions[nc[i]], positions[nb[i]], positions[na[i]], positions[i])

    return na, nb, nc, geo

def write_zmatrix(filename, atoms, na, nb, nc, geo):
    with open(filename, 'w') as f:
        nat = len(atoms)
        for i in range(nat):
            sym = atoms[i].symbol
            bl = geo[i, 0]
            ang = np.degrees(geo[i, 1])
            dih = np.degrees(geo[i, 2])
            if dih > 180.0:
                dih -= 360.0
            na_out = na[i] + 1 if na[i] is not None else 0
            nb_out = nb[i] + 1 if nb[i] is not None else 0
            nc_out = nc[i] + 1 if nc[i] is not None else 0

            if i == 0:
                f.write(f"{sym}\n")
            elif i == 1:
                f.write(f"{sym} {na_out} {bl:.6f}\n")
            elif i == 2:
                f.write(f"{sym} {na_out} {bl:.6f} {nb_out} {ang:.6f}\n")
            else:
                f.write(f"{sym} {na_out} {bl:.6f} {nb_out} {ang:.6f} {nc_out} {dih:.6f}\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python xyztozmat.py molecule.xyz")
        sys.exit(1)
    atoms = read(sys.argv[1], format='xyz')
    na, nb, nc, geo = build_zmatrix(atoms)
    outname = "zmatrix_" + sys.argv[1].split('.')[0] + ".zmat"
    write_zmatrix(outname, atoms, na, nb, nc, geo)
    print(f"Wrote Z-Matrix to {outname}")
