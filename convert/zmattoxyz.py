import numpy as np

def gmetry(natoms, geo, na, nb, nc):
    fail = False
    coord = np.zeros((3, natoms), dtype=np.float64)

    coord[:, 0] = [0.0, 0.0, 0.0]
    coord[0, 1] = geo[0, 1]
    coord[1, 1] = 0.0
    coord[2, 1] = 0.0

    if natoms == 2:
        return coord.T, fail

    ccos = np.cos(geo[1, 2])
    if na[2] == 0:
        coord[0, 2] = coord[0, 0] + geo[0, 2] * ccos
    else:
        coord[0, 2] = coord[0, 1] - geo[0, 2] * ccos
    coord[1, 2] = geo[0, 2] * np.sin(geo[1, 2])
    coord[2, 2] = 0.0

    for i in range(3, natoms):
        cosa = np.cos(geo[1, i])
        mb = nb[i]
        mc = na[i]

        xb = coord[0, mb] - coord[0, mc]
        yb = coord[1, mb] - coord[1, mc]
        zb = coord[2, mb] - coord[2, mc]

        rbc = 1.0 / np.sqrt(xb * xb + yb * yb + zb * zb)

        if abs(cosa) >= 0.9999999999:
            rbc2 = geo[0, i] * rbc * cosa
            coord[0, i] = coord[0, mc] + xb * rbc2
            coord[1, i] = coord[1, mc] + yb * rbc2
            coord[2, i] = coord[2, mc] + zb * rbc2
            continue

        ma = nc[i]

        xa = coord[0, ma] - coord[0, mc]
        ya = coord[1, ma] - coord[1, mc]
        za = coord[2, ma] - coord[2, mc]

        xyb = np.sqrt(xb * xb + yb * yb)
        k = -1
        if xyb <= 0.1:
            xpa = za
            za = -xa
            xa = xpa
            xpb = zb
            zb = -xb
            xb = xpb
            xyb = np.sqrt(xb * xb + yb * yb)
            k = 1

        costh = xb / xyb
        sinth = yb / xyb

        xpa = xa * costh + ya * sinth
        ypa = ya * costh - xa * sinth

        sinph = zb * rbc
        cosph = np.sqrt(abs(1.0 - sinph * sinph))

        xqa = xpa * cosph + za * sinph
        zqa = za * cosph - xpa * sinph

        yza = np.sqrt(ypa * ypa + zqa * zqa)
        if yza < 1e-4:
            fail = True
            return coord.T, fail

        coskh = ypa / yza
        sinkh = zqa / yza

        sina = np.sin(geo[1, i])
        sind = -np.sin(geo[2, i])
        cosd = np.cos(geo[2, i])

        xd = geo[0, i] * cosa
        yd = geo[0, i] * sina * cosd
        zd = geo[0, i] * sina * sind

        ypd = yd * coskh - zd * sinkh
        zpd = zd * coskh + yd * sinkh
        xpd = xd * cosph - zpd * sinph
        zqd = zpd * cosph + xd * sinph
        xqd = xpd * costh - ypd * sinth
        yqd = ypd * costh + xpd * sinth

        if k < 1:
            coord[0, i] = xqd + coord[0, mc]
            coord[1, i] = yqd + coord[1, mc]
            coord[2, i] = zqd + coord[2, mc]
        else:
            coord[0, i] = -zqd + coord[0, mc]
            coord[1, i] = yqd + coord[1, mc]
            coord[2, i] = xqd + coord[2, mc]

    return coord.T, fail


def read_zmatrix(filename):
    import numpy as np
    symbols = []
    na, nb, nc = [], [], []
    bond_len, bond_ang, dihed_ang = [], [], []

    with open(filename) as f:
        lines = [line.strip() for line in f if line.strip()]
    for i, line in enumerate(lines):
        parts = line.split()
        symbols.append(parts[0])
        if i == 0:
            na.append(0); nb.append(0); nc.append(0)
            bond_len.append(0.0); bond_ang.append(0.0); dihed_ang.append(0.0)
        elif i == 1:
            na.append(int(parts[1]))
            bond_len.append(float(parts[2]))
            nb.append(0)
            nc.append(0)
            bond_ang.append(0.0)
            dihed_ang.append(0.0)
        elif i == 2:
            na.append(int(parts[1]))
            bond_len.append(float(parts[2]))
            nb.append(int(parts[3]))
            bond_ang.append(np.deg2rad(float(parts[4])))
            nc.append(0)
            dihed_ang.append(0.0)
        else:
            na.append(int(parts[1]))
            bond_len.append(float(parts[2]))
            nb.append(int(parts[3]))
            bond_ang.append(np.deg2rad(float(parts[4])))
            nc.append(int(parts[5]))
            dihed_ang.append(np.deg2rad(float(parts[6])))

    # Umwandlung in 0-basierte Indizes für Python
    na = [x - 1 if x > 0 else 0 for x in na]
    nb = [x - 1 if x > 0 else 0 for x in nb]
    nc = [x - 1 if x > 0 else 0 for x in nc]

    return symbols, na, nb, nc, np.array(bond_len), np.array(bond_ang), np.array(dihed_ang)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python zmat_to_xyz_xtb.py zmatrix.zmat")
        sys.exit(1)

    fname = sys.argv[1]
    symbols, na, nb, nc, bond_len, bond_ang, dihed_ang = read_zmatrix(fname)
    xyz, fail = gmetry(len(symbols), np.vstack([bond_len, bond_ang, dihed_ang]), na, nb, nc)

    if fail:
        print("Fehler bei Rücktransformation")
    else:
        outname = "expected.xyz"
        with open(outname, "w") as f:
            f.write(f"{len(symbols)}\n")
            f.write("Reconstructed from Z-Matrix using xtb gmetry\n")
            for s, coord in zip(symbols, xyz):
                f.write(f"{s} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
        print(f"Wrote reconstructed XYZ to {outname}")
