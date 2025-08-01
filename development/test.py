from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

smiles = "CCO" # Ethanol
smiles = "C1CCCCC1" # Cyclohexan
smiles = "CCCO" # 1 propanol
smiles = "CC(O)C" # 2 propanol
smiles = "CCC" # propan
smiles = "C"

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

AllChem.EmbedMolecule(mol)
AllChem.UFFOptimizeMolecule(mol)

conf = mol.GetConformer()
atoms = []
for atom in mol.GetAtoms():
    idx = atom.GetIdx()
    pos = conf.GetAtomPosition(idx)
    atoms.append({
        "index": idx,
        "element": atom.GetSymbol(),
        "x": pos.x,
        "y": pos.y,
        "z": pos.z,
    })
[print(a) for a in atoms]

# Extract edges
bonds = []
for bond in mol.GetBonds():
    bonds.append({
        "begin": bond.GetBeginAtomIdx(),
        "end": bond.GetEndAtomIdx()
    })

bonds = [list(i.values()) for i in bonds]

# bond length 0.5 - 2.5 Angström
# bond angles 0°-180°
# Dihedral angles (torsion) -180° - 180°
# For each edge --> 1 bond length
# Wherever two bonds at one atom --> 1 bond angle
# Wherever 4 atoms are connected linearly (3 bonds) --> 1 dihedral

# Create graph
class Node:
    def __init__(self, id):
        self.id = id
        self.connected_to = []
        self.virtual = False
    @classmethod
    def virtual(cls, id):
        node = cls(id)
        node.virtual = True
        return node
    def add_connection(self, to):
        self.connected_to.append(to)
    def __repr__(self):
        if self.virtual:
            return "Virtual Node"
        return "Node"
    def __str__(self):
        id_str = "ID: "
        if self.virtual:
            id_str = "Virtual ID: "
        return id_str + str(self.id) + "; To: " + str(self.connected_to)


class Graph:
    def __init__(self):
        self.nodes = []
        self.id = -1
    def add_node(self, id):
        if self.id < id:
            self.nodes.append(Node(id))
            self.id = self.id + 1
    def add_connection(self, id, to):
        self.nodes[id].add_connection(to)
    def __repr__(self):
        return "Graph"
    def __str__(self):
        [print(n) for n in self.nodes]
        return ""
    def __len__(self):
        return len(self.nodes)

edges = bonds

def create_graph(edges):
    g = Graph()
    for e in edges:
        g.add_node(e[0])
        g.add_node(e[1])
    for e in edges:
        g.add_connection(e[0], e[1])
        g.add_connection(e[1], e[0])
    return(g)

graph = create_graph(edges)
print(graph)

# Retrieve parameters
def calc_lengths(edges):
    return len(edges), edges
bond_lengths = calc_lengths(edges)
print("Lengths:", bond_lengths)

def calc_angles(graph):
    count = 0
    angles = []
    for n in graph.nodes:
        connections = n.connected_to
        if len(connections) > 1:
            for i in range(len(connections)):
                for j in range(i + 1, len(connections)):
                    a = connections[i]
                    c = connections[j]
                    angles.append([a, n.id, c])
                    count += 1
    return count, angles

angles = calc_angles(graph)
print("Angles: ", angles)

def calc_dihedrals(graph):
    count = 0
    dihedrals = []
    seen = set()
    for b in graph.nodes:
        for c_id in b.connected_to:
            c = graph.nodes[c_id]
            key = tuple(sorted((b.id, c.id)))
            if key in seen:
                continue
            seen.add(key)
            for a_id in b.connected_to:
                if a_id == c.id:
                    continue
                for d_id in c.connected_to:
                    if d_id == b.id:
                        continue
                    dihedrals.append([a_id, b.id, c.id, d_id])
                    count += 1
    return count, dihedrals

dihedrals = calc_dihedrals(graph)
print("Dihedrals: ", dihedrals)

def needs_virtual(path):
    return len(path) < 3

def insert_virtual(id, graph, path):
    virtual_id = max(n.id for n in graph.nodes) + 1
    graph.nodes.append(Node.virtual(virtual_id))
    graph.nodes[virtual_id].add_connection(id)
    graph.nodes[id].add_connection(virtual_id)
    return path[-2:] + [virtual_id]

class Order:
    def __init__(self):
        self.order = []
    def add(self, elem):
        self.order.append(elem)
    def __len__(self):
        return len(self.order)

def traverse(graph, root = 0):
    visited = set()
    order = Order()
    counter = 0
    def dfs(current_id, path):
        if current_id in visited:
            return
        visited.add(current_id)
        if len(order) == 0:
            order.add((current_id, []))
        elif len(order) == 1:
            order.add((current_id, [path[0]]))
        elif len(order) > 1:
            if needs_virtual(path):
                path = insert_virtual(current_id, graph, path)
            order.add((current_id, path[-3:]))
        neighbors = list(graph.nodes[current_id].connected_to
        for neighbor_id in neighbors:
            if neighbor_id not in visited:
                dfs(neighbor_id, path + [current_id])
    dfs(root, [])
    return order

order = traverse(graph, 1)
[print(o) for o in order.order]
print(graph)

# Why standard SMILES → 3D tools don't help
# Tools like RDKit or Open Babel:
#   - Use random seeds or force-field optimization
#   - May yield different geometries for same SMILES
#   - Do not guarantee correspondence between parameters and resulting geometry
#   - Are not invertible: XYZ → parameter is not well-defined

# Backcalculation to xyz
# Atom 0 at origin: 0,0,0 or better choose atom with most connections
# Atom 1 placed along x-axis at given bond length
# Atom 2 placed in xy plane using bond length & bond angle
# Atom 3+ placed using bond length, and angle
# Afterwards apply dihedrals to construct the z coordinates

# Build placement tree:
# traverse graph:
# for each new atom record:
# its parents for bond length
# its grandparents for angle
# its great-grandparents for dihedral

class Atom:
    def __init__(self, id):
        self.id = id
        self.length = None
        self.angle = None
        self.dihedral = None
        self.x = Na
        self.y = Na
        self.z = Na
    def set_length(self, length):
        self.length = length
    def set_angle(self, angle):
        self.angle = angle
    def set_dihedral(self, dihedral):
        self.dihedral = dihedral

# Create atom instance for each entry in traversal order
# find for each instance the length, angle, and dihedral from the input parameters
# traverse and reconstruct xyz

def find_length_idx(l, r, bond_lengths):
    for i in range(len(bond_lengths[1])):
        length = bond_lengths[1][i]
        if (length[0] == l and length[1] == r) or (length[0] == r and length[1] == l):
            return(i)
    raise Exception("Couldn't find a bond length")

def within(elem, l):
    for i in l:
        if elem == i:
            return(True)
    return(False)

def which(elem, l):
    for i in range(len(l)):
        if elem == l[i]:
            return(i)
    return(Na)

def find_angle_idx(a, b, c, angles):
    for i in range(len(angles[1])):
        current_angle = angles[1][i]
        if within(a, current_angle) and within(b, current_angle) and within(c, current_angle):
            a_idx = which(a, current_angle)
            b_idx = which(b, current_angle)
            c_idx = which(c, current_angle)
            if a_idx != b_idx and a_idx != c_idx and b_idx != c_idx:
                return(i)
    raise Exception("Couldn't find an suitable angle")

def translate(parameters, bond_lengths, angles, dihedrals, order, n_atoms):
    n_lengths = bond_lengths[0]
    n_angles = angles[0]
    n_dihedrals = dihedrals[0]

    length_parameters = parameters[0:n_lengths]
    angle_parameters = parameters[n_lengths:n_lengths + n_angles]
    dihedral_parameters = parameters[n_lengths + n_angles:n_lengths + n_angles + n_dihedrals]

    atoms = []
    for i in range(len(order)):
        if i == 0:
            a = Atom(order[i][0])
            a.x = 0.0
            a.y = 0.0
            a.z = 0.0
            atoms.append(a)
        if i == 1:
            ref_atoms = order[i][1]
            idx = find_length_idx(order[i][0], ref_atoms[0], bond_lengths) # assume 2 entries in order[i]
            res = Atom(order[i][0])
            res.x = length_parameters[idx]
            res.y = 0.0
            res.z = 0.0
            atoms.append(res)
        if i == 2:
            a = order[i][0]
            ref_atoms = order[i][1]
            b = ref_atoms[len(ref_atoms) - 1]
            c = ref_atoms[len(ref_atoms) - 2]
            idx_length = find_length_idx(a, b, bond_lengths)
            idx_angle = find_angle_idx(a, b, c, angles)
            r = length_parameters[idx_length]
            theta = angle_parameters[idx_angle]
            res = Atom(a)
            atom_b = atoms[which(b, [i.id for i in atoms])]
            atom_c = atoms[which(c, [i.id for i in atoms])]
            res.x = atom_b.x + r * cos(theta)
            res.y = atom_b.y + r * sin(theta)
            res.z = 0.0
            atoms.append(res)








