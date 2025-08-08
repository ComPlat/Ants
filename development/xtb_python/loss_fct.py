import numpy as np
import random
from pso import pso
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tblite.interface import Calculator

class Atom:
    def __init__(self, atom_name, x, y, z):
        assert isinstance(atom_name, str), "Name has to be a string"
        assert len(atom_name) <= 2, "The atom name has a wrong length (<= 2 is expected)"
        assert isinstance(x, float), "X coordinate is not a real number"
        assert isinstance(y, float), "Y coordinate is not a real number"
        assert isinstance(z, float), "Z coordinate is not a real number"
        self.atom_name = atom_name
        self.x = x
        self.y = y
        self.z = z
    def stringify_for_file(self):
        return self.atom_name + '\t' + str(self.x) + '\t' + str(self.y) + '\t' + str(self.z)
    def stringify(self):
        return self.atom_name + ' ' + str(self.x) + ' ' + str(self.y) + ' ' + str(self.z)
    def set_coordinates(self, coordinates):
        self.x = coordinates[0]
        self.y = coordinates[1]
        self.z = coordinates[2]

class Molecule:
    def __init__(self):
        self.atoms = []
    def add(self, atom):
        assert isinstance(atom, Atom), "Expected an argument of type Atom"
        self.atoms.append(atom)
    def create_file(self, file_name):
        table = '\n'.join(atom.stringify_for_file() for atom in self.atoms)
        content = str(self.get_num_parameter() / 3) + "\n\n" + table
        with open(file_name, "w") as text_file:
            text_file.write(content)
    def stringify(self):
        return ';'.join(atom.stringify() for atom in self.atoms)
    def get_num_parameter(self):
        return len(self.atoms) * 3
    def set_coordinates(self, vector):
        assert isinstance(vector, np.ndarray), "Expected an argument of type numpy array"
        assert len(vector) % 3 == 0, "% 3 != 0"
        groups = vector.reshape(-1, 3)
        assert len(groups) == len(self.atoms)
        for i in range(len(groups)):
            self.atoms[i].set_coordinates(groups[i])

def evaluate_energy(vector, atomic_numbers):
    try:
        coordinates = vector.reshape(-1, 3)
        calc = Calculator(method="GFN2-xTB", numbers=atomic_numbers, positions=coordinates)
        result = calc.singlepoint()
        return(result.get("energy"))
    except:
        return float("inf")

random.seed(1234)
npar = 9
atomic_numbers = [8, 1, 1]
res = pso(npar, atomic_numbers, -0.9, 0.9, 10, 10, -100000, evaluate_energy)
print(res)
