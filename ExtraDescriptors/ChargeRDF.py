from rdkit import Chem
import numpy as np
from DescriptorElement import DescriptorElement
from DistanceMatrix import DistanceMatrix

__author__ = 'modlab'

def get_num_h_before_atom(atom):
    mol = atom.GetOwningMol()
    a_idx = atom.GetIdx()
    h_count = 0
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx == a_idx:
            break
        if atom.GetAtomicNum() == 1:
            h_count += 1
    return h_count


class ChargeRDF(DescriptorElement):
    """
    Implementation of Formula according to formula (54) on page 1027 in Handbook of Cheminformatics Vol. 3
    """
    description = "Radial distribution function of query atom on grid"

    def __init__(self, **options):
        self.charge_type = options.pop('charge_type')
        self.beta = options.pop('beta', 20)
        self.sf = options.pop('scaling', 10)
        self.r_min = options.pop('rmin', 1)
        self.r_max = options.pop('rmax', 11)
        self.step_size = options.pop('step_size', 0.1)
        self.no_h = options.pop('no_h', False)
        self.r_values = [self.r_min + val * self.step_size for val in range(int((self.r_max - self.r_min) / self.step_size))]
        super(ChargeRDF, self).__init__(**options)

    def calculate_elements(self, atom):
        self._parent_mol = atom.GetOwningMol()
        self._a_idx = atom.GetIdx()
        if self.no_h:
            hs_before_atom = get_num_h_before_atom(atom)
            self._a_idx = self._a_idx - hs_before_atom
            self._parent_mol = Chem.RemoveHs(self._parent_mol)
        try:
            self._dm = self._parent_mol._dm
        except AttributeError:
            self._parent_mol._dm = DistanceMatrix(self._parent_mol)
        self._dm = self._parent_mol._dm
        self._all_atoms = self._parent_mol.GetAtoms()
        self._distances = self._dm.distance_matrix[:, self._a_idx]
        rdf = []
        for val in self.r_values:
            rdf.append(self._radial_distribution_function(val, self._a_idx, self._all_atoms))
        return rdf

    def _radial_distribution_function(self, R, atom_idx, all_atoms):
        acc = 0.0
        for atm in all_atoms:
            a_idx = atm.GetIdx()
            if a_idx is not self._a_idx:
                acc += self.sf * float(atm.GetProp(self.charge_type)) * float(all_atoms[atom_idx].GetProp(
                    self.charge_type)) * \
                    np.exp(-self.beta * np.power(R - self._dm.distance_matrix[a_idx, self._a_idx], 2.))
        return acc
