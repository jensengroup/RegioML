import itertools
from rdkit import Chem
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


class ChargeAutocorrelation(DescriptorElement):
    description = 'Autocorrelationfunction for partial charges on grid'

    def __init__(self, **options):
        self.charge_type = options.pop('charge_type')
        self.rmin = options.pop('rmin', 1.0)
        self.rmax = options.pop('rmax', 11.0)
        self.no_h = options.pop('no_h', False)
        self.step_size = options.pop('step_size', 0.2)
        self.split_pos_neg = options.pop('split_pos_neg', True)
        self.variant2 = options.pop('variant2', False)
        self.range = self.rmax - self.rmin
        assert (self.range/self.step_size).is_integer()
        self.range_values = [self.step_size * step for step in range(int(self.range/self.step_size))]
        super(ChargeAutocorrelation, self).__init__(**options)

    def calculate_elements(self, atom):
        """ Overridden base method that is interace to calculate and return the descriptor's x-vector
        :param atom: RDKit atom for which the descriptor vector should be calculated
        :return: x-vector of descriptor
        :rtype: list
        """
        atom_idx = atom.GetIdx()
        parent_mol = atom.GetOwningMol()
        if self.no_h:
            hs_before_atom = get_num_h_before_atom(atom)
            atom_idx = atom_idx - hs_before_atom
            parent_mol = Chem.RemoveHs(parent_mol)
        dm = DistanceMatrix(parent_mol)
        all_atoms = parent_mol.GetAtoms()
        if self.split_pos_neg:
            ac_vector_p = []
            ac_vector_m = []
            for val in self.range_values:
                if not self.variant2:
                    ac_result = self._auto_correlation(self.rmin + val, self.rmin + val + self.step_size, atom_idx,
                                                               all_atoms, dm)
                    ac_vector_p.append(ac_result[0])
                    ac_vector_m.append(ac_result[1])
                else:
                    ac_result = self._auto_correlation_variant(self.rmin + val, self.rmin + val + self.step_size,
                                                               atom_idx, all_atoms, dm)
                    ac_vector_p.append(ac_result[0])
                    ac_vector_m.append(ac_result[1])
            return ac_vector_p + ac_vector_m
        else:
            ac_vector = []
            for val in self.range_values:
                if not self.variant2:
                    ac_vector.append(self._auto_correlation2(self.rmin + val, self.rmin +val + 1.5, atom_idx, all_atoms,
                                                             dm))
                else:
                    ac_vector.append(self._auto_correlation2_variant(self.rmin + val, self.rmin +val + 1.5, atom_idx,
                                                                     all_atoms, dm))
            return ac_vector

    def _auto_correlation(self, r_min, r_max, atom_idx, all_atoms, dm):
        """ Autocorrelation function with separating negative from positive contributions
        :param r_min: lower bound distance of bin
        :param r_max: upper bound distance of bin
        :param atom_idx: Index of atom for which the autocorrelation if calculated
        :param all_atoms: list of all atoms in the molecule
        :param dm distance matrix (cartesian space) of molecule
        :return: descriptor vector
        :rtype: list
        """
        atoms_in_region = []
        accum_plus = 0.0
        accum_minus = 0.0
        distances = dm.distance_matrix[:, atom_idx]
        for idx, distance in enumerate(distances):
            if r_min < distance <= r_max:
                atoms_in_region.append(all_atoms[idx])
        for i, j in itertools.combinations(atoms_in_region, 2):
            val = float(i.GetProp(self.charge_type)) * float(j.GetProp(self.charge_type))
            if val >= 0:
                accum_plus += val
            else:
                accum_minus += val
        return accum_plus, accum_minus

    def _auto_correlation_variant(self, r_min, r_max, atom_idx, all_atoms, dm):
        """ Autocorrelation function with separating negative from positive contributions
        :param r_min: lower bound distance of bin
        :param r_max: upper bound distance of bin
        :param atom_idx: index of atom for which the autocorrelation if calculated
        :param all_atoms: list of all atoms in the molecule
        :param dm distance matrix (cartesian space) of molecule
        :return: descriptor vector
        :rtype: list
        """
        atoms_in_region = []
        accum_plus = 0.0
        accum_minus = 0.0
        a_charge = float(all_atoms[atom_idx].GetProp(self.charge_type))
        distances = dm.distance_matrix[:, atom_idx]
        for idx, distance in enumerate(distances):
            if r_min < distance <= r_max:
                atoms_in_region.append(all_atoms[idx])
        for atom in atoms_in_region:
            val = a_charge * float(atom.GetProp(self.charge_type))
            if val >= 0:
                accum_plus += val
            else:
                accum_minus += val
        return accum_plus, accum_minus

    def _auto_correlation2(self, r_min, r_max, atom_idx, all_atoms, distance_matrix):
        """ Autocorrelation function without separating negative from positive contributions
        :param r_min: lower bound distance of bin
        :param r_max: upper bound distance of bin
        :param atom_idx: Index of atom for which the autocorrelation if calculated
        :param all_atoms: list of all atoms in the molecule
        :param distance_matrix: distance matrix (cartesian space) of molecule
        :return: descriptor vector
        :rtype: list
        """
        atoms_in_region = []
        accum = 0.0
        distances = distance_matrix.distance_matrix[:, atom_idx]
        for idx, distance in enumerate(distances):
            if r_min < distance < r_max:
                atoms_in_region.append(all_atoms[idx])
        for i, j in itertools.combinations(atoms_in_region, 2):
            # see : http://stackoverflow.com/questions/942543/operation-on-every-pair-of-element-in-a-list
            accum += float(i.GetProp(self.charge_type)) * float(j.GetProp(self.charge_type))
        return accum

    def _auto_correlation2_variant(self, r_min, r_max, atom_idx, all_atoms, dm):
        """ Autocorrelation function with separating negative from positive contributions
        :param r_min: lower bound distance of bin
        :param r_max: upper bound distance of bin
        :param atom:  atom for which the autocorrelation if calculated
        :param all_atoms: list of all atoms in the molecule
        :param dm distance matrix (cartesian space) of molecule
        :return: descriptor vector
        :rtype: list
        """
        atoms_in_region = []
        accum = 0.0
        a_charge = float(all_atoms[atom_idx].GetProp(self.charge_type))
        distances = dm.distance_matrix[:, atom_idx]
        for idx, distance in enumerate(distances):
            if r_min < distance <= r_max:
                atoms_in_region.append(all_atoms[idx])
        for atom in atoms_in_region:
            val = a_charge * float(atom.GetProp(self.charge_type))
            accum += val
        return accum