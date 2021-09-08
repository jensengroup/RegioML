import itertools
from rdkit.Chem import rdmolops
from DescriptorElement import DescriptorElement
from rdkit import Chem

__author__ = 'modlab'


class ChargeGraphAutocorrelation(DescriptorElement):
    description = 'Autocorrelationfunction for partial charges based on molecular graph'

    def __init__(self, **options):
        self.charge_type = options.pop('charge_type')
        self.rmin = options.pop('rmin', 1)
        self.rmax = options.pop('rmax', 10)
        self.no_h = options.pop('no_h', False)
        self.range_values = range(self.rmin, self.rmax+1)
        super(ChargeGraphAutocorrelation, self).__init__(**options)

    def calculate_elements(self, atom):
        parent_mol = atom.GetOwningMol()
        if self.no_h:
            parent_mol = Chem.RemoveHs(parent_mol)
        tdm = rdmolops.GetDistanceMatrix(parent_mol)  # topological distance matrix
        all_atoms = parent_mol.GetAtoms()
        ac_vector_p = []
        ac_vector_m = []
        for val in self.range_values:
            # ac_vector.append(self._auto_correlation(self.rmin + val, self.rmin +val + 1.5, atom, all_atoms, dm))
            ac_vector_p.append(self._auto_correlation(val, atom, all_atoms, tdm)[0])
            ac_vector_m.append(self._auto_correlation(val, atom, all_atoms, tdm)[1])
        return ac_vector_p + ac_vector_m

    def _auto_correlation(self, d, atom, all_atoms, tdm):
        atoms_in_region = []
        accum_plus = 0.0
        accum_minus = 0.0
        a_idx = atom.GetIdx()
        distances = tdm[:, a_idx]
        for idx, distance in enumerate(distances):
            if distance == d:
                atoms_in_region.append(all_atoms[idx])
        for i, j in itertools.combinations(atoms_in_region, 2):
            val = float(i.GetProp(self.charge_type)) * float(j.GetProp(self.charge_type))
            if val >= 0:
                accum_plus += val
            else:
                accum_minus += val
        return [accum_plus, accum_minus]