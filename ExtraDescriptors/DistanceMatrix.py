from numpy import zeros
from numpy import array

__author__ = 'modlab'


# TODO: Distance Matrix is crappy this way. It should be a matrix as such, i.e. child of numpy array with overriden
# TODO: methods to calculate the right elements or something. But it is strange to address an extra field when indexing
class DistanceMatrix(object):
    """ Class to store the pairwise interatomic distances of molecule
    """
    # yet specific to RDkit structures, I should keep that in mind ...
    # Could be speeded up by factor two by using symmetry
    def __init__(self, mol):
        self.mol = mol
        self.distance_matrix = self._calculate_dm()

    def _calculate_dm(self):
        atom_coordinates = [self.mol.GetConformer().GetAtomPosition(i) for i in range(len(self.mol.GetAtoms()))]
        num_atoms = len(atom_coordinates)
        dm = zeros((num_atoms, num_atoms))
        for i, ac_i in enumerate(atom_coordinates):
            for j, ac_j in enumerate(atom_coordinates):
                dm[i,j] = ac_i.Distance(ac_j)
        return dm


    # Other possibility from RDKIT: Chem.Get3DDistanceMatrix(mol)