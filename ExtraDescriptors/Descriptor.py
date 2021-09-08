__author__ = 'modlab'


class Descriptor(object):
    def __init__(self):
        self.header = None
        self.descriptor_vectors = []
        self.x = []
        self.y = None

    def __repr__(self):
        return '{} ::: {} {}'.format(str(self.header), str(self.x), str(self.y))


class AtomDescriptor(Descriptor):
    def __init__(self, mol_id, atom_idx):
        super(AtomDescriptor, self).__init__()
        self.header = dict(mol_id=mol_id, atom_idx=atom_idx)


class MolDescriptor(Descriptor):
    def __init__(self, mol_id):
        super(MolDescriptor, self).__init__()
        self.header = dict(mol_id=mol_id)