from numpy import array
from rdkit import Chem
from sklearn import preprocessing
import random

from ChargeAutocorrelation import ChargeAutocorrelation
from ChargeGraphAutocorrelation import ChargeGraphAutocorrelation
from ChargeRDF import ChargeRDF
from ChargeShell import ChargeShell
from ChargeShell import MassShell
from Descriptor import AtomDescriptor, MolDescriptor
from GraphChargeShell import GraphChargeShell
from PartialChargeDict import PartialChargeDict as PCD

__author__ = 'modlab'

"""Based on RDKit structures at the moment. On a later stage this should be replaced with my own molecule class"""


class DescriptorCreator(object):
    """
    Abstract base class for Atom and Molecule descriptor creators.
    Abstract descriptor (property) creator factory.
    Implements the common functionality.
    """
    def __init__(self, data_file_path):
        """
        :param data_file_path: path to file where the data for descriptor calculation is located
        :return:
        """
        self.data_file_path = data_file_path
        self.suppl = self._load_data(self.data_file_path)
        self.descriptor_properties = []
        self.descriptors = []
        # Descriptor creator class stores Matrix of all X  and Y numerical values and v
        self._X = array([])
        self._Y = array([])

    def reset(self):
        self.descriptor_properties = []
        self.descriptors = []
        # Descriptor creator class stores Matrix of all X  and Y numerical values and v
        self._X = array([])
        self._Y = array([])

    @property
    def X(self):
        # Test of size is necessary for numpy arrays. "if not self._X" does not work
        if self._X.size == 0:
            self._X = array([dv.x for dv in self.descriptors])
        return self._X

    @property
    def y(self):
        if self._Y.size == 0:
            self._Y = array([dv.y for dv in self.descriptors])
        return self._Y

    def scale_x(self):
        X_scaled = preprocessing.scale(self.X)
        for idx, dv in enumerate(self.descriptors):
            dv.x = list(X_scaled[idx])
        self._X = self._y = array([])

    def _load_data(self, data_file):
        """
        Hook that derived (child) classes have to implement to individualize data loading.
        Load the data from a data file and store in self.suppl. Generators might be used.
        :param data_file: File to load data from
        :return: list of data
        """
        raise NotImplementedError

    def create_property(self, prop_name, **options):
        """
        Creates properties and adds them to self.descriptor_properties, which then can be used by
        self.create_descriptor_vectors() to build the descriptor vectors.
        Property factory method.
        :param prop_name:
        :param options:
        :return: None
        """
        raise NotImplementedError

    def create_descriptor_vectors(self):
        """
        This function has to create the descriptor vectors and store them to self.descriptor.descriptors
        (header + descriptor)
        Note: the descriptor vectors should be set to zero before appending to them. Otherwise their
        state is messed up.
        :return: None
        """
        raise NotImplementedError


class ChemDescriptorCreator(DescriptorCreator):
    """ Abstract base class for molecular and atomic descriptors
        Implements the common methods for both classes.
    """
    def __init__(self, sdf_file_path):
        """
        :param sdf_file_path:
        :return:
        """
        super(ChemDescriptorCreator, self).__init__(sdf_file_path)

    def _load_data(self, sdf_file_path):
        """ Hook that can be overridden when an other molecule source than sd_file is required.
        Now the data is loaded in a list fully stored in memory, which is not optimal but necessary because
        the "generator-like" solution to directly work with RDKit resulted in crashes due to RDKit bugs
        :return: Iterable containing all structures
        """
        # TODO NOTE: crashes if stay with RDkit. Thus list has to be stored explicitly
        supplier = Chem.SDMolSupplier(sdf_file_path, removeHs=False)
        molecule_list = [mol for mol in supplier]
        return molecule_list

    def write_sdf(self, sd_file_path):
        """ Writes the molecues to SD file using RDkit
        :param sd_file_path: str path to sd file to write to
        :return: None
        """
        for mol in self.suppl:
            id = mol.GetProp('Mol_ID')
            for descriptor in self.descriptors:
                if descriptor.header['mol_id'] == id:
                    mol.SetProp(str(descriptor.header), str(descriptor.x))
        w = Chem.SDWriter(sd_file_path)
        for mol in self.suppl:
            w.write(mol)
        w.close()


class AtomDescriptorCreator(ChemDescriptorCreator):
    """ Concrete class that creates atom descriptors.
        create_property() is the factory method for atom descriptor properties.
    """
    def __init__(self, sdf_file_path):
        super(AtomDescriptorCreator, self).__init__(sdf_file_path)
        self._cast_sdf_mol_props_to_atoms()
        self._mol_index_list = []

    def _cast_sdf_mol_props_to_atoms(self):
        """ Assign atomic properties stored as molecular property to atom instances.

        A lot of atomic properties are stored as list of number as molecular properties in an SD File.
        This method casts them on atoms, such that later functions only needs to loop over atoms.

        :return: None
        """
        for molecule in self.suppl:
            # print molecule.GetProp('ID')
            atoms = molecule.GetAtoms()
            # make properties with plus and minus charge as well for Fukui calculations
            properties = []
            for pc_name in PCD.defined_partial_charges:
                properties.append(pc_name)
                properties.append(pc_name + '_minus')
                properties.append(pc_name + '_plus')
            properties += ['p orbital occupations', 'p orbital energies']
            for prop in properties:
                try:
                    if '\n' in molecule.GetProp(prop).rstrip('\n'):
                        # print  molecule.GetProp(prop)
                        props = molecule.GetProp(prop).split('\n')
                    else:
                        props = molecule.GetProp(prop).split()
                except KeyError:
                    'WARNING: Property {} not available'.format(prop)
                    continue
                for idx, a_prop in enumerate(props):
                    if prop in ('p orbital occupations', 'p orbital energies'):
                        tmp_a_prop = a_prop.split()[1:]
                        a_prop = ' '.join(tmp_a_prop)
                    atoms[idx].SetProp(prop, a_prop)
            # Set Atom property if it is SOM:
            som_indices = AtomDescriptorCreator._get_som_indices(molecule)
            for atom in atoms:
                if atom.GetIdx() in som_indices:
                    atom.SetProp('som', '1')
                else:
                    atom.SetProp('som', '0')

    @classmethod
    def _get_som_indices(cls, molecule):
        # Note that indices are already here converted to python standard, i.e. starting with 0 ...
        prop_names = molecule.GetPropNames()
        som_prop_names = [name for name in prop_names if 'SOM' in name]
        som_indices = set()
        for spn in som_prop_names:
            indices = molecule.GetProp(spn)
            for index in indices.split():
                som_indices.add(int(index)-1)
        return som_indices

    def create_property(self, prop_name, **options):
        # Das ist nicht wirklich factory pattern, da eine Factory das create als class-method hat.
        # das liesse sich so natuerlich auch implementieren, aber ich lasse das jetzt erstmal so
        # Man koennte anstatt konkrete strings auch class.__class__.__name__ benutzen
        """
        One could maybe put all into one statement of the sort:
        self.descriptor_properties.append(eval(prop_name({}).format(options))
        But I really want to avoid the use of eval at the moment ... However, it would result in really nicer
        encapsulation.

        Another possibility would be to delegate the creation process to the actual property class or property
        base (parent) class that could inspect the dynamic type and create an instance of "itself"
        """
        if prop_name == 'ChargeShell':
            self.descriptor_properties.append(ChargeShell(**options))
        elif prop_name == 'MassShell':
            self.descriptor_properties.append(MassShell(**options))
        elif prop_name == 'GraphChargeShell':
            self.descriptor_properties.append(GraphChargeShell(**options))
        elif prop_name == 'ChargeAutocorrelation':
            self.descriptor_properties.append(ChargeAutocorrelation(**options))
        elif prop_name == 'ChargeGraphAutocorrelation':
            self.descriptor_properties.append(ChargeGraphAutocorrelation(**options))
        elif prop_name == 'ChargeRDF':
            self.descriptor_properties.append(ChargeRDF(**options))
        else:
            raise Exception('Unknown descriptor element:  {}'.format(prop_name))

    def create_descriptor_vectors(self):
        if not self.descriptor_properties:
            raise Exception('No properties defined in DescriptorCreator. '
                            'Descriptor creation not possible')
        # create descriptor vectors only for possible EAS sites
        descriptor_index = 0
        for mol in self.suppl:
            self._mol_index_list.append([])
            eas_sites = mol.GetProp('EAS_atom_index').split('\n')
            eas_sites= [int(x)-1 for x in eas_sites]

            for atom in mol.GetAtoms():
                if atom.GetIdx() in eas_sites: 
                    # print(atom.GetIdx(), atom.GetSymbol())
                    self._mol_index_list[-1].append(descriptor_index)
                    descriptor_index += 1
                    atom_descriptor = AtomDescriptor(mol.GetProp('Mol_ID'), atom.GetIdx()+1)
                    atom_descriptor.header['elements'] = [prop.description for prop in self.descriptor_properties]
                    
                    # try:
                    #     atom_descriptor.y = list(map(int, mol.GetProp('EAS_atom_measured').split('\n')))[descriptor_index-1]
                    # except:
                    #     pass

                    for prop in self.descriptor_properties:
                        atom_descriptor.x += prop.calculate_elements(atom)
                    self.descriptors.append(atom_descriptor)


    def write_sdf(self, sd_file_path, descriptor_id='Descriptor'):
        for mol in self.suppl:
            ident = mol.GetProp('Mol_ID')
            prop_name = descriptor_id
            prop_content = ''
            for descriptor in self.descriptors:
                if descriptor.header['mol_id'] == ident:
                    prop_content += '{}'.format(descriptor.header['atom_idx'])
                    for ele in descriptor.x:
                        prop_content += ';{}'.format(ele)
                    prop_content += '\n'
            # print prop_name, prop_content
            mol.SetProp(prop_name, prop_content)
        w = Chem.SDWriter(sd_file_path)
        for mol in self.suppl:
            w.write(mol)
        w.close()
