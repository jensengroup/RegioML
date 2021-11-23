import os
import subprocess
import datetime

from rdkit import Chem
from rdkit.Chem import AllChem

from DescriptorCreator.locate_EAS_sites import locate_sites
from DescriptorCreator.GraphChargeShell import GraphChargeShell

# xTB path and calc setup
path = os.getcwd()
XTBHOME = os.path.join(path, 'dep/xtb-6.4.0')
XTBPATH = os.path.join(path, 'dep/xtb-6.4.0/share/xtb')
MANPATH = os.path.join(path, 'dep/xtb-6.4.0/share/man')
LD_LIBRARY_PATH = os.path.join(path, 'dep/xtb-6.4.0/lib')

OMP_NUM_THREADS = '1'
MKL_NUM_THREADS = '1'



class EASMolPreparation():
    """
    Class to implement the preparation of finding 
    Electrophilic Aromatic Substitution (EAS) sites in molecules.
    """

    def __init__(self):
       
        # Make seperate directory for descritptor calculations
        self.SQMroot = self._make_SQMroot()
        if not os.path.exists(self.SQMroot):
            os.mkdir(self.SQMroot)

        # Set env parameters for xTB        
        global XTBHOME
        global XTBPATH
        global MANPATH
        global LD_LIBRARY_PATH
        global OMP_NUM_THREADS
        global MKL_NUM_THREADS
        os.environ['XTBHOME'] = XTBHOME
        os.environ['XTBPATH'] = XTBPATH
        os.environ['MANPATH'] = MANPATH
        os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
        os.environ["OMP_NUM_THREADS"] = OMP_NUM_THREADS
        os.environ['MKL_NUM_THREADS'] = MKL_NUM_THREADS


    def _make_SQMroot(self):
        """
        makes a pathname for the SQM calculations (xTB 6.4.0)
        :return: SQMroot
        """
        cwd = os.getcwd()
        SQMroot = cwd + '/' + str(datetime.datetime.now()).split(' ')[0] + '-charges-xtb_6.4.0-calculations-to-descriptors'
        print(f'SQM folder is: \n{SQMroot}')
        return SQMroot


    def generate_3Dxyz(self, smi, name):

        # Smiles to RDKit mol object
        self.rdkit_mol = Chem.MolFromSmiles(smi)
        self.rdkit_mol = Chem.AddHs(self.rdkit_mol)

        # Embed mol object to get cartesian coordinates
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 123
        ps.useSmallRingTorsions=True
        if AllChem.EmbedMolecule(self.rdkit_mol, ps) == -1:
            print(f'1st embed failed for {name} with SMILES: {smi}; will try useRandomCoords=True')
            ps = AllChem.ETKDGv3()
            ps.randomSeed = 123
            ps.useSmallRingTorsions=True
            ps.useRandomCoords=True #added 21/6 - 2021
            # ps.maxIterations=1000 #added 21/6 - 2021
            if AllChem.EmbedMolecule(self.rdkit_mol, ps) == -1:
                raise Exception(f'2nd embed failed for {name} with SMILES: {smi}')

        # Make seperate directory for mol
        self.mol_calc_path = f'{self.SQMroot}/{name}'
        if not os.path.exists(self.mol_calc_path):
            os.mkdir(self.mol_calc_path)
        
        # Mol object to xyz file
        self.xyz_file_path = f'{self.mol_calc_path}/{name}.xyz'
        Chem.rdmolfiles.MolToXYZFile(self.rdkit_mol, self.xyz_file_path)


    def calc_CM5_charges(self, smi, name='pred_mol', optimize=False, save_output=False):
        """
        function to handling the quantum chemistry calculation
        carries out xTB 6.4.0 single point localized orbitals (Foster-Boys) calculation
        Linux OS is mandatory
        :parameter: optimize: if set to true, a GFN1-xTB (xTB version 6.4.0) geometry optimization is triggered.
        """

        # Generate xyz file from SMILES
        self.generate_3Dxyz(smi, name)

        # Get molecule properties
        chrg = Chem.GetFormalCharge(self.rdkit_mol)
        spin = 0 #spin hardcoded to zero

        # Run xTB calc
        if optimize:
            cmd = f'{XTBHOME}/bin/xtb --gfn 1 {self.xyz_file_path} --opt --lmo --chrg {chrg} --uhf {spin}' #TODO! add connectivity check!
        else:
            cmd = f'{XTBHOME}/bin/xtb --gfn 1 {self.xyz_file_path} --lmo --chrg {chrg} --uhf {spin}'
        
        proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=self.mol_calc_path)
        output = proc.communicate()[0]
        
        # Save calc output
        if save_output:
            with open(f'{self.mol_calc_path}/xtb.out', 'w') as f:
                f.write(output)

        # Get CM5 charges from output and append CM5 charges to RDKit mol object 
        cm5_list = []
        natoms = int(self.rdkit_mol.GetNumAtoms())
        for line_idx, line in enumerate(output.split('\n')):
            if 'Mulliken/CM5' in line:
                start = line_idx + 1
                endindex = start + natoms
                for i in range(start, endindex):
                    line = output.split('\n')[i]
                    cm5_atom = float(line.split()[2])
                    cm5_list.append(cm5_atom)
                break
        
        for i, atom in enumerate(self.rdkit_mol.GetAtoms()):
            atom.SetProp('cm5', str(cm5_list[i]))
        
        return cm5_list


    def create_descriptor_vector(self, prop_name, **options):
        """
        Create the GraphChargeShell descriptor 
        for all unique EAS sites in a molecule.
        :parameter: prop_name example: 'GraphChargeShell'
        :parameter: options example: {'charge_type': 'cm5', 'n_shells': 5, 'use_cip_sort': True}
        """
        
        if prop_name == 'GraphChargeShell':
            self.descriptor_properties = GraphChargeShell(**options)
        else:
            raise Exception(f'Unknown descriptor element: {prop_name}')

        # Create descriptor vector only for possible EAS sites     
        eas_sites = sorted(locate_sites(self.rdkit_mol)) # only unique atoms are considered 
        
        atom_indices = []
        descriptor_vector = []
        for atom in self.rdkit_mol.GetAtoms():
            if atom.GetIdx() in eas_sites: 
                atom_indices.append(atom.GetIdx())

                atom_descriptor = self.descriptor_properties.calculate_elements(atom)
                descriptor_vector.append(atom_descriptor)
        
        return atom_indices, descriptor_vector
