import os
import sys
import numpy as np
import subprocess
import datetime

from rdkit import Chem
from rdkit.Chem import AllChem

from DescriptorCreator import AtomDescriptorCreator
from locate_EAS_sites import locate_sites
from find_atoms import find_identical_atoms

# xTB path and calc setup
# XTBHOME = '/home/Ree/anaconda/envs/latest_rdenv' #apparently this is not working for all molecules, fx. comp. with Reaction ID: 33222961 fails
# XTBPATH = '/home/Ree/anaconda/envs/latest_rdenv/share/xtb'
# MANPATH = '/home/Ree/anaconda/envs/latest_rdenv/share/man'
# LD_LIBRARY_PATH = '/home/Ree/anaconda/envs/latest_rdenv/lib'
XTBHOME = '/opt/xtb/6.4.0'
XTBPATH = '/opt/xtb/6.4.0/share/xtb'
MANPATH = '/opt/xtb/6.4.0/share/man'
LD_LIBRARY_PATH = '/opt/xtb/6.4.0/lib'

OMP_NUM_THREADS = '1'
MKL_NUM_THREADS = '1'


def generate_3Dsdf(smi, name='mol'):

    rdkit_mol = Chem.MolFromSmiles(smi)
    rdkit_mol = Chem.AddHs(rdkit_mol)

    ps = AllChem.ETKDGv3()
    ps.randomSeed = 123
    ps.useSmallRingTorsions=True
    AllChem.EmbedMolecule(rdkit_mol,ps)
    # AllChem.EmbedMolecule(rdkit_mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True, ETversion=2, randomSeed=90)
    # AllChem.UFFOptimizeMolecule(rdkit_mol)
    # AllChem.MMFFOptimizeMolecule(rdkit_mol)

    rdkit_mol.SetProp('Input_SMILES', smi)
    rdkit_mol.SetProp('Mol_ID', name)

    filename = name + '.3D.sdf'
    w = Chem.SDWriter(filename)
    w.write(rdkit_mol)
    w.close()

    return filename


class EASMolPreparator(object):
    """ Class to implement the preparation of finding Electrophilic Aromatic Substitution (EAS) sites in molecules 
    and write the corresponding properties to a SD file.
    """

    def __init__(self, sdf_file_path, identifier):
        """

        :param data_file_path: path to sdf file
        """
        self.sdf_file_path = sdf_file_path
        self.identifier = identifier
        # CM5/Mulliken Charges
        self.cm5 = []
        self.Mulliken = []


class EASMolPreparation(EASMolPreparator):
    """
    Class that implements the handling of the sdf files that already have the properties:
    3D coordinates in mol block and some identifyer property
    """

    def __init__(self, sdf_file_path, identifier):
        """

        :param sdf_file_path: path to sdf file
        :param mol_identifier: rdkit identifier property
        """
        super(EASMolPreparation, self).__init__(sdf_file_path,identifier)
        self.suppl = self._load_data(sdf_file_path)
        self.sdf_file_path = sdf_file_path
        self.identifier = identifier
        self.QCroot = self._make_QCroot()

        if not os.path.exists(self.QCroot):
            os.mkdir(self.QCroot)

        self.sdf_file_path = os.path.join(self.QCroot, sdf_file_path)
        os.replace(sdf_file_path, self.sdf_file_path) # move the initial SD file into the calc folder


    def _load_data(self, sdf_file_path):
        """
        Function to load molecules as .sdf mol supplier
        """
        suppl = Chem.SDMolSupplier(sdf_file_path, removeHs=False)
        print(f'{len(suppl)} mol(s) in supplier')
        molecule_list = []
        failed_mols = []
        for mol_idx, mol in enumerate(suppl):
            if not mol:
                failed_mols.append(mol_idx)
                print('Molecule number {} could not be processed.'.format(mol_idx + 1))
            else:
                molecule_list.append(mol)
            return suppl

    def _make_QCroot(self):
        """
        makes a pathname for the QC (xTB 6.4.0) calculations
        :return: QCroot
        """
        cwd = os.getcwd()
        QCroot = cwd + '/' + str(datetime.datetime.now()).split(' ')[0] + '-charges-xtb_6.4.0-calculations-to-descriptors'
        print(f'QC folder is: {QCroot}')
        return QCroot

    def sdf_to_qc_preparation(self,optimize=False):
        """
        function to split a sdf file and prepare run folders for quantum chemistry
        carries out xTB 6.4.0 single point localized orbitals (Foster-Boys) calculation
        Linux OS is mandatory
        :parameter: optimize: if set to true, a GFN1-xTB (xTB version 6.4.0) geometry optimization is triggered.
        """

        global XTBHOME
        global XTBPATH
        global MANPATH
        global LD_LIBRARY_PATH

        global OMP_NUM_THREADS
        global MKL_NUM_THREADS

        # Set env parameters for xTB
        os.environ['XTBHOME'] = XTBHOME
        os.environ['XTBPATH'] = XTBPATH
        os.environ['MANPATH'] = MANPATH
        os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH

        os.environ["OMP_NUM_THREADS"] = OMP_NUM_THREADS
        os.environ['MKL_NUM_THREADS'] = MKL_NUM_THREADS

        # Convert .sdf to .xyz
        mols = [x for x in self.suppl]
        print(f'{len(mols)} molecule(s) in supplier.')
        
        for mol in self.suppl:
            molid = mol.GetProp(self.identifier)
            #print "now doing molid", molid
            idx = molid
            path = f'{self.QCroot}/{idx}'
            if not os.path.exists(path):
                os.mkdir(path)
            nat = mol.GetNumAtoms()
            x = []
            y = []
            z = []
            aid = []
            xyz = []
            for atom in mol.GetAtoms():
                sym = atom.GetSymbol()
                aid.append(sym)
            for i in range(0, nat):
                position = mol.GetConformer().GetAtomPosition(i)
                x.append(position.x)
                y.append(position.y)
                z.append(position.z)
            xyz = np.column_stack((x, y, z))
            fout = open(f'{path}/coord.xyz', "w")
            fout.write(str(nat))
            fout.write("\n ")
            fout.write("\n ")
            for i in range(0, nat):
                fout.write("{} {:06.8f} {:06.8f} {:06.8f} ".format(aid[i], xyz[i, 0], xyz[i, 1], xyz[i, 2]))
                fout.write("\n ")
            fout.close()

            # Run xTB calc
            if optimize:
                cmd = f'{XTBHOME}/bin/xtb --gfn 1 coord.xyz --opt --lmo'
            else:
                cmd = f'{XTBHOME}/bin/xtb --gfn 1 coord.xyz --lmo'
            
            proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=path)
            output = proc.communicate()[0]
            
            # Save calc output
            with open(f'{path}/xtb.out', 'w') as f:
                f.write(output)


    def add_charges_and_eas_sites_sdf(self, measured):
        # Function to add Mulliken and cm5 charges as properties to SDF files
        output_sdf = f'{self.sdf_file_path[:-7]}-xtbcharges-eassites.3D.sdf'
        w = Chem.SDWriter(output_sdf)
        print(output_sdf)
        for mol in self.suppl:
            # Set charges
            mulliken = ''
            cm5 = ''
            idx = mol.GetProp(self.identifier)
            #print(f'Now doing {idx}')
            nat = int(mol.GetNumAtoms())
            fmulliken = open(self.QCroot + '/' + idx + '/charges')
            for line in fmulliken.readlines():
                mulliken = mulliken + line
            fxtb = open(self.QCroot + '/' + idx + '/xtb.out')
            lines = fxtb.readlines()
            for line in lines:
                if line.__contains__('Mulliken/CM5'):
                    mindex = lines.index(line)
            start = mindex + 1
            endindex = start + nat
            for i in range(start, endindex):
                cm5 = cm5 + lines[i].split()[2] + '\n'
            mol.SetProp('cm5', cm5)
            mol.SetProp('Mulliken', mulliken)

            # Set EAS sites (index start 1 instead of 0)
            eas_sites = sorted(locate_sites(mol)) # only unique atoms are considered 
            atom_symbols = []
            for atom in mol.GetAtoms():
                if atom.GetIdx() in eas_sites: 
                    atom_symbols.append(atom.GetSymbol())
            
            mol.SetProp('Atom', '\n'.join(atom_symbols))

            eas_sites_index_one = [str(elem + 1) for elem in eas_sites]
            mol.SetProp('EAS_atom_index', '\n'.join(eas_sites_index_one))

            EAS_atom_measured = []
            for site in eas_sites:
                if site in measured:
                    EAS_atom_measured.append('1')
                else:
                    EAS_atom_measured.append('0')
            mol.SetProp('EAS_atom_measured', '\n'.join(EAS_atom_measured))

            w.write(mol)
        w.close()


# Examples for descriptor types :
descriptordict = {'GACF_3_cm5': 
                                [('ChargeGraphAutocorrelation', {'charge_type': 'cm5', 'rmin': 1, 'rmax': 3})],
                  'CS_3_cm5': 
                                [('ChargeShell', {'charge_type': 'cm5', 'n_shells': 3})],
                  'CRDF_6_0.2_cm5': 
                                [('ChargeRDF', {'charge_type': 'cm5', 'rmin':1, 'rmax':6, 'beta':20, 'step_size': 0.2})],
                  'CACF_8_0.5_cm5':
                                [('ChargeAutocorrelation', {'charge_type': 'cm5', 'rmin': 1, 'rmax': 8, 'step_size': 0.5})],
                  'MS_2': 
                                [('MassShell', {'n_shells':2})],
                  'CGS_3_cm5':  
                                [('GraphChargeShell', {'charge_type': 'cm5', 'n_shells': 3, 'use_cip_sort': True})],
                  'combined_CACF_10_0.2_CGS_2_cm5_CS_7_cm5': 
                                [('GraphChargeShell', {'charge_type': 'cm5', 'n_shells': 2, 'use_cip_sort': True}),
                                ('ChargeAutocorrelation', {'charge_type': 'cm5', 'rmin': 1, 'rmax': 10, 'step_size': 0.2}),
                                ('ChargeShell', {'charge_type': 'cm5', 'n_shells': 7})]
                  }


if __name__ == "__main__":  
    
    smiles_filename = sys.argv[1]
    exam = True

    # Read data from input file
    f = open(smiles_filename, "r")

    comps_name = []
    comps_smiles = []
    comps_measured_atoms = []
    for line in f:
        words = line.split()
        comps_name.append(words[0])
        comps_smiles.append(words[1])

        if exam:
            measured_atoms = [int(x) for x in words[2].split(",")]
            comps_measured_atoms.append(find_identical_atoms(words[1], measured_atoms))
        else:
            comps_measured_atoms.append(None)
    f.close() 

    for i in range(len(comps_name)):
        try:
            name = comps_name[i]
            smi = comps_smiles[i]
            measured = comps_measured_atoms[i]

            filename = generate_3Dsdf(smi, name)
            identifier = 'Mol_ID'

            prepar = EASMolPreparation(filename,identifier) 
            prepar.sdf_to_qc_preparation(optimize=False) # optimize = True triggers GFN-xTB coordinate optimization.
            prepar.add_charges_and_eas_sites_sdf(measured)

            # Calculate the descriptor
            sdf = f'{prepar.QCroot}/{name}-xtbcharges-eassites.3D.sdf'
            sd_file_path = f'{sdf}.descriptors.sdf'
            for des_type in descriptordict:
                #print(os.path.isfile(sdf))
                adc = AtomDescriptorCreator(sdf)
                descriptor_elements = descriptordict[des_type]
                for des in descriptor_elements:
                    adc.create_property(des[0], **des[1])
                adc.create_descriptor_vectors()
                print(f'{des_type} descriptor has {len(adc.X[0, :])} dimensions.')
                adc.write_sdf(sd_file_path, descriptor_id=des_type)
                
                sdf = sd_file_path
                # xtrain = np.asarray(adc.X)
                # print(xtrain)
                # print(adc.y)
            print('\n')
        except:
            continue

