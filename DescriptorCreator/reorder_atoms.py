from rdkit import Chem

def single_bonded_mol(molecule):
    """ 
    Removes all double bonds, 
    but keeps the atoms in order. 
    """

    rd_mol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    for atom in molecule.GetAtoms():
        rd_atom = Chem.rdchem.Atom(atom.GetSymbol())
        rd_atom.SetAtomMapNum(atom.GetIdx())
        rd_mol.AddAtom(rd_atom)

    for atom in molecule.GetAtoms():
        idx1 = atom.GetIdx()
        for atomNeighbor in atom.GetNeighbors():
            idx2 = atomNeighbor.GetIdx()
            if idx1 < idx2:
                rd_mol.AddBond(idx1, idx2, Chem.rdchem.BondType.SINGLE)
    return rd_mol.GetMol()


def get_atoms_in_order(ref_mol, mols):
    """ 
    Takes a list of mol objects and 
    reorder the atoms according to a reference mol. 
    """
    
    reordered_mols = []
    single_bonded_ref = single_bonded_mol(ref_mol)
    for m in mols:
        atoms = single_bonded_mol(m).GetSubstructMatch(single_bonded_ref)
        reordered_mols.append(Chem.rdmolops.RenumberAtoms(m, atoms))
    
    return reordered_mols


if __name__ == "__main__":

    ref = Chem.MolFromSmiles('c1c(c2cc(sc2)C)n[nH]c1')
    mol_taut_1 = Chem.MolFromSmiles('c1c(-c2cc(C)sc2)n[nH]c1')
    mol_taut_2 = Chem.MolFromSmiles('c1c(-c2cc(C)sc2)[nH]nc1')

    mols = [mol_taut_1, mol_taut_2]

    new_taut_1, new_taut_2 = get_atoms_in_order(ref, mols)
    new_smiles = [Chem.MolToSmiles(new_taut_1, canonical=False), Chem.MolToSmiles(new_taut_2, canonical=False)]
    print(new_smiles)