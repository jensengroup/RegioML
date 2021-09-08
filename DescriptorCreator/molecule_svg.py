
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import PrepareMolForDrawing
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)

from DescriptorCreator.reorder_atoms import get_atoms_in_order
from DescriptorCreator.find_atoms import find_identical_atoms, find_identical_atoms_with_scores
from collections import defaultdict

# Drawing Options
color_predicted = (0.2, 1, 0.0) # Green
color_loseicted = (1.0, 0.1, 0.3) # Red
color_measured = (0.0, 0.0, 0.0) # Black
arad = 0.25
subImgSize = (300,300)


def draw2d(mol, name, atom_scores, subImgSize, highlight_predicted, measure=None):
    
    global color_predicted
    global color_loseicted
    global color_measured
    global arad
    
    d2d = rdMolDraw2D.MolDraw2DSVG(subImgSize[0], subImgSize[1])
    dos = d2d.drawOptions()
    dos.atomHighlightsAreCircles = False
    dos.fillHighlights = False

    #label atoms with scores
    atomHighlighs = defaultdict(list)
    highlightRads = {}

    highlight_predicted, atom_scores = find_identical_atoms_with_scores(mol, highlight_predicted, atom_scores)
    for i, idx in enumerate(highlight_predicted):
        atom_score = int(round(atom_scores[i]*100))
        if atom_score >= 50:
            dos.atomLabels[idx] = '{}%'.format(atom_score)
            atomHighlighs[idx].append(color_predicted)
            highlightRads[idx] = arad
        elif atom_score > 5:
            dos.atomLabels[idx] = '{}%'.format(atom_score)
            atomHighlighs[idx].append(color_loseicted)
            highlightRads[idx] = arad

    #highlight measured atoms
    if measure:
        for idx in measure:
            atomHighlighs[idx].append(color_measured)
            highlightRads[idx] = arad
    
    [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
    mol = PrepareMolForDrawing(mol)
    d2d.DrawMoleculeWithHighlights(mol, name, dict(atomHighlighs), {}, highlightRads, {})
    d2d.FinishDrawing()

    return d2d.GetDrawingText()


def generate_structure(ref_smi, smiles, names, atom_scores, highlight_predicted, highlight_measure=None, molsPerRow=4):

    global subImgSize

    if names == None:
        names = ['' for i in range(len(smiles))]

    nRows = len(smiles) // molsPerRow
    if len(smiles) % molsPerRow:
        nRows += 1
    if nRows == 1:
        molsPerRow = len(smiles)
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])

    header = """<svg version='1.1' baseProfile='full'
                xmlns='http://www.w3.org/2000/svg'
                        xmlns:rdkit='http://www.rdkit.org/xml'
                        xmlns:xlink='http://www.w3.org/1999/xlink'
                    xml:space='preserve'
    width='{0}px' height='{1}px' viewBox='0 0 {0} {1}'>
    <!-- END OF HEADER -->""".format(fullSize[0],fullSize[1])

    spacer = '<g transform="translate({0},{1})">\n{2}</g>'

    ### Make sure the atoms are in order ###
    ref_mol = Chem.MolFromSmiles(ref_smi)
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    mols = get_atoms_in_order(ref_mol, mols)

    if highlight_measure:
        highlight_measure = [int(x) for x in highlight_measure.split(',')]
        highlight_measure = find_identical_atoms(ref_mol, highlight_measure)
    
    cwidth = 0
    cheight = 0
    drawed_mols = []
    for i in range(len(smiles)):
        res = draw2d(mols[i], names[i], atom_scores[i], subImgSize, highlight_predicted[i], highlight_measure)
        res = res.split("\n")
        end_of_header = res.index("<!-- END OF HEADER -->") + 1 
        res = "\n".join(res[end_of_header:-2])
        
        res = "".join(spacer.format(int(cwidth*subImgSize[0]), int(cheight*subImgSize[1]), res))
        drawed_mols.append(res)
        
        if int(i+1) % molsPerRow == 0 and i != 0:
            cheight += 1
            cwidth = 0
        elif molsPerRow == 1:
            cheight += 1
            cwidth = 0
        else:
            cwidth += 1

    svg = header + "\n" + "\n".join(drawed_mols) + "\n</svg>"

    return svg


if __name__ == "__main__":

    smiles = ['CCCN1C2=C(C(OC)=CC=C2)[C+](C3=C(OC)C=CC=C3OC)C4=C1C=CC=C4OC']
    names = ['test_mol']
    highlight_predicted = [[9, 10, 11, 17, 18]]
    atom_scores = [[2.01352495e-01, 7.17073891e-05, 4.52574807e-01, 9.77132978e-01,
 2.00897332e-04]]
    highlight_measure = [17]
    
    result_svg = generate_structure(smiles[0], smiles, names, atom_scores, highlight_predicted, highlight_measure)
    
    fd = open('test.svg','w')
    fd.write(result_svg)
    fd.close()
    


