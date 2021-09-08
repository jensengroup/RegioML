
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import PrepareMolForDrawing
from rdkit.Chem import rdDepictor
#rdDepictor.SetPreferCoordGen(True)

from reorder_atoms import get_atoms_in_order
from collections import defaultdict

# Drawing Options
color_predicted = (0.2, 1, 0.0) # Green
color_loseicted = (1.0, 0.1, 0.3) # Red
color_measured = (0.0, 0.0, 0.0) # Black
arad = 0.25
#molsPerRow = 4
subImgSize = (300,300)


def draw2d(mol, name, atom_scores_predicted, atom_scores_loseicted, subImgSize, highlight_predicted, highlight_loseicted, measure=None):
    
    global color_predicted
    global color_loseicted
    global color_measured
    global arad
    
    d2d = rdMolDraw2D.MolDraw2DSVG(subImgSize[0], subImgSize[1])
    dos = d2d.drawOptions()
    dos.atomHighlightsAreCircles = False
    dos.fillHighlights = False

    #label atoms with scores
    for i, idx in enumerate(highlight_predicted):
        dos.atomLabels[idx] = '{:.0f}%'.format(atom_scores_predicted[i])

    for i, idx in enumerate(highlight_loseicted):
        dos.atomLabels[idx] = '{:.0f}%'.format(atom_scores_loseicted[i])
    
    #highlight atoms
    atomHighlighs = defaultdict(list)
    highlightRads = {}
    for idx in highlight_predicted:
        atomHighlighs[idx].append(color_predicted)
        highlightRads[idx] = arad

    # did threshold find some predictions?
    # find ones not in predicted list
    highlight_loseicted = list(set(highlight_loseicted)-set(highlight_predicted))
    if len(highlight_loseicted):
        for idx in highlight_loseicted:
            atomHighlighs[idx].append(color_loseicted)
            highlightRads[idx] = arad

    if measure:
        for idx in measure:
            atomHighlighs[idx].append(color_measured)
            highlightRads[idx] = arad
    
    [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
    mol = PrepareMolForDrawing(mol)
    d2d.DrawMoleculeWithHighlights(mol, name, dict(atomHighlighs), {}, highlightRads, {})
    d2d.FinishDrawing()

    return d2d.GetDrawingText()


def generate_structure(ref_smi, smiles, names, atom_scores, predicted, highlight_measure=None):
    
    global subImgSize
    #global molsPerRow
    molsPerRow = 4

    highlight_predicted, highlight_loseicted = predicted
    atom_scores_predicted, atom_scores_loseicted = atom_scores

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
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    mols = get_atoms_in_order(Chem.MolFromSmiles(ref_smi), mols)

    cwidth = 0
    cheight = 0
    drawed_mols = []
    for i in range(len(smiles)):
        res = draw2d(mols[i], names[i], atom_scores_predicted[i], atom_scores_loseicted[i], subImgSize, highlight_predicted[i], highlight_loseicted[i], highlight_measure)
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

    svg = generate_structure('n1ccc[nH]1', ['n1ccc[nH]1'], None, [[[0.3]],[[0.2]]], [[[1]],[[0]]], highlight_measure=None)
    print(svg)