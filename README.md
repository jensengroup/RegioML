# RegioML
RegioML is an atom-based machine learning model for predicting the regioselectivities of electrophilic aromatic substitution reactions. The model relies on CM5 atomic charges computed using semiempirical tight binding (GFN1-xTB) combined with the ensemble decision tree variant light gradient boosting machine (LightGBM).

More information is available at the [RegioML paper](https://doi.org/).

# Installation

We recommend using anaconda to install the Python 3 environment:
```conda create -n regioml python=3.8.8 && conda activate sites```

OR

```conda create -n regioml python=3.8.8 && source activate sites```

Then install dependencies

```conda install -n regioml -c rdkit rdkit && conda install -n regioml pip```
```pip install -r requirements.txt```

# Usage

An example of the command line use of RegioML:

    # Run predictions:

    python regioml.py -s 'c1(ccno1)C'


    # Run predictions, specify name, and highlight the experimentally observed reaction sites (black circles):

    python regioML.py -s 'c1ccc(cc1C(F)(F)F)c1cccc(n1)OC' -n 'mol-1' -o '13,11'

The results are then printed and viewable by 2D structures with regioselective indicators (in .svg format).


