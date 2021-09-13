# RegioML
RegioML is an atom-based machine learning model for predicting the regioselectivities of electrophilic aromatic substitution reactions. The model relies on CM5 atomic charges computed using semiempirical tight binding (GFN1-xTB) combined with the ensemble decision tree variant light gradient boosting machine (LightGBM).

More information is available in the [RegioML paper](https://doi.org/). 
Furthermore, additional code for dataset curation, descriptors, and machine learning etc. are available [here](https://sid.erda.dk/sharelink/HypB1igzDl).

# Installation

We recommend using anaconda to install the Python 3 environment:

    conda env create -f environment.yml && conda activate regioml

Then download the binaries of xtb version 6.4.0:

    mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.4.0/xtb-210201.tar.xz; tar -xvf ./xtb-210201.tar.xz; cd ..


# Usage

An example of the command line use of RegioML:

    # Run predictions:

    python regioML.py -s 'c1(ccno1)C'


    # Run predictions, specify name, and highlight observed reaction sites (black circles):

    python regioML.py -s 'c1ccc(cc1C(F)(F)F)c1cccc(n1)OC' -n 'mol-1' -o '13,11'

The results are then printed and viewable as a 2D structure with regioselective indicators (in .svg format).

The atom scores are predicted by the best LightGBM classification model from training on the entire collection of data using 10-fold cross-validation ('models/LGBM_measured_allData_final_model.txt'), where values above 50% indicates that an atom should be reactive (green circles). However, atoms with scores above 5% are also highlighted (red circles). The predicted low, medium, or high reactivity are based on the highest proton affinity within the molecule obtained by the best LightGBM regression model from training on the entire collection of data using 10-fold cross-validation ('models/LGBM_regressor_GFN1_allData_final_model.txt'). 

The classification model used in the paper to obtain the results in Table 1 can be accessed using the following command line:

    python regioML.py -s 'c1(ccno1)C' -m 'models/LGBM_random_noTaut_spCM5_measured_final_best_model.txt'


HAPPY PREDICTING :-)
