# RegioML
RegioML is an atom-based machine learning model for predicting the regioselectivities of electrophilic aromatic substitution reactions. The model relies on CM5 atomic charges computed using semiempirical tight binding (GFN1-xTB) combined with the ensemble decision tree variant light gradient boosting machine (LightGBM).

More information is available in the [RegioML paper](https://doi.org/).

# Installation

We recommend using anaconda to install the Python 3 environment:

```conda env create -f environment.yml && conda activate regioml```

Then download the binaries of xtb version 6.4.0:

```mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.4.0/xtb-210201.tar.xz; tar -xvf ./xtb-210201.tar.xz; cd ..```


# Usage

An example of the command line use of RegioML:

    # Run predictions:

    python regioML.py -s 'c1(ccno1)C'


    # Run predictions, specify name, and highlight observed reaction sites (black circles):

    python regioML.py -s 'c1ccc(cc1C(F)(F)F)c1cccc(n1)OC' -n 'mol-1' -o '13,11'

The results are then printed and viewable as a 2D structure with regioselective indicators (in .svg format).

<svg version='1.1' baseProfile='full'
                xmlns='http://www.w3.org/2000/svg'
                        xmlns:rdkit='http://www.rdkit.org/xml'
                        xmlns:xlink='http://www.w3.org/1999/xlink'
                    xml:space='preserve'
    width='300px' height='300px' viewBox='0 0 300 300'>
    <!-- END OF HEADER -->
<g transform="translate(0,0)">
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300' height='300' x='0' y='0'> </rect>
<ellipse cx='266.915' cy='153.986' rx='19.4483' ry='6.27598'  style='fill:none;stroke:#33FF00;stroke-width:8.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<ellipse cx='206.313' cy='178.395' rx='15.1784' ry='6.27598'  style='fill:none;stroke:#FF194C;stroke-width:8.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 252.708,158.703 L 241.205,178.351' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 246.066,158.772 L 238.014,172.525' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-20 atom-5 atom-0' d='M 241.503,129.016 L 252.975,149.166' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 241.205,178.351 L 217.758,178.209' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 209.913,173.245 L 198.629,153.425' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 213.171,167.454 L 205.272,153.58' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 198.629,153.425 L 204.435,143.507' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 204.435,143.507 L 210.242,133.589' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-3 atom-8' d='M 198.629,153.425 L 186.273,153.351' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-3 atom-8' d='M 186.273,153.351 L 173.917,153.276' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 216.791,128.867 L 229.147,128.941' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 229.147,128.941 L 241.503,129.016' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 220.463,134.586 L 229.113,134.638' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 229.113,134.638 L 237.762,134.69' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 241.503,129.016 L 247.259,119.184' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 247.259,119.184 L 253.016,109.351' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 260.309,104.461 L 272.343,104.534' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 272.343,104.534 L 284.377,104.607' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 167.444,148.508 L 161.749,138.504' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 161.749,138.504 L 156.053,128.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-21 atom-13 atom-8' d='M 155.755,177.835 L 161.561,167.917' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-21 atom-13 atom-8' d='M 161.561,167.917 L 167.368,157.999' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-10' d='M 156.053,128.5 L 127.57,128.328' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-10 atom-11' d='M 127.57,128.328 L 121.763,138.246' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-10 atom-11' d='M 121.763,138.246 L 115.957,148.164' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-11 atom-12' d='M 115.881,157.655 L 121.576,167.659' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-11 atom-12' d='M 121.576,167.659 L 127.272,177.663' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-11 atom-14' d='M 109.408,152.887 L 97.0518,152.812' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-11 atom-14' d='M 97.0518,152.812 L 84.6958,152.737' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-12 atom-13' d='M 127.272,177.663 L 155.755,177.835' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 82.2377,151.298 L 76.4448,161.194' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 76.4448,161.194 L 70.6518,171.089' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 87.1539,154.177 L 81.361,164.072' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 81.361,164.072 L 75.5681,173.967' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-16' d='M 84.6958,152.737 L 79.0489,142.819' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-16' d='M 79.0489,142.819 L 73.4021,132.901' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-16 atom-17' d='M 66.1881,127.958 L 54.1539,127.885' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-16 atom-17' d='M 54.1539,127.885 L 42.1197,127.812' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-17 atom-17 atom-18' d='M 42.1197,127.812 L 42.2917,99.3289' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-18 atom-17 atom-19' d='M 42.1197,127.812 L 41.9477,156.296' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-19 atom-17 atom-20' d='M 42.1197,127.812 L 13.6364,127.64' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-0' d='M 253.608 157.079
L 255.374 157.079
L 255.374 151.052
L 253.425 151.656
L 253.163 150.984
L 255.636 149.878
L 256.445 150.015
L 256.445 157.079
L 258.028 157.079
L 258.028 157.991
L 253.608 157.991
L 253.608 157.079
' fill='#000000'/>
<path  class='atom-0' d='M 261.674 158.082
Q 260.193 158.082, 259.453 156.988
Q 258.723 155.894, 258.723 153.946
Q 258.723 151.998, 259.453 150.915
Q 260.182 149.833, 261.674 149.833
Q 263.167 149.833, 263.896 150.915
Q 264.625 151.998, 264.625 153.946
Q 264.625 155.894, 263.885 156.988
Q 263.155 158.082, 261.674 158.082
M 261.674 157.17
Q 262.552 157.17, 263.019 156.361
Q 263.486 155.541, 263.486 153.946
Q 263.486 152.362, 263.019 151.553
Q 262.552 150.744, 261.674 150.744
Q 260.808 150.744, 260.33 151.553
Q 259.863 152.362, 259.863 153.946
Q 259.863 155.541, 260.33 156.361
Q 260.808 157.17, 261.674 157.17
' fill='#000000'/>
<path  class='atom-0' d='M 268.271 158.082
Q 266.79 158.082, 266.049 156.988
Q 265.32 155.894, 265.32 153.946
Q 265.32 151.998, 266.049 150.915
Q 266.779 149.833, 268.271 149.833
Q 269.764 149.833, 270.493 150.915
Q 271.222 151.998, 271.222 153.946
Q 271.222 155.894, 270.482 156.988
Q 269.752 158.082, 268.271 158.082
M 268.271 157.17
Q 269.148 157.17, 269.616 156.361
Q 270.083 155.541, 270.083 153.946
Q 270.083 152.362, 269.616 151.553
Q 269.148 150.744, 268.271 150.744
Q 267.405 150.744, 266.927 151.553
Q 266.46 152.362, 266.46 153.946
Q 266.46 155.541, 266.927 156.361
Q 267.405 157.17, 268.271 157.17
' fill='#000000'/>
<path  class='atom-0' d='M 273.455 157.968
L 278.32 149.548
L 279.118 150.004
L 274.253 158.424
L 273.455 157.968
M 273.808 154.538
Q 272.92 154.538, 272.464 153.912
Q 272.02 153.274, 272.02 152.169
Q 272.02 151.063, 272.464 150.448
Q 272.92 149.821, 273.808 149.821
Q 274.709 149.821, 275.153 150.448
Q 275.597 151.063, 275.597 152.169
Q 275.597 153.285, 275.153 153.912
Q 274.709 154.538, 273.808 154.538
M 273.808 153.718
Q 274.23 153.718, 274.458 153.331
Q 274.686 152.932, 274.686 152.169
Q 274.686 151.405, 274.458 151.018
Q 274.23 150.63, 273.808 150.63
Q 273.387 150.63, 273.159 151.029
Q 272.931 151.417, 272.931 152.169
Q 272.931 152.932, 273.159 153.331
Q 273.387 153.718, 273.808 153.718
M 278.879 158.093
Q 277.99 158.093, 277.534 157.467
Q 277.09 156.828, 277.09 155.723
Q 277.09 154.618, 277.534 154.003
Q 277.99 153.376, 278.879 153.376
Q 279.779 153.376, 280.223 154.003
Q 280.667 154.618, 280.667 155.723
Q 280.667 156.84, 280.223 157.467
Q 279.779 158.093, 278.879 158.093
M 278.879 157.273
Q 279.3 157.273, 279.528 156.885
Q 279.756 156.487, 279.756 155.723
Q 279.756 154.96, 279.528 154.573
Q 279.3 154.185, 278.879 154.185
Q 278.457 154.185, 278.229 154.584
Q 278.001 154.971, 278.001 155.723
Q 278.001 156.487, 278.229 156.885
Q 278.457 157.273, 278.879 157.273
' fill='#000000'/>
<path  class='atom-2' d='M 199.454 178.15
Q 200.24 178.378, 200.616 178.891
Q 201.004 179.392, 201.004 180.19
Q 201.004 180.873, 200.662 181.409
Q 200.32 181.933, 199.693 182.229
Q 199.067 182.514, 198.246 182.514
Q 197.381 182.514, 196.731 182.218
Q 196.093 181.91, 195.58 181.295
L 196.23 180.634
Q 196.731 181.181, 197.153 181.397
Q 197.574 181.602, 198.246 181.602
Q 198.976 181.602, 199.42 181.215
Q 199.864 180.816, 199.864 180.178
Q 199.864 179.358, 199.397 178.993
Q 198.941 178.617, 197.95 178.617
L 197.369 178.617
L 197.369 177.797
L 197.882 177.797
Q 198.759 177.786, 199.226 177.41
Q 199.693 177.022, 199.693 176.304
Q 199.693 175.78, 199.306 175.473
Q 198.919 175.154, 198.258 175.154
Q 197.586 175.154, 197.164 175.393
Q 196.754 175.632, 196.435 176.236
L 195.649 175.814
Q 195.934 175.142, 196.617 174.698
Q 197.301 174.242, 198.258 174.242
Q 199.443 174.242, 200.138 174.8
Q 200.833 175.359, 200.833 176.304
Q 200.833 176.954, 200.48 177.421
Q 200.126 177.888, 199.454 178.15
' fill='#000000'/>
<path  class='atom-2' d='M 204.65 182.491
Q 203.168 182.491, 202.428 181.397
Q 201.699 180.304, 201.699 178.355
Q 201.699 176.407, 202.428 175.325
Q 203.157 174.242, 204.65 174.242
Q 206.142 174.242, 206.871 175.325
Q 207.601 176.407, 207.601 178.355
Q 207.601 180.304, 206.86 181.397
Q 206.131 182.491, 204.65 182.491
M 204.65 181.58
Q 205.527 181.58, 205.994 180.771
Q 206.461 179.95, 206.461 178.355
Q 206.461 176.772, 205.994 175.963
Q 205.527 175.154, 204.65 175.154
Q 203.784 175.154, 203.305 175.963
Q 202.838 176.772, 202.838 178.355
Q 202.838 179.95, 203.305 180.771
Q 203.784 181.58, 204.65 181.58
' fill='#000000'/>
<path  class='atom-2' d='M 209.834 182.377
L 214.699 173.957
L 215.496 174.413
L 210.631 182.833
L 209.834 182.377
M 210.187 178.948
Q 209.298 178.948, 208.842 178.321
Q 208.398 177.683, 208.398 176.578
Q 208.398 175.473, 208.842 174.857
Q 209.298 174.231, 210.187 174.231
Q 211.087 174.231, 211.531 174.857
Q 211.976 175.473, 211.976 176.578
Q 211.976 177.694, 211.531 178.321
Q 211.087 178.948, 210.187 178.948
M 210.187 178.127
Q 210.608 178.127, 210.836 177.74
Q 211.064 177.341, 211.064 176.578
Q 211.064 175.814, 210.836 175.427
Q 210.608 175.04, 210.187 175.04
Q 209.765 175.04, 209.537 175.439
Q 209.31 175.826, 209.31 176.578
Q 209.31 177.341, 209.537 177.74
Q 209.765 178.127, 210.187 178.127
M 215.257 182.503
Q 214.368 182.503, 213.913 181.876
Q 213.468 181.238, 213.468 180.133
Q 213.468 179.027, 213.913 178.412
Q 214.368 177.786, 215.257 177.786
Q 216.157 177.786, 216.601 178.412
Q 217.046 179.027, 217.046 180.133
Q 217.046 181.249, 216.601 181.876
Q 216.157 182.503, 215.257 182.503
M 215.257 181.682
Q 215.679 181.682, 215.906 181.295
Q 216.134 180.896, 216.134 180.133
Q 216.134 179.369, 215.906 178.982
Q 215.679 178.595, 215.257 178.595
Q 214.835 178.595, 214.608 178.993
Q 214.38 179.381, 214.38 180.133
Q 214.38 180.896, 214.608 181.295
Q 214.835 181.682, 215.257 181.682
' fill='#000000'/>
<path  class='atom-4' d='M 211.237 124.811
L 213.88 129.083
Q 214.142 129.505, 214.564 130.268
Q 214.985 131.032, 215.008 131.077
L 215.008 124.811
L 216.079 124.811
L 216.079 132.877
L 214.974 132.877
L 212.137 128.206
Q 211.806 127.659, 211.453 127.032
Q 211.111 126.406, 211.009 126.212
L 211.009 132.877
L 209.961 132.877
L 209.961 124.811
L 211.237 124.811
' fill='#0000FF'/>
<path  class='atom-6' d='M 252.191 104.457
Q 252.191 102.521, 253.148 101.438
Q 254.105 100.356, 255.894 100.356
Q 257.683 100.356, 258.64 101.438
Q 259.597 102.521, 259.597 104.457
Q 259.597 106.417, 258.628 107.534
Q 257.66 108.639, 255.894 108.639
Q 254.116 108.639, 253.148 107.534
Q 252.191 106.429, 252.191 104.457
M 255.894 107.727
Q 257.124 107.727, 257.785 106.907
Q 258.457 106.075, 258.457 104.457
Q 258.457 102.874, 257.785 102.076
Q 257.124 101.267, 255.894 101.267
Q 254.663 101.267, 253.991 102.065
Q 253.33 102.862, 253.33 104.457
Q 253.33 106.087, 253.991 106.907
Q 254.663 107.727, 255.894 107.727
' fill='#FF0000'/>
<path  class='atom-8' d='M 168.363 149.22
L 171.006 153.493
Q 171.268 153.914, 171.69 154.678
Q 172.111 155.441, 172.134 155.487
L 172.134 149.22
L 173.205 149.22
L 173.205 157.287
L 172.1 157.287
L 169.263 152.615
Q 168.932 152.068, 168.579 151.442
Q 168.237 150.815, 168.135 150.621
L 168.135 157.287
L 167.087 157.287
L 167.087 149.22
L 168.363 149.22
' fill='#0000FF'/>
<path  class='atom-11' d='M 111.396 148.876
L 114.039 153.149
Q 114.301 153.57, 114.723 154.334
Q 115.145 155.097, 115.167 155.143
L 115.167 148.876
L 116.238 148.876
L 116.238 156.943
L 115.133 156.943
L 112.296 152.271
Q 111.966 151.725, 111.613 151.098
Q 111.271 150.471, 111.168 150.278
L 111.168 156.943
L 110.12 156.943
L 110.12 148.876
L 111.396 148.876
' fill='#0000FF'/>
<path  class='atom-15' d='M 66.6023 177.342
Q 66.6023 175.405, 67.5594 174.322
Q 68.5164 173.24, 70.3052 173.24
Q 72.094 173.24, 73.051 174.322
Q 74.0081 175.405, 74.0081 177.342
Q 74.0081 179.301, 73.0396 180.418
Q 72.0712 181.523, 70.3052 181.523
Q 68.5278 181.523, 67.5594 180.418
Q 66.6023 179.313, 66.6023 177.342
M 70.3052 180.612
Q 71.5357 180.612, 72.1965 179.791
Q 72.8687 178.959, 72.8687 177.342
Q 72.8687 175.758, 72.1965 174.96
Q 71.5357 174.151, 70.3052 174.151
Q 69.0747 174.151, 68.4025 174.949
Q 67.7416 175.746, 67.7416 177.342
Q 67.7416 178.971, 68.4025 179.791
Q 69.0747 180.612, 70.3052 180.612
' fill='#FF0000'/>
<path  class='atom-16' d='M 66.9002 128.007
Q 66.9002 126.07, 67.8572 124.988
Q 68.8143 123.905, 70.6031 123.905
Q 72.3918 123.905, 73.3489 124.988
Q 74.306 126.07, 74.306 128.007
Q 74.306 129.967, 73.3375 131.083
Q 72.3691 132.188, 70.6031 132.188
Q 68.8257 132.188, 67.8572 131.083
Q 66.9002 129.978, 66.9002 128.007
M 70.6031 131.277
Q 71.8336 131.277, 72.4944 130.457
Q 73.1666 129.625, 73.1666 128.007
Q 73.1666 126.423, 72.4944 125.626
Q 71.8336 124.817, 70.6031 124.817
Q 69.3726 124.817, 68.7003 125.614
Q 68.0395 126.412, 68.0395 128.007
Q 68.0395 129.636, 68.7003 130.457
Q 69.3726 131.277, 70.6031 131.277
' fill='#FF0000'/>
<path  d='M 116.06 284.744
L 116.06 283.592
L 117.628 283.592
L 117.884 281.192
L 119.004 281.192
L 119.004 283.592
L 121.516 283.592
L 121.516 284.744
L 119.004 284.744
L 119.004 289.288
Q 119.004 290.712, 120.204 290.712
Q 120.684 290.712, 121.388 290.488
L 121.644 291.544
Q 120.764 291.912, 119.932 291.912
Q 118.86 291.912, 118.172 291.272
Q 117.5 290.632, 117.5 289.368
L 117.5 284.744
L 116.06 284.744
' fill='#000000'/>
<path  d='M 122.38 287.688
Q 122.38 285.672, 123.388 284.568
Q 124.396 283.448, 126.236 283.448
Q 128.044 283.448, 128.844 284.536
Q 129.66 285.608, 129.66 287.64
L 129.66 287.896
L 123.916 287.896
Q 123.948 289.288, 124.572 290.024
Q 125.196 290.76, 126.364 290.76
Q 127.004 290.76, 127.596 290.616
Q 128.188 290.456, 128.908 290.136
L 129.34 291.16
Q 128.524 291.576, 127.788 291.768
Q 127.052 291.944, 126.284 291.944
Q 124.428 291.944, 123.404 290.824
Q 122.38 289.704, 122.38 287.688
M 126.236 284.632
Q 125.292 284.632, 124.716 285.176
Q 124.156 285.72, 123.98 286.776
L 128.076 286.776
Q 127.964 285.672, 127.516 285.16
Q 127.068 284.632, 126.236 284.632
' fill='#000000'/>
<path  d='M 130.924 290.12
Q 131.644 290.424, 132.156 290.568
Q 132.684 290.712, 133.244 290.712
Q 133.964 290.712, 134.364 290.408
Q 134.78 290.104, 134.78 289.512
Q 134.78 289.144, 134.556 288.92
Q 134.348 288.68, 134.076 288.568
Q 133.82 288.456, 133.1 288.216
Q 132.988 288.184, 132.252 287.928
Q 131.516 287.672, 131.116 287.128
Q 130.732 286.584, 130.732 285.784
Q 130.732 284.792, 131.5 284.12
Q 132.268 283.448, 133.756 283.448
Q 134.38 283.448, 134.924 283.592
Q 135.484 283.736, 136.092 283.976
L 135.692 285.16
Q 135.164 284.952, 134.7 284.824
Q 134.252 284.696, 133.756 284.696
Q 133.004 284.696, 132.62 285
Q 132.236 285.288, 132.236 285.768
Q 132.236 286.216, 132.54 286.44
Q 132.844 286.648, 133.516 286.888
Q 133.692 286.936, 133.836 287
L 134.22 287.144
Q 134.892 287.368, 135.292 287.592
Q 135.708 287.816, 135.996 288.28
Q 136.284 288.728, 136.284 289.48
Q 136.284 290.664, 135.436 291.32
Q 134.604 291.96, 133.26 291.96
Q 132.476 291.96, 131.804 291.8
Q 131.132 291.624, 130.46 291.336
L 130.924 290.12
' fill='#000000'/>
<path  d='M 136.748 284.744
L 136.748 283.592
L 138.316 283.592
L 138.572 281.192
L 139.692 281.192
L 139.692 283.592
L 142.204 283.592
L 142.204 284.744
L 139.692 284.744
L 139.692 289.288
Q 139.692 290.712, 140.892 290.712
Q 141.372 290.712, 142.076 290.488
L 142.332 291.544
Q 141.452 291.912, 140.62 291.912
Q 139.548 291.912, 138.86 291.272
Q 138.188 290.632, 138.188 289.368
L 138.188 284.744
L 136.748 284.744
' fill='#000000'/>
<path  d='M 150.492 292.616
L 150.492 293.928
L 142.332 293.928
L 142.332 292.616
L 150.492 292.616
' fill='#000000'/>
<path  d='M 161.308 283.448
Q 162.7 283.448, 163.404 284.2
Q 164.108 284.936, 164.108 286.376
L 164.108 291.816
L 162.604 291.816
L 162.604 286.488
Q 162.604 285.512, 162.22 285.08
Q 161.852 284.632, 161.02 284.632
Q 160.3 284.632, 159.676 284.984
Q 159.068 285.32, 158.7 285.928
Q 158.716 286.072, 158.716 286.376
L 158.716 291.816
L 157.212 291.816
L 157.212 286.488
Q 157.212 285.512, 156.828 285.08
Q 156.46 284.632, 155.628 284.632
Q 154.908 284.632, 154.3 284.984
Q 153.692 285.32, 153.308 285.912
L 153.308 291.816
L 151.804 291.816
L 151.804 283.592
L 153.036 283.592
L 153.196 284.712
Q 154.268 283.448, 155.916 283.448
Q 157.852 283.448, 158.46 284.888
Q 159.548 283.448, 161.308 283.448
' fill='#000000'/>
<path  d='M 164.844 287.688
Q 164.844 285.688, 165.852 284.568
Q 166.876 283.448, 168.748 283.448
Q 170.62 283.448, 171.628 284.568
Q 172.652 285.688, 172.652 287.688
Q 172.652 289.688, 171.628 290.824
Q 170.62 291.944, 168.748 291.944
Q 166.876 291.944, 165.852 290.824
Q 164.844 289.688, 164.844 287.688
M 166.38 287.688
Q 166.38 289.192, 166.988 289.976
Q 167.612 290.76, 168.748 290.76
Q 169.884 290.76, 170.492 289.976
Q 171.116 289.192, 171.116 287.688
Q 171.116 286.184, 170.492 285.416
Q 169.884 284.632, 168.748 284.632
Q 167.612 284.632, 166.988 285.416
Q 166.38 286.184, 166.38 287.688
' fill='#000000'/>
<path  d='M 174.028 279.672
L 175.5 279.672
L 175.5 291.816
L 174.028 291.816
L 174.028 279.672
' fill='#000000'/></g>
</svg>

HAPPY PREDICTING :-)
