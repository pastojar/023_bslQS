This repository stores the code used for the analysis that led to writing my paper in the Journal of Hydrology (https://doi.org/10.1016/j.jhydrol.2022.129051).



It contains the following subfolders and files:

- R-script files (*.r) named bslQS_1***.r - bslQS_6***.r. The numbers represent the order in which the scripts are executed. The files and numberes represent logical units 
(1 - WAA calibration; 2 - R-R modelling; 3 - statistical evaluation; 4 - creating boxplots; 5 - creating scatterplots; 6 - uncertainty analysis). 

- \inst\awk - short awk scripts for SWMM files manipuation

- \R - scripts containing functions called by the scripts in the mother folder



In order to perfrom the analysis, additional files are necessery which are not stored here:

- \inst\extdata - a collection of rainfall and runoff data, including some basic statistics, used in the analyses  +  external apps (swmm5.exe and gawk.exe)

- \inst\swmm\ - files relevant for SWMM modelling (original input file, its modified version, rainfall time series files, output file)
