# adpProject
This repository includes codes related to a paper named "An approximate dynamic programming approach for sequential pig marketing decisions at herd level". 

The codes have been written in C++ programming lanuage envirenment and are compiled using MS VS2010 compiler. 

The compiler needs to the following external libraries: 
- Armadillo C++ Linear Algebra Library (http://arma.sourceforge.net)
- Boost C++ Library - Release 1.55.0 (http://www.boost.org/users/history/version_1_55_0.html)
- ILOG Concert Technology C++ library (http://www-eio.upc.edu/lceio/manuals/cplex75/doc/concert12/doc/userman/html/preface.html)

Notes: 
1- Only the main files (e.g. *.h and *.cpp ) have been commited in the repository. This has been done based on the standard format of gitigonore for visual studio (see https://github.com/github/gitignore/blob/master/VisualStudio.gitignore).

2- You need to copy folder "csv_fils" in a directory in your own pc. Next you need to specify the address of this directory in functions "readIniSlopes()", "readFinalSlopes()", "storeSlopes()"  in the file adp.cpp .
A same procedure should be done in the main file "main.cpp" for the csv file "policyInfo.csv" in folder "csv_fils". 

3- Results of numerical experiments in the paper have been reported as csv files in folder "paper_results".
