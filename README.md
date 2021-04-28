# Dynamics of SARS-CoV-2 with waning immunity in the UK population
Code to accompany an age-structured difference equation model on the dynamics of SARS-CoV-2 that incorporates waning immunity. Preprint currently hosted at https://www.medrxiv.org/content/10.1101/2020.07.24.20157982v2 (doi:10.1101/2020.07.24.20157982)

After cloning the repository, run `SEIIRRS-analysis.R` in _RStudio_. This script calls the functions in `SEIIRRS-functions.R`, which runs the age-structured discrete time model. Ensure that required _R_ packages `lubridate`, `testthat `& `rstudioapi` are installed prior to running the code. 

Questions about the code can be directed to the corresponding author, Thomas Crellen (thomas.crellen@bdi.ox.ac.uk). The age structured contact matrices for the UK population (in the `/data` folder) are the property of Dr Petra Klepac (petra.klepac@lshtm.ac.uk).

# Citing this code resource
As these analysis notes are an integral part of the paper and do not have their own DOI, we kindly request that you cite the full paper as follows:

> @ARTICLE{Crellen2020,
 author={Crellen, T and Pi, L and Davis, EL and Pollington, TM and Lucas, TCDL and Ayabina, D and Borlase, A and Toor, J and Prem, K and Medley, GF and Klepac, P and Hollingsworth, TD},
 year={2020},  
 title={{Dynamics of SARS-CoV-2 with waning immunity in the UK population}},  
 journal={Phil. Trans. R. Soc. B},  
 doi={10.1098/rstb.2020.0274},  
}

[![CC BY 4.0][cc-by-shield]][cc-by]  

This work is licensed under a [Creative Commons Attribution 4.0 International License][cc-by].  

[cc-by]: http://creativecommons.org/licenses/by/4.0/  
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
