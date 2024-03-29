# covid-19-immunity
Age structured difference equation model on the dynamics of SARS-CoV-2 that incorporates waning immunity. Code to accompany medRxiv preprint "Dynamics of SARS-CoV-2 with Waning Immunity in the UK Population" by Crellen et al. 2020. https://www.medrxiv.org/content/10.1101/2020.07.24.20157982v1

After cloning the repository, run SEIIRRS-analysis.R in RStudio. This script calls the functions in SEIIRRS-functions.R, which runs the age-structured discrete time model. Ensure that required R packages lubridate, testthat & rsudioapi are installed prior to running the code. 

The code is the property of Thomas Crellen (thomas.crellen@bdi.ox.ac.uk). The age structured contact matrices for the UK population (in the data folder) are the property of Dr Petra Klepac (petra.klepac@lshtm.ac.uk).
