# covid-19-immunity
Age structured difference equation model on the dynamics of SARS-CoV-2 that incorporates waning immunity. Code to accompany medRxiv preprint "Dynamics of SARS-CoV-2 with Waning Immunity in the UK Population" by Crellen et al. 2020.

After cloning the repository, run SEIIRRS-scenarios.R in RStudio. This script calls the functione in SEIIRRS-functions.R, which runs the age-structured discrete time model. Ensure that required R packages lubridate and testthat are installed prior to running the code. 

The code is the property of Thomas Crellen (thomas.crellen@bdi.ox.ac.uk). The age structured contact matrices for the UK population (in the data folder) are the property of Dr Petra Klepac (petra.klepac@lshtm.ac.uk).
