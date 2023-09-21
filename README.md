# microbiomeR
Theoretical models for antibiotic resistance epidemiology in the healthcare setting. 

This R code and Mathematica notebook  support Smith et al. (2021), Microbiome-pathogen interactions drive epidemiological dynamics of antibiotic resistance: a modeling study applied to nosocomial pathogen control. eLife 2021. 
https://elifesciences.org/articles/68764
# about
We present a suite of ODE models describing bacterial colonization dynamics in the healthcare setting. Each model accounts for different within-host interactions:
* Model 1: bacterial colonization
* Model 2: exclusive colonization strain competition
* Model 3: microbiome-pathogen competition
* Model 4: strain-microbiome competition
* Model 5: strain-microbiome competition with interspecific horizontal gene transfer

Models are evaluated using the same parameter set corresponding to a generic pathogen (C^R). ODEs are integrated numerically. Epidemiological outcomes at population dynamic equilibrium are evaluated using univariate and bivariate analysis. 

# repository files (R)
* microbiomeR.Rproj
  * associated R project
* ODEs.R
  * ODEs for each model
* parameters.R
  * initial state variable vectors for each model
  * parameter vector for generic C^R
  * alternative parameter vectors for various analyses
* functions.R
  * functions returning epidemiological outcomes at population dynamic equilibrium: prevalence, incidence, resistance rate
* solve.R
  * execute analyses
* figures.R
  * render figures 

# repository files (Mathematica)
* microbiome_ecology.nb
  * notebook containing model ODEs, R0 expressions, numerical solutions,  figures

# contact
David Smith \
Institut Pasteur / Inserm / UVSQ \
david.smith@pasteur.fr \
davidrobertmundysmith@gmail.com

# citation
Smith DRM, Temime L, Opatowski L. Microbiome-pathogen interactions drive epidemiological dynamics of antibiotic resistance: a modelling study applied to nosocomial pathogen control. Elife. 2021 Sep 14;10:e68764. doi: 10.7554/eLife.68764. PMID: 34517942.
