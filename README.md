# microbiomeR
Theoretical models for antibiotic resistance epidemiology. R code and a Mathematica notebook to support Smith et al., (2021). Microbiome-pathogen interactions drive epidemiological dynamics of antibiotic resistance: modelling insights for infection control

# about
We present a suite of ODE models describing bacterial colonization dynamics in the healthcare setting. Each model accounts for different within-host interactions:
* Model 1: bacterial colonization
* Model 2: exclusive colonization strain competition
* Model 3: microbiome-pathogen competition
* Model 4: strain-microbiome competition
* Model 5: strain-microbiome competition with interspecific horizontal gene transfer

Models are evaluated using the same parameter set corresponding to a generic pathogen (C^R). ODEs are integrated numerically. Epidemiological outcomes at population dynamic equilibrium are evaluated using univariate and bivariate analysis. 

# repository files (R)
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
