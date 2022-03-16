# COVID

This repository contains two script files that simulate the initial spread of COVID-19 in a population, focused on Aotearoa New Zealand. The code uses a stochastic/branching process model and is appropriate for the initial stages of an outbreak when there are a small number of cases and random chance plays a role in determining if an outbreak spreads or self-extinguishes. There are two script files:

* **stochastModelDelta.m**,  which considers the Delta variant, accounts for vaccinated and unvaccinated individuals, calculates asymptomatic and symptomatic infections and hospitalization rates. This version can consider different population level controls. 
* **stochastModelOmicron.m**, which considers the Omicron variant, accounts for boosted, vaccinated and unvaccinated individuals, and calculates symptomatic infections. 

The script files can be run for different parameters to consider the impact of different vaccination rates, population controls (for Delta), and vaccine effectiveness (i.e., to examine waning vaccine effectiveness). 

The code is written in MATLAB and, for small numbers of realizations, the code can be run efficiently on a standard desktop/laptop computer. Computational cost scales linearly with the number of realizations, so for large numbers of realizations (i.e., > 10,000) it can be more practical to run on a computing cluster, particular if considering multiple different vaccination rates/initial seed conditions. 

For more details about the mathematical model and results, see:
* Watson, L. M. (2022) Simulating the impact of vaccination rates on the initial stages of a COVID-19 outbreak in New Zealand (Aotearoa) with a stochastic model, _New Zealand Medical Journal_. Available from medRxiv: [https://doi.org/10.1101/2021.11.22.21266721](https://doi.org/10.1101/2021.11.22.21266721)
* Watson, L. M. (2022) Likelihood of infecting or getting infected with COVID-19 as a function of vaccination status, as investigated with a stochastic model for Aotearoa New Zealand for Delta and Omicron variants, _New Zealand Medical Journal_. Available from medRxiv: [https://doi.org/10.1101/2021.11.28.21266967](https://doi.org/10.1101/2021.11.28.21266967).

The first release, corresponding to the version used in the papers above, is archived at Zenodo with a DOI: [![DOI](https://zenodo.org/badge/176619362.svg)](https://zenodo.org/badge/latestdoi/176619362)

### How do I get set up? ###
* Clone this respository to your local directory and run the script file
* Codes to analyze simulation results, and previous simulation results, are availabe at [https://osf.io/gnc9v/](https://osf.io/gnc9v/). 

### Who do I talk to? ###
* Leighton Watson: leighton.watson@canterbury.ac.nz
