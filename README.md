# kinship
# English

Code corresponding to the simulations used for the manuscript "SIMULATING THE EFFECTS OF KINSHIP AND POSTMARITAL RESIDENCE PATTERNS ON MITOCHONDRIAL DNA DIVERSITY IN MORTUARY CONTEXTS" (AJPA-2023-00149). They simulate the behavior of mitochondrial DNA variants in human populations within mortuary contexts, following different kinship rules.

All scripts must be saved in the same folder, and the parameters of the simulations are set from the "mother simulation" script. Four parameters can be changed:

FERTILITY: normal distribution of live births, high (tfr=1) or low (tfr=2).

JUVENILE MORTALITY: normal distribution of prereproductive mortality, high (cm=1) or low (cm=2).

RELATIONSHIP NORM ADHERENCE: can be 0 (implying zero deviation from the norm) or a variable number from a normal distribution.

MIGRATION RATE: can be 0 or a variable number from a normal distribution.

Examples of simulation output are available in the "baseline results" and "mod results" folders. The output files are readable as dataframes in R (using the read.csv2 command) and contain the following variables:

run: simulation run number (1 to 500).

gen: generation number within each run (1 to 15)

nf (1 to 6): initial size (in number of females) of the group

tfr: fertility rate of the simulation (1 or 2)

cm: juvenile mortality rate of the simulation (1 or 2)

mr: migration rate, expressed as the proportion of females in the group coming from a population outside the set of 6 groups under consideration

ad: adherence to the kinship norm, expressed as the proportion of females in the group that DO NOT meet the prescribed kinship norms 

fnf (1 to 6): final size (in number of women) of the group at the end of the generation

divm (1 to 6): haplotype diversity of the male individuals in the group at the end of the generation

divf (1 to 6): haplotype diversity of the female individuals of the group at the end of the generation

divt (1 to 6): total haplotype diversity of the group at the end of the generation.

In the "samples" folders are the estimated diversities for samples of 10, 20, 30 and 50 individuals at the end of each run for the estimation of their usefulness in archaeological contexts.

Comments in the scripts are in Spanish. My apologies for not translating them; I believe online translating engines are good enough by now.
