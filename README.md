# kinship
# Español

Código correspondiente a las simulaciones utilizadas para el manuscrito "SIMULATING THE EFFECTS OF KINSHIP AND POSTMARITAL RESIDENCE PATTERNS ON MITOCHONDRIAL DNA DIVERSITY IN MORTUARY CONTEXTS" (DOI: 10.1002/ajpa.24910). Simulan el comportamiento de variantes de ADN mitocondrial en poblaciones humanas en contextos mortuorios, obedeciendo a distintas normas de parentesco.

Todos los scripts deben guardarse en la misma carpeta, y los parámetros de las simulaciones se establecen a partir del script "simulacion madre". Se pueden cambiar cuatro parámetros:

FERTILIDAD: distribución normal de nacidos vivos, alta (tfr=1) o baja (tfr=2)

MORTALIDAD JUVENIL: distribución normal de mortalidad prerreproductiva, alta (cm=1) o baja (cm=2)

ADHESION A LA NORMA DE PARENTESCO: puede ser 0 (lo que implica cero desviación de la norma) o un número variable a partir de una distribución normal.

TASA DE MIGRACIÓN: puede ser 0 o un número variable a partir de una distribución normal.

Ejemplos de la salida de las simulaciones están disponibles en las carpetas "baseline results" y "mod results". Los archivos de salida son legibles como dataframes en R (usando el comando read.csv2) y contienen las siguientes variables:

run: número de corrida de la simulación (1 a 500)

gen: número de generación de la corrida (1 a 15)

nf (1 a 6): tamaño inicial (en número de mujeres) del grupo

tfr: tasa de fertilidad de la simulación (1 o 2)

cm: tasa de mortalidad juvenil de la simulación (1 o 2)

mr: tasa migratoria, expresada como la proporción de mujeres del grupo provenientes de una población externa al conjuto de 6 grupos considerado

ad: adhesión a la norma de parentesco, expresada como la proporción de mujeres del grupo que NO cumplen las normas de parentesco prescriptas 

fnf (1 a 6): tamaño final (en número de mujeres) del grupo al cabo de la generación

divm (1 a 6): diversidad de haplotipos de los individuos masculinos del grupo al cabo de la generación

divf (1 a 6): diversidad de haplotipos de los individuos femeninos del grupo al cabo de la generación

divt (1 a 6): diversidad total de haplotipos del grupo al cabo de la generación

En las carpetas "samples" se encuentran las diversidades estimadas para muestras de 10, 20, 30 y 50 individuos del final de cada corrida para la estimación de su utilidad en contextos arqueológicos.

# English

Code corresponding to the simulations used for the manuscript "SIMULATING THE EFFECTS OF KINSHIP AND POSTMARITAL RESIDENCE PATTERNS ON MITOCHONDRIAL DNA DIVERSITY IN MORTUARY CONTEXTS" (DOI: 10.1002/ajpa.24910). They simulate the behavior of mitochondrial DNA variants in human populations within mortuary contexts, following different kinship rules.

All scripts must be saved in the same folder, and the parameters of the simulations are set from the script in the "simulacion madre" file. Four parameters can be changed:

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
