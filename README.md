# kinship
# Español

Código correspondiente a las simulaciones utilizadas para el manuscrito "SIMULATING THE EFFECTS OF KINSHIP AND POSTMARITAL RESIDENCE PATTERNS ON MITOCHONDRIAL DNA DIVERSITY IN MORTUARY CONTEXTS" (AJPA-2023-00149). Simulan el comportamiento de variantes de ADN mitocondrial en poblaciones humanas en contextos mortuorios, obedeciendo a distintas normas de parentesco.

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
