#### Simulacion de parentesco matrilineal, generacion 1

# PREVIOS:

# Uso el sufijo $matrilineal en vez del nro de indice (1)
# para las estructuras de esta simulacion
# Es de lectura algo mas larga pero mas facil que sucesivos corchetes

# Lista de grupos de parentesco y sus nros de hijxs
grupos$matrilineal<-list()
gruposhijxs$matrilineal<-list()

# Lista de individuos en la nueva generacion
generacion$matrilineal<-list()

# El vector cementerios ya fue definido en la simulacion madre

# Lista de generaciones de cementerios
piladehuesos$matrilineal<-list()

# Los haplos internos (definidos en la simulacion madre) se reparten
# al azar con reposicion en seis vectores de largo definido por nf
	for(i in 1:6)
	{
	grupos$matrilineal[[i]]<-sample(haplos_int,nf[i],replace=T)
	}

# Muestreo con reposicion de la fracción migrante del vector de haplotipos externos
	for(i in 1:6)
		{
		migrants<-sample(haplos_ext,
		                 round(length(grupos$matrilineal[[i]])*mr),replace=T)
		locals<-sample(grupos$matrilineal[[i]],
		               length(grupos$matrilineal[[i]])-length(migrants),replace=F)
		grupos$matrilineal[[i]]<-append(locals, migrants)
		}

# Muestreo al azar equivalente a la sobrevivencia en vida reproductiva
	for(i in 1:6)
		{
		survives<-runif(length(grupos$matrilineal[[i]]),0,1)
		grupos$matrilineal[[i]]<-grupos$matrilineal[[i]][survives<fs]
		}

# Largamos la generacion 1:
# Generamos un vector de valores normales
# que determinara el numero de hijos que tendra cada mujer.
# los parametros se establecen en la simulacion madre

	for(i in 1:6)
		{
		gruposhijxs$matrilineal[[i]]<-
		  round(abs(rnorm(length(grupos$matrilineal[[i]]),
		                  fertility[tfr,1],fertility[tfr,2])))
		}

# multiplicacion por un vector de supervivencia juvenil
# definido en la simulacion madre
	for(i in 1:6)
		{
		gruposhijxs$matrilineal[[i]]<-
		  round(gruposhijxs$matrilineal[[i]]*
		          (1-abs(rnorm(length(gruposhijxs$matrilineal[[i]]),
		                       mortality[cm,1],mortality[cm,2])))) 
		}

# Asignacion de haplos a las progenies

for(i in 1:6)
{ 
  generacion$matrilineal[[i]]<-vector()
  # se multiplica cada haplo por cada vector de hijes
  # por ejemplo, si grupos[[1]] es c("A","B","C","D")
  # y si gruposhijxs[[1]] es c(2,3,2,3)
  # generacion[[1]] es c("A","A","B","B","B","C","C","D","D","D")
  for(j in 1:length(grupos$matrilineal[[i]]))
  {
    progenie<-rep(grupos$matrilineal[[i]][j],gruposhijxs$matrilineal[[i]][j]) 
    generacion$matrilineal[[i]]<-append(generacion$matrilineal[[i]],progenie)
  } #cierre del bucle j  

  # se asigna sexo M y F al azar a les integrantes de cada generacion
  # y se divide cada generacion en subcomponentes $varon y $mujer en forma acorde
muestra<-sample(c("varon","mujer"),
                length(generacion$matrilineal[[i]]), replace = T)
generacion$matrilineal[[i]]<-split(generacion$matrilineal[[i]],muestra)
} #cierre del bucle i

# Pasaje de la progenie a otro vector de varones o mujeres
# segun el patron de residencia en un vector "cementerio"

# MATRILINEALIDAD/MATRILOCALIDAD:
# Este es sin dudas el mas sencillito:
# Entierro en cementerio de matrilinaje, en cuyo caso el cementerio
# repite el patron de la generacion
for(i in 1:6)
{
cementerios$matrilineal[[i]]<-generacion$matrilineal[[i]]
}
