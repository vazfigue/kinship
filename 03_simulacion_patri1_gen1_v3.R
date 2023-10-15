#### Simulacion de parentesco patrilineal con reclutamiento, generacion 1

# PREVIOS:

# Uso el sufijo $patrilineal1 en vez del nro de indice (3)
# para las estructuras de esta simulacion
# Es de lectura algo mas larga pero mas facil que sucesivos corchetes

# Lista de grupos de parentesco y sus nros de hijxs
grupos$patrilineal1<-list()
gruposhijxs$patrilineal1<-list()

# Lista de individuos en la nueva generacion
generacion$patrilineal1<-list()

# El vector cementerios ya fue definido en la simulacion madre

# Lista de generaciones de cementerios
piladehuesos$patrilineal1<-list()

# Los haplos internos (definidos en la simulacion madre) se reparten
# al azar con reposicion en seis vectores de largo definido por nf
	for(i in 1:6)
	{
	grupos$patrilineal1[[i]]<-sample(haplos_int,nf[i],replace=T)
	}

# Muestreo con reposicion de la fracción migrante del vector de haplotipos externos
	for(i in 1:6)
		{
		migrants<-sample(haplos_ext,
		                 round(length(grupos$patrilineal1[[i]])*mr),replace=T)
		locals<-sample(grupos$patrilineal1[[i]],
		               length(grupos$patrilineal1[[i]])-length(migrants),replace=F)
		grupos$patrilineal1[[i]]<-append(locals, migrants)
		}

# Muestreo al azar equivalente a la sobrevivencia en vida reproductiva
	for(i in 1:6)
		{
		survives<-runif(length(grupos$patrilineal1[[i]]),0,1)
		grupos$patrilineal1[[i]]<-grupos$patrilineal1[[i]][survives<fs]
		}

# Largamos la generacion 1:
# Generamos un vector de valores normales
# que determinara el numero de hijos que tendra cada mujer.
# los parametros se establecen en la simulacion madre

	for(i in 1:6)
		{
		gruposhijxs$patrilineal1[[i]]<-
		  round(abs(rnorm(length(grupos$patrilineal1[[i]]),
		                  fertility[tfr,1],fertility[tfr,2])))
		}

# multiplicacion por un vector de supervivencia juvenil
# definido en la simulacion madre
	for(i in 1:6)
		{
		gruposhijxs$patrilineal1[[i]]<-
		  round(gruposhijxs$patrilineal1[[i]]*
		          (1-abs(rnorm(length(gruposhijxs$patrilineal1[[i]]),
		                       mortality[cm,1],mortality[cm,2])))) 
		}

# Asignacion de haplos a las progenies

for(i in 1:6)
{ 
  generacion$patrilineal1[[i]]<-vector()
  # se multiplica cada haplo por cada vector de hijes
  # por ejemplo, si grupos[[1]] es c("A","B","C","D")
  # y si gruposhijxs[[1]] es c(2,3,2,3)
  # generacion[[1]] es c("A","A","B","B","B","C","C","D","D","D")
  for(j in 1:length(grupos$patrilineal1[[i]]))
  {
    progenie<-rep(grupos$patrilineal1[[i]][j],gruposhijxs$patrilineal1[[i]][j]) 
    generacion$patrilineal1[[i]]<-append(generacion$patrilineal1[[i]],progenie)
  } #cierre del bucle j  

  # se asigna sexo M y F al azar a les integrantes de cada generacion
  # y se divide cada generacion en subcomponentes $varon y $mujer en forma acorde
muestra<-sample(c("varon","mujer"),
                length(generacion$patrilineal1[[i]]), replace = T)
generacion$patrilineal1[[i]]<-split(generacion$patrilineal1[[i]],muestra)
} #cierre del bucle i

# Pasaje de la progenie a otro vector de varones o mujeres
# segun el patron de residencia en un vector "cementerio"

# PATRILINEALIDAD/PATRILOCALIDAD CON RECLUTAMIENTO/COGNATICIO CON PATRILOCALIDAD:
# Es similar a la inversa del cognaticio,
# ya que son las mujeres nacidas las que se entierran en otro sitio.
# en un loop anidado donde el grupo de pertenencia de la mujer
# queda excluido como destino
for(i in 1:6)
{ 
  cementerios$patrilineal1[[i]]$mujer<-vector() #seteo de vector vacio
  cementerios$patrilineal1[[i]]$varon<-vector() #seteo de vector vacio
}

a<-1:6 #indice generico

for(i in 1:6)
{
  largos<-NULL
  for(j in a[-i])
  {
    largos<-append(largos, length(generacion$patrilineal1[[j]]$mujer)) #vector de conteo de mujeres por grupo
  }#cierre bucle j (ojo, hay otro abajo) 
  largos<-largos/sum(largos) #conversion de los largos a proporciones del total
  
  sizes<-round(largos*length(generacion$patrilineal1[[i]]$mujer)) 
  sizes<-replace(sizes,sizes==0,1) #mecanismo de seguridad por si algun elemento tiene valor cero
  generacion$patrilineal1[[i]]$mujer<-sample(generacion$patrilineal1[[i]]$mujer) #barajado de mujeres
  suppressWarnings(reparto<-split(generacion$patrilineal1[[i]]$mujer,rep(seq_along(sizes), sizes)))
  #division en 5 partes, proporcionales al grupo de llegada
  #suprimo los warnings temporalmente porque a veces el reparto no suma igual al largo total
   k<-1
   for(j in a[-i])
  {
		cementerios$patrilineal1[[j]]$mujer<-
		  append(cementerios$patrilineal1[[j]]$mujer,reparto[[k]]) #reparto en los demas grupos
		k<-k+1
	} # cierre bucle j
} # cierre bucle i

# Y con los varones es pan comido ya que se quedan en el cementerio de su patrilinaje...
for(i in 1:6)
{
  cementerios$patrilineal1[[i]]$varon<-generacion$patrilineal1[[i]]$varon
}

