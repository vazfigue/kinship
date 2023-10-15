#### Simulacion de parentesco cognaticio con matrilocalidad, generacion 1

# PREVIOS:

# Uso el sufijo $bilateral en vez del nro de indice (2)
# para las estructuras de esta simulacion
# Es de lectura algo mas larga pero mas facil que sucesivos corchetes

# Lista de grupos de parentesco y sus nros de hijxs
grupos$bilateral<-list()
gruposhijxs$bilateral<-list()

# Lista de individuos en la nueva generacion
generacion$bilateral<-list()

# El vector cementerios ya fue definido en la simulacion madre

# Lista de generaciones de cementerios
piladehuesos$bilateral<-list()

# Los haplos internos (definidos en la simulacion madre) se reparten
# al azar con reposicion en seis vectores de largo definido por nf
	for(i in 1:6)
	{
	grupos$bilateral[[i]]<-sample(haplos_int,nf[i],replace=T)
	}

# Muestreo con reposicion de la fracción migrante del vector de haplotipos externos
	for(i in 1:6)
		{
		migrants<-sample(haplos_ext,
		                 round(length(grupos$bilateral[[i]])*mr),replace=T)
		locals<-sample(grupos$bilateral[[i]],
		               length(grupos$bilateral[[i]])-length(migrants),replace=F)
		grupos$bilateral[[i]]<-append(locals, migrants)
		}

# Muestreo al azar equivalente a la sobrevivencia en vida reproductiva
	for(i in 1:6)
		{
		survives<-runif(length(grupos$bilateral[[i]]),0,1)
		grupos$bilateral[[i]]<-grupos$bilateral[[i]][survives<fs]
		}

# Largamos la generacion 1:
# Generamos un vector de valores normales
# que determinara el numero de hijos que tendra cada mujer.
# los parametros se establecen en la simulacion madre

	for(i in 1:6)
		{
		gruposhijxs$bilateral[[i]]<-
		  round(abs(rnorm(length(grupos$bilateral[[i]]),
		                  fertility[tfr,1],fertility[tfr,2])))
		}

# multiplicacion por un vector de supervivencia juvenil
# definido en la simulacion madre
	for(i in 1:6)
		{
		gruposhijxs$bilateral[[i]]<-
		  round(gruposhijxs$bilateral[[i]]*
		          (1-abs(rnorm(length(gruposhijxs$bilateral[[i]]),
		                       mortality[cm,1],mortality[cm,2])))) 
		}

# Asignacion de haplos a las progenies

for(i in 1:6)
{ 
  generacion$bilateral[[i]]<-vector()
  # se multiplica cada haplo por cada vector de hijes
  # por ejemplo, si grupos[[1]] es c("A","B","C","D")
  # y si gruposhijxs[[1]] es c(2,3,2,3)
  # generacion[[1]] es c("A","A","B","B","B","C","C","D","D","D")
  for(j in 1:length(grupos$bilateral[[i]]))
  {
    progenie<-rep(grupos$bilateral[[i]][j],gruposhijxs$bilateral[[i]][j]) 
    generacion$bilateral[[i]]<-append(generacion$bilateral[[i]],progenie)
  } #cierre del bucle j  

  # se asigna sexo M y F al azar a les integrantes de cada generacion
  # y se divide cada generacion en subcomponentes $varon y $mujer en forma acorde
muestra<-sample(c("varon","mujer"),
                length(generacion$bilateral[[i]]), replace = T)
generacion$bilateral[[i]]<-split(generacion$bilateral[[i]],muestra)
} #cierre del bucle i

# Pasaje de la progenie a otro vector de varones o mujeres
# segun el patron de residencia en un vector "cementerio"

# COGNATICIO CON MATRILOCALIDAD:
# Los varones que nacieron ahi se entierran en otra parte
# en un loop anidado donde el grupo de pertenencia del varon
# queda excluido como destino
 for(i in 1:6)
{ 
 cementerios$bilateral[[i]]$mujer<-vector(length = 0) #seteo de vector vacio
 cementerios$bilateral[[i]]$varon<-vector(length = 0) #seteo de vector vacio
}

a<-1:6 #indice generico

for(i in 1:6)
{
  largos<-NULL
  for(j in a[-i])
  {
    largos<-append(largos, length(generacion$bilateral[[j]]$varon)) #vector de conteo de varones por grupo
  }#cierre bucle j (ojo, hay otro abajo) 
  largos<-largos/sum(largos) #conversion de los largos a proporciones del total
  
  sizes<-round(largos*length(generacion$bilateral[[i]]$varon)) 
  sizes<-replace(sizes,sizes==0,1) #mecanismo de seguridad por si algun elemento tiene valor cero
  generacion$bilateral[[i]]$varon<-sample(generacion$bilateral[[i]]$varon) #barajado de varones
  suppressWarnings(reparto<-split(generacion$bilateral[[i]]$varon,rep(seq_along(sizes), sizes)))
  #division en 5 partes, proporcionales al grupo de llegada
  #suprimo los warnings temporalmente porque a veces el reparto no suma igual al largo total
   k<-1
   for(j in a[-i])
  {
    cementerios$bilateral[[j]]$varon<-
      append(cementerios$bilateral[[j]]$varon,reparto[[k]]) #reparto en los demas grupos
    k<-k+1
  } # cierre bucle j
} # cierre bucle i

# con las mujeres es facil porque al ser matrilocal van al mismo cementerio...
  for(i in 1:6)
  {
    cementerios$bilateral[[i]]$mujer<-generacion$bilateral[[i]]$mujer
  }
