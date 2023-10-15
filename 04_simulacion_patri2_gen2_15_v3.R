#### Simulacion de parentesco patrilineal sin reclutamiento, generaciones 2 a 15

# PREVIOS:

# Las estructuras de datos relevantes ya fueron definidas en la
# simulacion madre y en la generacion 1

# Uso el sufijo $patrilineal2 en vez del nro de indice (4)
# para las estructuras de esta simulacion
# Es de lectura algo mas larga pero mas facil que sucesivos corchetes

# Repartimos las mujeres de la generacion previa en la lista grupos
# ATENTO: aca es donde a menudo crashea el script con un subscript out of bounds
# por no haber mujeres
# Para ver si lo evito agrego una llave de seguridad
# que agregue UN haplo al azar si alguno de los vectores $mujer esta vacio
for(i in 1:6)
{
  if (length(generacion$patrilineal2[[i]]$mujer)==0){
    grupos$patrilineal2[[i]]<-sample(haplos_int,1)
} else {
      grupos$patrilineal2[[i]]<-generacion$patrilineal2[[i]]$mujer
}
}

# Muestreo con reposicion de la fracción migrante del vector de haplotipos externos
	for(i in 1:6)
		{
		migrants<-sample(haplos_ext,
		                 round(length(grupos$patrilineal2[[i]])*mr),replace=T)
		locals<-sample(grupos$patrilineal2[[i]],
		               length(grupos$patrilineal2[[i]])-length(migrants),replace=F)
		grupos$patrilineal2[[i]]<-append(locals, migrants)
		}

# Muestreo al azar equivalente a la sobrevivencia en vida reproductiva
	for(i in 1:6)
		{
		survives<-runif(length(grupos$patrilineal2[[i]]),0,1)
		grupos$patrilineal2[[i]]<-grupos$patrilineal2[[i]][survives<fs]
		}

# Largamos la generacion 1:
# Generamos un vector de valores normales
# que determinara el numero de hijos que tendra cada mujer.
# los parametros se establecen en la simulacion madre

	for(i in 1:6)
		{
		gruposhijxs$patrilineal2[[i]]<-
		  round(abs(rnorm(length(grupos$patrilineal2[[i]]),
		                  fertility[tfr,1],fertility[tfr,2])))
		}

# multiplicacion por un vector de supervivencia juvenil
# definido en la simulacion madre
	for(i in 1:6)
		{
		gruposhijxs$patrilineal2[[i]]<-
		  round(gruposhijxs$patrilineal2[[i]]*
		          (1-abs(rnorm(length(gruposhijxs$patrilineal2[[i]]),
		                       mortality[cm,1],mortality[cm,2])))) 
		}

### ACA ES DONDE LA COSA SE PONE LINDA ### porque

# PATRILINEALIDAD/PATRILOCALIDAD SIN RECLUTAMIENTO
# La mujer da a luz hijos que se quedan en el patrilinaje del esposo
# y ella es enterrada con su patrilinaje
# Por lo tanto, lo logico es que la descendencia de CADA MUJER de CADA grupo sea redirigida
# a otro sitio sin mas tramite
# La forma ideal de operacionalizar eso es:
# 1) hacer un sorteo de a que cementerio va la descendencia de cada mujer (asumimos azar)
# 2) mandar a toda la descendencia para ese cementerio
# Ese sorteo se hace no al momento de asignacion de haplos, sino con el vector de hijos vivos

for(i in 1:6)
{ 
  preasignacion$madres[[i]]<-vector() #seteo de vector vacio
  preasignacion$nrohijos[[i]]<-vector() #seteo de vector vacio
}

a<-1:6 #indice generico
semillas<-c(3,11,22,36,44,57) #vector de seeds para el barajado
for(i in 1:6)
{
  largos<-NULL
  for(j in a[-i])
  {
    largos<-append(largos, length(grupos$patrilineal2[[j]])) #vector de conteo mujeres por grupo
  } #cierre bucle j (ojo, hay otro abajo) 
  
  largos<-largos/sum(largos) #conversion de los largos a proporciones del total
  sizes<-round(largos*length(grupos$patrilineal2[[i]]))
  sizes<-replace(sizes,sizes==0,1) #mecanismo de seguridad por si algun elemento tiene valor cero
  # (el valor cero genera un vector reparto con un elemento menos)  

  set.seed(semillas[i]) #fijado del seed para el barajado asi mujeres y sus haplos se barajan en el mismo orden
	# (la maravilla de los seudoaleatorios)
  grupos$patrilineal2[[i]]<-sample(grupos$patrilineal2[[i]]) #barajado de madres
  set.seed(semillas[i]) 
  gruposhijxs$patrilineal2[[i]]<-sample(gruposhijxs$patrilineal2[[i]]) #barajado de numeros de progenie
  suppressWarnings(reparto<-split(grupos$patrilineal2[[i]],rep(seq_along(sizes), sizes)))
  #division en 5 partes, proporcionales al grupo de llegada, de las madres
  #suprimo los warnings temporalmente porque a veces el reparto no suma igual al largo total
  suppressWarnings(reparto2<-split(gruposhijxs$patrilineal2[[i]],rep(seq_along(sizes), sizes)))
  #division en 5 partes, proporcionales al grupo de llegada, de los numeros de progenie
  
  k<-1
   for(j in a[-i])
  {
    preasignacion$madres[[j]]<-
      append(preasignacion$madres[[j]],reparto[[k]]) #reparto de madres
	preasignacion$nrohijos[[j]]<-
      append(preasignacion$nrohijos[[j]],reparto2[[k]]) #reparto del n de hijos de cada madre
	k<-k+1
  } # cierre bucle j
  
  rm(.Random.seed, envir=globalenv()) #para quitar la semilla prefijada
  #(que puede sesgar la aleatorizacion posterior)

} # cierre bucle i

# AHORA SI: Asignacion de haplos a las progenies

for(i in 1:6)
{ 
  generacion$patrilineal2[[i]]<-vector()
  # se multiplica cada haplo por cada vector de hijes
  # por ejemplo, si preasignacion$madres[[1]] es c("A","B","C","D")
  # y si preasignacion$nrohijos[[1]] es c(2,3,2,3)
  # generacion[[1]] es c("A","A","B","B","B","C","C","D","D","D")
  for(j in 1:length(preasignacion$madres[[i]]))
  {
    progenie<-rep(preasignacion$madres[[i]][[j]],preasignacion$nrohijos[[i]][j]) 
    generacion$patrilineal2[[i]]<-append(generacion$patrilineal2[[i]],progenie)
  } #cierre del bucle j  

  # se asigna sexo M y F al azar a les integrantes de cada generacion
  # y se divide cada generacion en subcomponentes $varon y $mujer en forma acorde
muestra<-sample(c("varon","mujer"),
                length(generacion$patrilineal2[[i]]), replace = T)
generacion$patrilineal2[[i]]<-split(generacion$patrilineal2[[i]],muestra)
} #cierre del bucle i

# Pasaje de la progenie a otro vector de varones o mujeres
# segun el patron de residencia en un vector "cementerio"

# Y en este caso van todos para el cementerio!!

  for(i in 1:6)
  {
    cementerios$patrilineal2[[i]]<-generacion$patrilineal2[[i]]
  }
