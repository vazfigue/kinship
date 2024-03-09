#########

## Seteo de los grandes parametros para las simulaciones
## Van aca arriba porque son faciles de encontrar

## FERTILIDAD Y MORTALIDAD JUVENIL

# tfr: tasa de fertilidad total (1=alta o 2=baja)
tfr<-1 #sustituir por 2 si se quiere que el modelo corra fertilidad baja
fertility<-matrix(c(6.1,1.1,5.5,1.2), 2, 2, byrow = T)

# cm: mortalidad juvenil (1=alta o 2=baja)
cm<-1 #sustituir por 2 si se quiere que el modelo corra mortalidad baja
mortality<-matrix(c(0.434,0.111,0.381,0.104), 2, 2, byrow = T) 

## IMPORTANTE:
## Dada la formulacion de la simulacion, sin mecanismos de seguridad
## en caso de que una poblacion alcance poblacion cero,
## los escenarios de alta tfr + alta mortalidad y baja tfr + baja mortalidad
## colapsan en ocasiones. Paciencia.
## El escenario de baja tfr + alta mortalidad colapsa repetidamente
## y no es viable.

# fs: supervivencia femenina en edad reproductiva
fs<-0.7 #se mantiene fijo para todas las simulaciones, pero puede cambiarse a voluntad
# 0.75 es de Ache (Hill y Hurtado)
# Antes probe con 0.46 de Semai (Fix) y crasheaba por extincion de la poblacion

# TENER EN CUENTA QUE BASELINE VS MOD SE DEFINE PARA CADA CORRIDA
# (ABAJO)

## NOMBRES DE ARCHIVO PARA EL GUARDADO DE LOS RESULTADOS 

# Definicion de una matriz de nombres para el guardado
# de los resultados de las simulaciones baseline
# El nombre de archivo a guardar va a depender de los valores de tfr y cm
# Descomentar para las simulaciones baseline, comentar si se activan mr y ad
files<-matrix(c("HF_HM_base_","HF_LM_base_",
                "LF_HM_base_","LF_LM_base_"), 2, 2, byrow = T)

# Definicion de una matriz de nombres para el guardado
# de los resultados de las simulaciones con migracion y (no) adhesion
# El nombre de archivo a guardar va a depender de los valores de tfr y cm
# Comentar para las simulaciones baseline, descomentar si se activan mr y ad
# files<-matrix(c("HF_HM_mod_","HF_LM_mod_",
                # "LF_HM_mod_","LF_LM_mod_"), 2, 2, byrow = T)

#########

# Inicializacion del cronometro
start.time <- Sys.time()

# Generacion de todas las estructuras de datos que se van a emplear en las subrutinas

# Funcion in da house para calcular diversidad haplotipica
# (no saco de un paquete porque las formulas, p.ej. de Pegas usan secuencias de ADN)
hapdiv <- function(mivector){
  frecabs<-table(mivector)
  resultado<-(sum(frecabs)/(sum(frecabs)-1))*(1-sum((frecabs/sum(frecabs))^2))  
  return(resultado)
}

# Etiquetas de las listas que se integraran a listas de listas
tags<-c("matrilineal","bilateral","patrilineal1","patrilineal2")

# Lista de listas de grupos de parentesco y su progenie
grupos<-vector("list", length=4)
names(grupos)<-tags
gruposhijxs<-vector("list", length=4)
names(gruposhijxs)<-tags

# Lista de individuos en la nueva generacion
generacion<-vector("list", length=4)
names(generacion)<-tags

# Lista de vectores de cementerios de grupos de parentesco, con subgrupos varon y mujer
cementerios<-vector("list", length=4)
names(cementerios)<-tags
for(i in 1:4)
{
	for (j in 1:6)
	{
	cementerios[[i]][[j]]<-list()
	}
}

# Lista de generaciones de cementerios
piladehuesos<-vector("list", length=4)
names(piladehuesos)<-tags

# Lista de matrices de datos para TODOS los datos de TODAS las corridas
# y de los remuestreos de 10, 20, 30 y 50 individuos
data<-vector("list", length=4)
names(data)<-tags
resamples<-vector("list", length=4)
names(resamples)<-tags

for(i in 1:4)
{
 data[[i]]<-data.frame(matrix(nrow = 0,ncol = 36))
 resamples[[i]]<-data.frame(matrix(nrow = 0,ncol = 88))
}

# etiquetas de data (se pegan mas abajo)
# run: numero de corrida
# gen: numero de generacion
# nf: numero de mujeres inicial iniciando la corrida
# tfr: tasa de fertilidad total (alta o baja)
# cm: mortalidad juvenil (alta o baja)
# mr: tasa de migracion
# ad: (no) adhesion a la pauta de residencia
# fnf: numero de mujeres final en la generacion
# divM: diversidad de haplotipos masculina
# divF: diversidad de haplotipos femenina
# divT: diversidad de haplotipos total

# Lista de dataframes con distancias de Nei entre grupos
# de la generacion 15 de cada corrida
distances<-vector("list", length=4)
names(distances)<-tags
for(i in 1:4)
{
distances[[i]]<-list()	
	for(j in 1:3)
	{
	distances[[i]][[j]]<-data.frame(matrix(nrow = 0,ncol = 15))
	}
}

# Vector de 24 haplos "internos" a la poblacion (A,B,C,D, etc)
haplos_int<-LETTERS[seq(1,24)]

# Vector de otros 24 haplos "externos" a la poblacion (AA, AB, AC, AD, etc)
haplos_ext<-paste0("A",LETTERS[seq(1,24)])

# Vector de todos los haplos para la tabla de frecuencias de haplos
haplostodos<-c(haplos_int,haplos_ext)

# Inicializacion de una barra de avance
pb <- txtProgressBar(min = 1, max = 100, style = 3)

for(run in 1:500) #setear i en 1 si se quiere hacer una sola corrida o debugging
{ # apertura del loop run

  # Definicion de parametros para cada corrida
  
  # nf: numero de mujeres inicial en cada grupo
  nf<-sample(20:65,6,replace=T) #el minimo era de 12 (de vuelta, Semai) pero crasheaba la simulacion

  # mr: tasa de migracion (igual para todos los grupos en una corrida)
  # mr<-abs(rnorm(1, 0.189, 0.092))# comentar para las simulaciones "baseline"
  mr<-0 #commentar para las simulaciones con migracion
  
  # ad: (no) adhesion a la pauta de residencia (igual para todos los grupos en una corrida)
  # ad<-abs(rnorm(1, 0.2, 0.1))# comentar para las simulaciones "baseline"
  ad<-0 #comentar para las simulaciones con no adhesion a las reglas

	# GENERACION 1
	gen<-1
	source("01_simulacion_matri_gen1_v3.R")
	source("02_simulacion_cognat_gen1_v3.R")
	source("03_simulacion_patri1_gen1_v3.R")
	source("04_simulacion_patri2_gen1_v3.R")	

# modificacion de los cementerios en funcion de la no adhesion a las pautas de residencia
a<-1:4
for(j in 1:4)
{ # apertura del loop j
extkin<-sample(a[-j],6, replace=T) #muestra de 6 valores enteros entre 1 y 4 excluyendo al valor indicado por j 
	for(i in 1:6)
		{		
		nonstrictM<-sample(cementerios[[extkin[i]]][[i]]$varon,
			round(length(cementerios[[j]][[i]]$varon)*ad),replace=T)	
		nonstrictF<-sample(cementerios[[extkin[i]]][[i]]$mujer,
			round(length(cementerios[[j]][[i]]$mujer)*ad),replace=T)
		strictM<-sample(cementerios[[j]][[i]]$varon,
			length(cementerios[[j]][[i]]$varon)-length(nonstrictM),replace=F)
		strictF<-sample(cementerios[[j]][[i]]$mujer,
			length(cementerios[[j]][[i]]$mujer)-length(nonstrictF),replace=F)
		cementerios[[j]][[i]]$varon<-append(strictM, nonstrictM)
		cementerios[[j]][[i]]$mujer<-append(strictF, nonstrictF)
		}
# guardado del cementerio en piladehuesos
piladehuesos[[j]][[gen]]<-cementerios[[j]]

# calculos de diversidad haplotipica miscelanea
    hmasc<-vector(length=6)
    hfem<-vector(length=6)
    htotal<-vector(length=6)
    
for(i in 1:6)
  {
    hmasc[i]<-hapdiv(cementerios[[j]][[i]]$varon)
    hfem[i]<-hapdiv(cementerios[[j]][[i]]$mujer)
    htotal[i]<-hapdiv(unlist(cementerios[[j]][[i]]))
  }

# guardado de todos los resultados de la simulacion
# en una fila del dataframe de resultados
data[[j]]<-rbind(data[[j]], c(run,gen,nf,tfr,cm,mr,ad,
							length(cementerios[[j]][[1]]$mujer),length(cementerios[[j]][[2]]$mujer),
							length(cementerios[[j]][[3]]$mujer),length(cementerios[[j]][[4]]$mujer),
							length(cementerios[[j]][[5]]$mujer),length(cementerios[[j]][[6]]$mujer),
							hmasc,hfem,htotal))

names(data[[j]])<-c("run","gen","nf1","nf2","nf3","nf4","nf5","nf6","tfr","cm","mr",
                     "ad","fnf1","fnf2","fnf3","fnf4","fnf5","fnf6","divm1","divm2",
                     "divm3","divm4","divm5","divm6","divf1","divf2","divf3","divf4",
                     "divf5","divf6","divt1","divt2","divt3","divt4","divt5","divt6")
} # cierre del loop j

####aca van las generaciones 2 a 15
for(gen in 2:15)
{
	source("01_simulacion_matri_gen2_15_v3.R")
	source("02_simulacion_cognat_gen2_15_v3.R")
	source("03_simulacion_patri1_gen2_15_v3.R")
	source("04_simulacion_patri2_gen2_15_v3.R")
	
  # modificacion de los cementerios en funcion de la no adhesion a las pautas de residencia
	a<-1:4
	for(j in 1:4)
	{ # apertura del loop j
	extkin<-sample(a[-j],6, replace=T) #muestra de 6 valores enteros entre 1 y 4 excluyendo al valor indicado por j 
		for(i in 1:6)
		{		
		nonstrictM<-sample(cementerios[[extkin[i]]][[i]]$varon,
			round(length(cementerios[[j]][[i]]$varon)*ad),replace=T)	
		nonstrictF<-sample(cementerios[[extkin[i]]][[i]]$mujer,
			round(length(cementerios[[j]][[i]]$mujer)*ad),replace=T)
		strictM<-sample(cementerios[[j]][[i]]$varon,
			abs(length(cementerios[[j]][[i]]$varon)-length(nonstrictM)),replace=F)
		strictF<-sample(cementerios[[j]][[i]]$mujer,
			abs(length(cementerios[[j]][[i]]$mujer)-length(nonstrictF)),replace=F)
		cementerios[[j]][[i]]$varon<-append(strictM, nonstrictM)
		cementerios[[j]][[i]]$mujer<-append(strictF, nonstrictF)
		}
	# guardado del cementerio en piladehuesos
	piladehuesos[[j]][[gen]]<-cementerios[[j]]

	# calculos de diversidad haplotipica miscelanea
	hmasc<-vector(length=6)
	hfem<-vector(length=6)
	htotal<-vector(length=6)

	for(i in 1:6)
	{
	hmasc[i]<-hapdiv(cementerios[[j]][[i]]$varon)
	hfem[i]<-hapdiv(cementerios[[j]][[i]]$mujer)
	htotal[i]<-hapdiv(unlist(cementerios[[j]][[i]]))
	}
	
  # guardado de todos los resultados de la simulacion
  # en una fila del dataframe de resultados
	data[[j]]<-rbind(data[[j]], c(run,gen,nf,tfr,cm,mr,ad,
							length(cementerios[[j]][[1]]$mujer),length(cementerios[[j]][[2]]$mujer),
							length(cementerios[[j]][[3]]$mujer),length(cementerios[[j]][[4]]$mujer),
							length(cementerios[[j]][[5]]$mujer),length(cementerios[[j]][[6]]$mujer),
							hmasc,hfem,htotal))
	} # cierre del loop j

} # cierre del loop gen

###ACA VAN LOS CALCULOS DE DIVERSIDAD DE LOS REMUESTREOS

for(j in 1:4)
{ # apertura del loop j

  # calculos de diversidad haplotipica miscelanea de muestras (5M 5F, 10M 10F, 15M 15F, 25M 25F)
  hmasc10<-vector(length=6)
  hfem10<-vector(length=6)
  htotal10<-vector(length=6)
  hmasc20<-vector(length=6)
  hfem20<-vector(length=6)
  htotal20<-vector(length=6)
  hmasc30<-vector(length=6)
  hfem30<-vector(length=6)
  htotal30<-vector(length=6)
  hmasc50<-vector(length=6)
  hfem50<-vector(length=6)
  htotal50<-vector(length=6)
  
  for(i in 1:6)
  {
    if (length(cementerios[[j]][[i]]$varon)<5){
      hmasc10[i]<-NA
    } else {
      hmasc10[i]<-hapdiv(sample(cementerios[[j]][[i]]$varon,5))
    }
    if (length(cementerios[[j]][[i]]$mujer)<5){
      hfem10[i]<-NA
    } else {
      hfem10[i]<-hapdiv(sample(cementerios[[j]][[i]]$mujer,5))
    }
    if (length(unlist(cementerios[[j]][[i]]))<10){
      htotal10[i]<-NA
    } else {
      htotal10[i]<-hapdiv(sample(unlist(cementerios[[j]][[i]]),10))
    }
    
    if (length(cementerios[[j]][[i]]$varon)<10){
      hmasc20[i]<-NA
    } else {
      hmasc20[i]<-hapdiv(sample(cementerios[[j]][[i]]$varon,10))
    }
    if (length(cementerios[[j]][[i]]$mujer)<10){
      hfem20[i]<-NA
    } else {
      hfem20[i]<-hapdiv(sample(cementerios[[j]][[i]]$mujer,10))
    }
    if (length(unlist(cementerios[[j]][[i]]))<20){
      htotal20[i]<-NA
    } else {
      htotal20[i]<-hapdiv(sample(unlist(cementerios[[j]][[i]]),20))
    }
    
    if (length(cementerios[[j]][[i]]$varon)<15){
      hmasc30[i]<-NA
    } else {
      hmasc30[i]<-hapdiv(sample(cementerios[[j]][[i]]$varon,15))
    }
    if (length(cementerios[[j]][[i]]$mujer)<15){
      hfem30[i]<-NA
    } else {
      hfem30[i]<-hapdiv(sample(cementerios[[j]][[i]]$mujer,15))
    }
    if (length(unlist(cementerios[[j]][[i]]))<30){
      htotal30[i]<-NA
    } else {
      htotal30[i]<-hapdiv(sample(unlist(cementerios[[j]][[i]]),30))
    }
    
    if (length(cementerios[[j]][[i]]$varon)<25){
      hmasc50[i]<-NA
    } else {
      hmasc50[i]<-hapdiv(sample(cementerios[[j]][[i]]$varon,25))
    }
    if (length(cementerios[[j]][[i]]$mujer)<25){
      hfem50[i]<-NA
    } else {
      hfem50[i]<-hapdiv(sample(cementerios[[j]][[i]]$mujer,25))
    }
    if (length(unlist(cementerios[[j]][[i]]))<50){
      htotal50[i]<-NA
    } else {
      htotal50[i]<-hapdiv(sample(unlist(cementerios[[j]][[i]]),50))
    }
  }
  
# guardado de los resultados de los remuestreos
# en una fila de otro dataframe 
resamples[[j]]<-rbind(resamples[[j]], c(run,gen,mr,ad,hmasc10,mean(hmasc10),hfem10,mean(hfem10),htotal10,mean(htotal10),
                                        hmasc20,mean(hmasc20),hfem20,mean(hfem20),htotal20,mean(htotal20),
                                        hmasc30,mean(hmasc30),hfem30,mean(hfem30),htotal30,mean(htotal30),
                                        hmasc50,mean(hmasc50),hfem50,mean(hfem50),htotal50,mean(htotal50)))

names(resamples[[j]])<-c("run","gen","mr","ad","divM1_10","divM2_10","divM3_10","divM4_10","divM5_10","divM6_10","divM10_mean",
                         "divF1_10","divF2_10","divF3_10","divF4_10","divF5_10","divF6_10","divF10_mean",
                         "divT1_10","divT2_10","divT3_10","divT4_10","divT5_10","divT6_10","divT10_mean",
                         "divM1_20","divM2_20","divM3_20","divM4_20","divM5_20","divM6_20","divm20_mean",
                         "divF1_20","divF2_20","divF3_20","divF4_20","divF5_20","divF6_20","divF20_mean",
                         "divT1_20","divT2_20","divT3_20","divT4_20","divT5_20","divT6_20","divT20_mean",
                         "divM1_30","divM2_30","divM3_30","divM4_30","divM5_30","divM6_30","divM30_mean",
                         "divF1_30","divF2_30","divF3_30","divF4_30","divF5_30","divF6_30","divF30_mean",
                         "divT1_30","divT2_30","divT3_30","divT4_30","divT5_30","divT6_30","divT30_mean",
                         "divM1_50","divM2_50","divM3_50","divM4_50","divM5_50","divM6_50","divM50_mean",
                         "divF1_50","divF2_50","divF3_50","divF4_50","divF5_50","divF6_50","divF50_mean",
                         "divT1_50","divT2_50","divT3_50","divT4_50","divT5_50","divT6_50","divT50_mean")
} # cierre del loop j

# Suspension del sistema por 0,01s
  Sys.sleep(0.01)
  
  # Detalle de progreso
  setTxtProgressBar(pb, run/5)

} # cierre del loop run

# Cierre de la barra
close(pb)

# Guardado de los archivos de datos
for(i in 1:4)
{
filename<-paste(files[tfr,cm],tags[i],".csv",sep="")
filename2<-paste(files[tfr,cm],tags[i],"resamp.csv",sep="")
write.csv2(data[[i]],file=filename,row.names = F)
write.csv2(resamples[[i]],file=filename2,row.names = F)
# write.csv2(distances[[i]][[1]],file=paste("distM_",filename,sep=""),row.names = F)
# write.csv2(distances[[i]][[2]],file=paste("distF_",filename,sep=""),row.names = F)
# write.csv2(distances[[i]][[3]],file=paste("distT_",filename,sep=""),row.names = F)
}

# Guardado de los .RData por las dudas y para levantar en otros analisis
datafilename<-paste(files[tfr,cm],".RData",sep="")
save.image(datafilename)

# Finalizacion del cronometro
end.time <- Sys.time()

time.taken <- round(end.time - start.time,2)

# TIEMPOS (Intel Core i3, 3.40 GHz, 16 GB RAM):
# Baseline alta fertilidad alta mortalidad: 22,58 min
# Baseline alta fertilidad baja mortalidad: 30,4 min
# Baseline baja fertilidad alta mortalidad: no se simulo, las poblaciones crasheaban repetidamente
# Baseline baja fertilidad baja mortalidad: 20,5 min
