#################################################################
#FUNCION PARA OMPTIMIZACION INCLUYENDO NIEVE ABCD

#EJECUCION:

#"c:\Program Files\R\R-4.2.1\bin\Rscript.exe" ABCDmod.R 2 NSE 1 data3045.txt ABCDparam.txt

#Los argumentos son:
#ABCDmod.R: nombre del archivo de código
#2: Modo de ejecución (0 1 2 3)
#NSE: Función objetivo (MSE RMSE NSE R2 KGE)
#1: Si queremos guardar el modelo (0 o 1)
#data3045.txt: Archivo de datos
#ABCDparam.txt: Archivo de parámetros
#

#variables de estado que dejamos en 0 en el primer mes de la simulacion
#Sw0=0;Sg0=0;Spt0=0

#parámetros del modelo (estos sería por defecto, y os dejo las orquillas):

#abcd.param=c(0.98,350,0.3,0.5,0.5,6,-0.95)
#abcd.lower=c(0.00001,10,0,0.00001,0.00001,2,-6)
#abcd.upper=c(1,2000,1,100,1,11,1)

#c("a","b","c","d","Xg","a1","a2")

#Recordad, los cuatro primeros son del abcd, Xg el de intercambio y los dos últimos las temperaturas del modelo nival

#modifica controla el tipo de modelo:

#0: abcd original de 4 parámetros (a b c d)
#1: abcd y sub. (los cuatro + Xg)
#2: 4 de abcd + Xg + a1 + a2
#3: abcd y nival: parámetros a b c d + a1 + a2

#Los datos de entrada son dos archivos de texto:

#a) Parámetros del modelo


#b) Datos


#Y los parámetros textuales de ejecución son:


#################################################################


#0. librerias:##############################################################

require(readr)
require(hydroGOF)

#1. Parámetros propios del programa#########################################

args <- commandArgs(trailingOnly = TRUE)
print(args)

Sw0=0
Sg0=0
Spt0=0

mod = args[1]
f0=args[2]
save.data = args[3]

dataname = args[4]
paramname = args[5]


#1. lectura de los datos de entrada (archivos txt)###########################

#1.1. Lectura de los parámetros del modelo (son los que irán variando)

#parameters
param = read_csv(paramname, col_names = FALSE)

#datos
data <- read_csv("data3045.txt")
data$date = as.Date(paste0(data$date,"01"),format = "%Y%m%d")
data = data[order(data$date),]

#2. Ejecución del modelo#####################################################

#Para modificaciones con nival (modos 2 o 3)
if(mod >= 2){
  data[,"snt"]<-NA
  data[,"rt"]<-NA
  data[,"mt"]<-NA
  data[,"spt"]<-NA    
}
#y resto
data[,"Wi"]<-NA 
data[,"Yi"]<-NA
data[,"Swi"]<-NA
data[,"eti"]<-NA  
data[,"Rgi"]<-NA    
data[,"Roi"]<-NA    
data[,"Sgi"]<-NA      
data[,"Qgi"]<-NA
data[,"Qbi"]<-NA
data[,"Qmmsim"]<-NA


for(i in 1:nrow(data)){

  if(i == 1){
    Swi = Sw0
    Sgi = Sg0
    Spti = Spt0
  }
  
  #A) RUTINA WASMOD (modos 2 o 3)
  if(mod >= 2){
    #nieve
    data[i,"snt"] <- data[i,"pre"] * (1 - exp((data[i,"tas"] - param[6]) / (param[6] - param[7]))^2)
    data[i,"snt"] <- max(data[i,"snt"],0)
    #lluvia
    data[i,"rt"] <- data[i,"pre"] - data[i,"snt"]
    #deshielo
    data[i,"mt"] <- Spti * (1 - exp((-(data[i,"tas"] - param[7])) / (param[6] - param[7]))^2)
    data[i,"mt"] <- max(data[i,"mt"],0)
    #almacenamiento
    data[i,"spt"] <- Spti + data[i,"snt"] - data[i,"mt"]
    #Y voy a agregar a la precipitacion el deshielo
    liquido <- data[i,"rt"] + data[i,"mt"]
  }
  
  #B) PRIMER EMBALSE ABCD
  if(mod >= 2){
    data[i,"Wi"] <- liquido + Swi
  }else{
    data[i,"Wi"] <- Swi + data[i,"pre"]
  }
  
  #Componente Yi
  data[i,"Yi"] <- (data[i,"Wi"] + param[2])/(2 * param[1]) - (((data[i,"Wi"] + param[2]) / (2 * param[1]))^2 - ((data[i,"Wi"] * param[2]) / param[1]))^0.5
  #Componente Swi
  data[i,"Swi"] <- data[i,"Yi"] * exp(- data[i,"etp"] / param[2])
  data[i,"eti"] <- data[i,"Yi"] - data[i,"Swi"] # ETi
  #descomp salida
  data[i,"Rgi"] <- param[3] * (data[i,"Wi"] - data[i,"Yi"])
  data[i,"Roi"] <- (1 - param[3]) * (data[i,"Wi"] - data[i,"Yi"])
  
  #C) SEGUNDO EMBALSE
  data[i,"Sgi"] <- (data[i,"Rgi"] + Sgi) / (1 + param[4])
  data[i,"Qgi"] <- param[4] * data[i,"Sgi"]
  #mod si utilizamos Xg o Xg+Wasmod (modos 1 o 2) o el resto:
  if(mod == 0 || mod == 3){
    data[i,"Qbi"] <- data[i,"Qgi"]
  }else{
    data[i,"Qbi"] <- param[5] * data[i,"Qgi"]
  }
  
  #suma descarga y esc
  data[i,"Qmmsim"] <- data[i,"Roi"] + data[i,"Qbi"]
  
  #y asignamos nuevos valores para el resto del bucle
  Sgi <-  data[i,"Sgi"]
  Swi <- data[i,"Swi"]
  if(mod >= 2){
    Spti <- data[i,"spt"]
  }

}

if(save.data == 1){
  write_csv(data,file="ABCDresults.txt")
}

#3. Estadisticos de bondad de ajuste en calibración############

#f0 que vamos a establecer como el negativo de NSE (para minimizar)
adjust = suppressWarnings(gof(sim=data$Qmmsim,obs=data$qmmobs,digits=6)[c(3,4,6,9,17,19),1])

if(f0 %in% c("NSE","R2","KGE")){
  adjust[["f0"]] = - adjust[[f0]]
}else{
  adjust[["f0"]] = adjust[[f0]]
}

write_csv(data.frame(t(data.frame(adjust))),file="ABCDadjust.txt")
