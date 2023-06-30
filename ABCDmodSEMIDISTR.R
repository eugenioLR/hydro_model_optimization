#################################################################
#FUNCION PARA OMPTIMIZACION INCLUYENDO NIEVE ABCD

#OK, OS HAGO ESTE MAIN PORQUE EL PRIMER ARCHIVO ERA LA FUNCIÓN DEL MODELO, PERO AHORA SE DEBE EJECUTAR CONSIDERANDO QUE SON VARIAS SUBCUENCAS Y QUE LAS QUE HAY AGUAS ABAJO CONTIENEN COMO ENTRADA LAS DE AGUAS ARRIBA

#este archivo es una modificación del que os pasé, con los mismo parámetros más dos adicionales:

#archivo de características de las subcuencas
#1 si queremos calibrar y 0 si no

#Además, el archivo de parametros ahora tiene tres filas: valor por omisión y orquilla

#EJECUCION:

#"c:\Program Files\R\R-4.2.1\bin\Rscript.exe" ABCDmodSEMIDISTR.R 2 NSE 1 0 data.txt ABCDparamSEMIDISTR.txt basins.txt

#Los argumentos son:
#ABCDmodSEMIDISTR.R: nombre de este archivo de código
#2: Modo de ejecución (0 1 2 3)
#NSE: Función objetivo (MSE RMSE NSE R2 KGE)
#1: Si queremos guardar el modelo (0 o 1)
#0: SI queremos calibrar el modelo (0 no y 1 si)
#data.txt: Archivo de datos (en él aparecen todas las subcuencas)
#ABCDparam.txt: Archivo de parámetros. Ahora incluyo los valores centrales, y la segunda fila es la orquilla inferior y la tercera la superior
#basins.txt: caracteristicas de las cuencas
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

#SOBRE EL ARCHIVO DE SUBCUENCAS:
#basins.txt tiene 4 columnas:
#code: código de la subcuenca
#order: Orden de la subcuenca: las de orden 1 son las que cabecera (veréis más adelante que primero se ejecutan los modelos en éstas, después las de orden 2, y así sucesivamente)
#codedown: Importante, es el código de la cuenca que tiene cada una de ellas aguas abajo; es la columna que define la topología del esquema de cuenca
#supha: Necesario porque ahora trabajaremos con hm3 para el caudal y no mm/m2 como antes

#################################################################

#0. librerias:###################################################

require(readr)#lectura de datos

require(hydroGOF)#bondad de ajuste
require(rtop)#optimizacion SCE-UA

require(dplyr)

#0. funcion ABCD#################################################

#Lo que antes era el archivo que os pasaba lo he definido como función para ir llamándola de forma iterativa. Tiene los mismos parámetros que os dije en su momento de modo, datos, etc, más el parámero de superficie. El ndata son los datos adicionales que se agregarían si fuese una subcuenca con entradas de subcuencas aguas arriba

fun.abcd = function(mod,subdata,param,supha,subdata2=NA){
  
  #Variables de estado:
  Sw0=0
  Sg0=0
  Spt0=0
  #Para modificaciones con nival (modos 2 o 3)
  if(mod >= 2){
    subdata[,"snt"]<-NA
    subdata[,"rt"]<-NA
    subdata[,"mt"]<-NA
    subdata[,"spt"]<-NA    
  }
  #y resto
  subdata[,"Wi"]<-NA 
  subdata[,"Yi"]<-NA
  subdata[,"Swi"]<-NA
  subdata[,"eti"]<-NA  
  subdata[,"Rgi"]<-NA    
  subdata[,"Roi"]<-NA    
  subdata[,"Sgi"]<-NA      
  subdata[,"Qgi"]<-NA
  subdata[,"Qbi"]<-NA
  subdata[,"Qmmsim"]<-NA
  subdata[,"Qhmsim"]<-NA
  
  #iniciemos el bucle para cada mes:
  
  for(i in 1:nrow(subdata)){
    if(i == 1){
      Swi = Sw0
      Sgi = Sg0
      Spti = Spt0
    }
    #A) RUTINA WASMOD (modos 2 o 3)
    if(mod >= 2){
      #nieve
      subdata[i,"snt"] <- subdata[i,"pre"] * (1 - exp((subdata[i,"tas"] - param[7]) / (param[6] - param[7]))^2)
      subdata[i,"snt"] <- max(subdata[i,"snt"],0)
      #lluvia
      subdata[i,"rt"] <- subdata[i,"pre"] - subdata[i,"snt"]
      #deshielo
      subdata[i,"mt"] <- Spti * (1 - exp((-(subdata[i,"tas"] - param[7])) / (param[6] - param[7]))^2)
      subdata[i,"mt"] <- max(subdata[i,"mt"],0)
      #almacenamiento
      subdata[i,"spt"] <- Spti + subdata[i,"snt"] - subdata[i,"mt"]
      #Y voy a agregar a la precipitacion el deshielo
      liquido <- subdata[i,"rt"] + subdata[i,"mt"]
    }
    #B) PRIMER EMBALSE ABCD
    if(mod >= 2){
      subdata[i,"Wi"] <- liquido + Swi
    }else{
      subdata[i,"Wi"] <- Swi + subdata[i,"pre"]
    }
    #Componente Yi
    subdata[i,"Yi"] <- (subdata[i,"Wi"] + param[2])/(2 * param[1]) - (((subdata[i,"Wi"] + param[2]) / (2 * param[1]))^2 - ((subdata[i,"Wi"] * param[2]) / param[1]))^0.5
    #Componente Swi
    subdata[i,"Swi"] <- subdata[i,"Yi"] * exp(- subdata[i,"etp"] / param[2])
    subdata[i,"eti"] <- subdata[i,"Yi"] - subdata[i,"Swi"] # ETi
    #descomp salida
    subdata[i,"Rgi"] <- param[3] * (subdata[i,"Wi"] - subdata[i,"Yi"])
    subdata[i,"Roi"] <- (1 - param[3]) * (subdata[i,"Wi"] - subdata[i,"Yi"])
    
    #C) SEGUNDO EMBALSE
    subdata[i,"Sgi"] <- (subdata[i,"Rgi"] + Sgi) / (1 + param[4])
    subdata[i,"Qgi"] <- param[4] * subdata[i,"Sgi"]

    #mod si utilizamos Xg o Xg+Wasmod (modos 1 o 2) o el resto:
    if(mod == 0 || mod == 3){
      subdata[i,"Qbi"] <- subdata[i,"Qgi"]
    }else{
      subdata[i,"Qbi"] <- param[5] * subdata[i,"Qgi"]
    }
    #suma descarga y esc
    subdata[i,"Qmmsim"] <- subdata[i,"Roi"] + subdata[i,"Qbi"]
    #y asignamos nuevos valores para el resto del bucle
    Sgi <-  subdata[i,"Sgi"]
    Swi <- subdata[i,"Swi"]
    if(mod >= 2){
      Spti <- subdata[i,"spt"]
    }
  }
  #Ahora calculamos la aportación en hm3/mes:
  subdata[,"Qhmsim"] = subdata[,"Qmmsim"] * supha / 100000
  
  #Finalmente, si existen subcuencas aguas arriba se agrega al Q obtenido; esta parte procede del cogido dentro del que se ejecuta esta funcion
  if(length(subdata2) > 1){
    subdata$Qhmsim = subdata$Qhmsim + subdata2
  }
  return(subdata)
}  

#1. Parámetros propios del programa#########################################

args <- commandArgs(trailingOnly = TRUE)
print(args)

mod = args[1]
f0=args[2]
save.data = args[3]
runopt = args[4]
dataname = args[5]
paramname= args[6]
basinsname= args[7]

#mod = 2
#f0 = "NSE"
#save.data = 1
#runopt = 1
#dataname = "data.txt"
#paramname = "ABCDparamSEMIDISTR.txt"
#basinsname = "basins.txt"


#2. lectura de los datos de entrada (archivos txt)###########################

#parameters
param = as.data.frame(read_csv(paramname, col_names = FALSE))
#datos
data <- as.data.frame(read_csv(dataname))
data$date = as.Date(paste0(data$date,"01"),format = "%Y%m%d")
data = data[order(data$code,data$date),]

basins = as.data.frame(read_csv(basinsname))
#y le vamos a añadir a basins columnas para que nos incluya en ellos los estadisticos de bondad de ajuste (recordad, calcularemos para cada subcuenca y la última en realidad será la bondad de ajuste general)
#De paso le voy a añadir otra columna para almacenar los parámetros utilizados en cada subcuenca
basins = data.frame(basins,MSE=NA,RMSE=NA,PBIAS=NA,NSE=NA,R2=NA,KGE=NA,f0=NA,param=NA)

#3. OK COMENCEMOS CON LOS MODELOS###########################

#Bucle por subcuenca. Para asegurarnos de que se ejecutan para los diferentes ordenes:
basins = basins[order(basins$order),]

#en mi versión creo dobles bucles para poder paralelizar el proceso pero aquí mejor en secuencial y si tenéis que cambiarlo no es mayor problema

subbasins = basins$code

#bucle por subcuenca:
for(i in 1:length(subbasins)){
    
  subb = subbasins[i]
  subdata = data[data[,1] == subb,]
  supha = basins[basins[,1] == subb, 4]
    
  #este check es para saber si tenemos datos aguas arriba o no. Si existe creamos el vector a sumar en la función
  subb2=basins[basins$codedown == subb,]

  
  #
  if(nrow(subb2)>0){
    subdata2 = data.mod[data.mod$code %in% subb2$code,]
    subdata2 = aggregate(Qhmsim ~ date, data = subdata2, FUN = sum)[,"Qhmsim"]
  }else{
    subdata2 = NA
  }

  #ejecutemos el modelo. Ok, esta función es por si quiero optimizar (runopt=1) o no (runopt=0). En vuestro caso no utilizar el código del else:
  
  if(runopt == 0){
    subdata.mod = fun.abcd(mod = mod, subdata = subdata, param = as.numeric(param[1,]), supha = supha, subdata2 = subdata2)
    param.new = as.numeric(param[1,])
  }else if (runopt == 1){
    #OJO ESTE CODIGO ES PARA MI OPTIMIZACIÓN SCE-UA!!
    
    #voy a acotar parametros en funcion del mod:
    if(mod == 0){
      param2 = param[,1:4]
    }else if (mod == 1){
      param2 = param[,1:5]
    }else if (mod == 3){
      param2 = param[,c(1:4,6:7)]
    }else{
      param2 = param
    }

    abcd.run <- function(x) {
      if(mod == 0){
        pars<- c(x,-99,-99,-99)
      }else if (mod == 1){
        param2 = c(x,-1,-99)
      }else if (mod == 3){
        param2 = c(x[1:4],-99,x[5:6])
      }else{
        pars = x
      }
      modelico=fun.abcd(mod = mod, subdata = subdata, param = pars, supha = supha, subdata2 = subdata2)
      return(modelico)
    } 
    abcd.F0.nse <- function (x) {
      modelico=abcd.run(x)
      #F0 NSE:
      - NSeff(sim=modelico[,"Qhmsim"],obs=modelico[,"qhmobs"])
    }
    #optimiza VERSION sceua
    print(system.time(abcd.fit.nse <- sceua(OFUN = abcd.F0.nse, pars = as.numeric(param2[1,]) , lower = as.numeric(param2[2,]),upper = as.numeric(param2[3,]), maxn=500)))
    
    param.new = round(abcd.fit.nse$par,5)
    
    if(mod == 0){
      param.new<- c(round(abcd.fit.nse$par,5),-99,-99,-99)
    }else if (mod == 1){
      param.new = c(round(abcd.fit.nse$par,5),-1,-99)
    }else if (mod == 3){
      param.new = c(round(abcd.fit.nse$par[1:4],5),-99,abcd.fit.nse$par[5:6])
    }else{
      param.new = round(abcd.fit.nse$par,5)
    }
    subdata.mod = fun.abcd(mod = mod, subdata = subdata, param = param.new, supha = supha, subdata2 = subdata2)
  }
  
  #ahora evaluo la bondad de ajuste
  adjust = suppressWarnings(gof(sim=subdata.mod$Qhmsim,obs=subdata.mod$qhmobs,digits=6)[c(3,4,6,9,17,19),1])
  if(f0 %in% c("NSE","R2","KGE")){
    adjust[["f0"]] = - adjust[[f0]]
  }else{
    adjust[["f0"]] = adjust[[f0]]
  }
  basins[basins$code == subb,c("MSE","RMSE","PBIAS","NSE","R2","KGE","f0")] = adjust

  #y parametros:
  basins[basins$code == subb,c("param")] = paste(param.new,collapse = ";")
  
  #y aquí vamos generando el nuevo data frame  
  if(i == 1){
    data.mod = subdata.mod
  }else{
    data.mod = bind_rows(data.mod,subdata.mod)
  }
}


#exportar
if(save.data == 1){
  write_csv(data.mod,file="ABCDresultsSEMIDIST.txt")
}

write_csv(basins,file="ABCDadjustSEMIDIST.txt")

