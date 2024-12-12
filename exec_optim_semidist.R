require(readr)
require(hydroGOF)

# #datos
# data <- read_csv("data/CHGdataSIMPA.txt", show_col_types = FALSE)
# # data <- read_csv("data/dataSimpa.txt", show_col_types = FALSE)
# data$date = as.Date(paste0(data$date,"01"), format = "%Y%m%d")
# # data <- data[data[,1] == 5029 | data[,1] == 5054 | data[,1] == 5060 | data[,1] == 5071 | data[,1] == 5043,]
# data = data[order(data$date),]

# basins = as.data.frame(read_csv("data/CHGbasins.txt"))
# # basins = as.data.frame(read_csv("data/basinsSimpa.txt"))


# code,order,codedown,supha
# 5029,1,5043,16168.94
# 5054,1,5043,8257.84
# 5060,1,5043,15684.26
# 5071,1,5043,12170.17
# 5043,2,0,277148.04

data <- NA

init_global = function(data_file, basin_file){
    #datos
    data <- read_csv(data_file, show_col_types = FALSE)
    data$date <- as.Date(paste0(data$date,"01"), format = "%Y%m%d")
    data <<- data[order(data$date),]
    basins <<- as.data.frame(read_csv(basin_file))
}

get_basin_q = function(mod, param, basin_code, prev_q) {
    #2. Ejecución del modelo#####################################################

    subdata <- data[data[,1] == basin_code,]
    curr_basin <- basins[basins[,1] == basin_code,]
    supha <- as.numeric(curr_basin[4])

    #y resto
    #data[,"Qmmsim"]<-NA
    Qmmsim <- rep(0, nrow(subdata))
    
    wi <- 0
    yi <- 0
    eti <- 0
    rgi <- 0
    roi <- 0
    Sgi <- 0
    Swi <- 0
    Spti <- 0

    pre = subdata[,"pre"]
    tas = subdata[,"tas"]
    etp = subdata[,"etp"]

    for(i in 1:nrow(subdata)){
        #A) RUTINA WASMOD (modos 2 o 3)
        if(mod >= 2){

            #nieve
            snt <- max(pre[i,1] * (1 - exp((tas[i,1] - param[6])/(param[6] - param[7]))^2), 0)
            
            #lluvia
            rt <- pre[i,1]-snt

            #deshielo
            mt <- max(pre[i,1] * (1 - exp((tas[i,1] - param[7])/(param[6] - param[7]))^2), 0)

            #almacenamiento
            Spti <- Spti + snt - mt

            #Y voy a agregar a la precipitacion el deshielo
            liquido <- rt + mt

            wi <- liquido + Swi
        }else{
            wi <- Swi + pre[i,1]
        }
        
        #Componente Yi
        aux1 <- (wi + param[2]) / (2*param[1])
        yi <- aux1 - ((aux1)^2 - ((wi * param[2]) / param[1]))^0.5
        
        #Componente Swi
        Swi <- yi * exp(-etp[i,1]/param[2])
        eti <- yi - Swi

        #descomp salida
        aux <- wi - yi

        rgi <- param[3] * aux
        roi <- (1-param[3]) * aux


        #C) SEGUNDO EMBALSE
        Sgi <- (rgi+Sgi)/(1+param[4])
        Qgi <- param[4] * Sgi

        #mod si utilizamos Xg o Xg+Wasmod (modos 1 o 2) o el resto:
        if(mod == 0 || mod == 3){
            Qbi <- Qgi
        }else{
            Qbi <- param[5] * Qgi
        }
        
        #suma descarga y esc
        Qmmsim[i] <- roi + Qbi
    }
    #3. Estadisticos de bondad de ajuste en calibración############
    Qhmsim <- as.numeric(Qmmsim)*supha/100000 + as.numeric(prev_q)
    Qhmsim
}

eval_basin_param = function(mod, param, basin_code, prev_q) {
  Qhmsim <- get_basin_q(mod, param, basin_code, prev_q)
  subdata <- data[data[,1] == basin_code,]

  (subdata_filtered <- subdata[subdata$qhmobs != -100,])
  (Qhmsim_filtered <- Qhmsim[subdata$qhmobs != -100])

  suppressWarnings(gof(sim=Qhmsim_filtered,obs=subdata_filtered$qhmobs,digits=6)[c(3,4,6,9,17,19),1])
}

eval_basin = function(Qhmsim, basin_code) {
  subdata <- data[data[,1] == basin_code,]

  subdata_filtered <- subdata[subdata$qhmobs != -100, ]
  Qhmsim_filtered <- Qhmsim[subdata$qhmobs != -100]

  suppressWarnings(gof(sim=as.numeric(Qhmsim_filtered),obs=subdata_filtered$qhmobs,digits=6)[c(3,4,6,9,17,19),1])
}
