require(readr)
require(hydroGOF)

init_global_single = function(data_file, basin_file, basin_code){
    #datos
    data <- read_csv(data_file, show_col_types = FALSE)
    data$date <- as.Date(paste0(data$date,"01"), format = "%Y%m%d")
    data <- data[data[,"code"] == basin_code,]
    data_single <<- data[order(data$date),]
    basin_info_single <<- read_csv(basin_file, show_col_types = FALSE)
    supha_single <<- basin_info_single[basin_info_single[,"code"] == basin_code,"supha"]
    basin_code_single <<- basin_code
}

get_basin_q_single = function(mod, param) {
    #2. Ejecución del modelo#####################################################

    #y resto
    # data[,"Qhmsim"] <- NA
    
    
    wi <- 0
    yi <- 0
    eti <- 0
    rgi <- 0
    roi <- 0
    Sgi <- 0
    Swi <- 0
    Spti <- 0
    Qmmsim <- 0
    Qhmsim <- rep(0, nrow(data_single))

    pre = data_single[,"pre"]
    tas = data_single[,"tas"]
    etp = data_single[,"etp"]

    for(i in 1:nrow(data_single)){
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
        aux1 <- wi + param[2]
        aux2 <- 2*param[1]
        yi <- aux1/aux2 - ((aux1/aux2)^2 - ((wi * param[2]) / param[1]))^0.5
        
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
        Qmmsim <- roi + Qbi
        Qhmsim[i] <- Qmmsim*supha_single/100000
    }

    #3. Estadisticos de bondad de ajuste en calibración############
    
    adjust = suppressWarnings(gof(sim=as.numeric(Qhmsim), obs=data_single$qhmobs, digits=6)[c(3,4,6,9,17,19),1])
    
    adjust
}

eval_basin_param_single = function(mod, param) {
  Qhmsim <- get_basin_q_single(mod, param)

  subdata <- data_single[data_single[,1] == basin_code_single,]
  subdata_filtered <- subdata[subdata$qhmobs != -100,]
  Qhmsim_filtered <- Qhmsim[subdata$qhmobs != -100]

  suppressWarnings(gof(sim=Qhmsim_filtered,obs=subdata_filtered$qhmobs,digits=6)[c(3,4,6,9,17,19),1])
}