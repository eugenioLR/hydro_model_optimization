require(readr)
require(hydroGOF)

#datos
data <- read_csv("data/CHGdataSIMPA5043AG.txt", show_col_types = FALSE)
data$date = as.Date(paste0(data$date,"01"), format = "%Y%m%d")
data <- data[data[,1] == 5043,]
data <- data[order(data$date),]

basin_info <- read_csv("data/CHGbasins5043AG.txt", show_col_types = FALSE)
supha <- basin_info[,"supha"]

hydro_prob = function(mod, param) {
    #2. Ejecución del modelo#####################################################

    #y resto
    # data[,"Qhmsim"] <- NA
    
    
    Swi <- 0
    Sgi <- 0
    Spti <- 0
    wi <- 0
    yi <- 0
    eti <- 0
    rgi <- 0
    roi <- 0
    Sgi <- 0
    Swi <- 0
    Qmmsim <- 0
    Qhmsim <- rep(0, nrow(data))

    pre = data[,"pre"]
    tas = data[,"tas"]
    etp = data[,"etp"]

    for(i in 1:nrow(data)){
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
        # data[i,"Qhmsim"] <- roi + Qbi
        Qmmsim <- roi + Qbi
        Qhmsim[i] <- Qmmsim*supha/100000
    }

    #3. Estadisticos de bondad de ajuste en calibración############
    
    # adjust = suppressWarnings(gof(sim=data$Qmmsim, obs=data$qmmobs, digits=6)[c(3,4,6,9,17,19),1])
    adjust = suppressWarnings(gof(sim=as.numeric(Qhmsim), obs=data$qhmobs, digits=6)[c(3,4,6,9,17,19),1])
    
    adjust
}