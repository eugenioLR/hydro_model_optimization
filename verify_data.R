require(readr)
require(hydroGOF)

#datos
data <- read_csv("data3045.txt", show_col_types = FALSE)
data$date = as.Date(paste0(data$date,"01"),format = "%Y%m%d")
data = data[order(data$date),]
adjust = suppressWarnings(gof(sim=data$Qmmsim, obs=data$qmmobs, digits=6)[c(3,4,6,9,17,19),1])