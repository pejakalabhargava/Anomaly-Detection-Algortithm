install.packages("stringr")
install.packages("sna")
install.packages("igraph")

source("F:\\Fall 2014\\GDM\\P4\\P4_bkakran\\findAnamoly.R")
#feature = 'degree'
#feature = 'clustCoeff' 
#feature = 'egoNet' 
findAnamoly(fileInputPath="F:\\Fall 2014\\GDM\\P4\\data\\enron",
            fileOutPutPath="F:\\Fall 2014\\GDM\\P4\\output",
            featureToUse="egoNet",
            windowSizeLength=7)
