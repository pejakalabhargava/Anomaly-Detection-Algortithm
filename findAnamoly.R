
#This method is used to calculate the feature required for the anamoly detection algorithm
#This takes a day's data as the input and the total number of vertices in the entire time
#series graph and then the required feature to be calculated for the day's data
#Usage Example: 
#getFeatureDataForADay(dataList = data,numberOfVertices = 10,feature = c("degree","clustCoeff","egoNet"))
getFeatureDataForADay <- function(dataList,numberOfVertices,feature) {
  #Create a matrix to hold feature data since there will be entry for each vertex as a column
  featureMatrix = matrix(0,ncol=numberOfVertices,nrow = 1); 
  #Represents the number of vertices in the system
  v=seq(0,numberOfVertices-1)
  #Create a subgraph induced by the day's data.However there will 
  #be entry for all vertices of the time series graph
  #Note:Graph is undirected
  graph =  graph.data.frame(dataList, directed=FALSE, vertices=v)
  if(feature == 'degree') {
    #if feature is degree get the degree using degree function of igraph package
    featureMatrix = degree(graph)
    featureMatrix = as.matrix(featureMatrix)
    #Create the transopse of matrix so that we have a column and rows equal to number of vertices
    featureMatrix = t(featureMatrix)
  } else if(feature == 'clustCoeff'){
    # Clustering coefficient is a measure of the degree to which nodes in a graph tend to cluster together.
    # It is is the proportion of given node's neighbors that can be reached by other neighbors
    # gives the clustering coefficient of each node
    #Transitivity measures the probability that the adjacent vertices of a vertex are connected.
    #This is sometimes also called the clustering coefficient.The local transitivity of a vertex is
    #the ratio of the triangles connected to the vertex and the triples centered on the vertex.
    #source:R Documentation
    clustCoefficient = transitivity(graph,type="local")
    clustCoefficient = as.matrix(clustCoefficient)
    featureMatrix = t(clustCoefficient)
    #replace Na in matrix with 0
    featureMatrix[is.na(featureMatrix)] = 0
    } else if(feature == 'egoNet') {
      # The egonet of a node is the subgraph induced by the node and its neighbors
      #We need to find number of edges in the egonet of the node
      #get the adjacency data from the graph
      adj=get.adjacency(graph)
      #Convert the adjacency list to adjacent matrix
      adj=data.matrix(adj, rownames.force = NA)
      #Calculate the egonet for each of the vertex in the graph.
      #Program uses sna package for this
      egoNet<-ego.extract(adj,neighborhood="combined")
      #Get the egonet for each vertex
      for (k in 1:numberOfVertices){
        #get the subscript as index as string
        indexForego =toString(k-1)
        #Get the egonet matrix
        mat = data.matrix(egoNet[[indexForego]])
        #Number of edges is degree/2
        featureMatrix[1,k]=length(which(mat==1))/2
      }
      #detach("package:sna", unload=TRUE)
    }
  #return the results
  return(as.matrix(featureMatrix))
}

#This method is used to calculcate the correlation matrix for the given input data
getCorrelationMatrix <- function(featureDataMatrix) {
  #Get rows and columns
  rows = nrow(featureDataMatrix)
  column = ncol(featureDataMatrix)
  featureDataMatrix=as.data.frame(featureDataMatrix)
  #Getcovariance Matrix
  covarianceMatrix=cov(featureDataMatrix)
  #Get column mean
  colMeansFeatureData = colMeans(featureDataMatrix)
  #Get column standard Deviation
  set.seed(42)
  colStdFeatureData =apply(featureDataMatrix, 2, sd)
  correlationmatrix = matrix(1,nrow=column,ncol=column)
  featureDataMatrix=as.matrix(featureDataMatrix)
  #calculate the correlation matrix using the standard forumla
  #make sure the standard deviation is not zero
  for (i in 1:column){
    for (j in 1:column){
      if(colStdFeatureData[i] != 0 && colStdFeatureData[j] != 0) {
        correlationmatrix[i,j] = covarianceMatrix[i,j]/(colStdFeatureData[i] * colStdFeatureData[j])
      } 
    }
  }
  return(as.matrix(correlationmatrix))
}
#This method uses R inbuilt function cor to calculate the pearson correaltion matrix
getPearsonCor <- function(featureDataMatrix) {
  #calculate the correlation matrix
  correlationMatrixLocal=cor(featureDataMatrix, use = "pairwise.complete.obs",method = "pearson")
  #Repalce na of matrix with 0
  #This is because because over the time window W of 7 days, most of the
  #nodes have no activity -their W-length vectors are all 0's
  #and thus the pair-wise correlations of such 0 vectors are computed to be 1.
  correlationMatrixLocal[is.na(correlationMatrixLocal)] = 1
  return(as.matrix(correlationMatrixLocal))
}

#This method is used to calculate the moving range average
#Moving Range Average is : |X1-X2| +|X2-X3| +------+|Xk-1 - Xk| / k-1  
getmovingRangeAverage <- function(zvector) {
  zAvg=0;
  sum =0
  #calculate sum |X1-X2| +|X2-X3| +------+|Xk-1 - Xk|  
  for (k in 2:length(zvector)){
    sum = sum +abs((zvector[k-1])-(zvector[k]))
  }
  #Find average
  zAvg = sum/(length(zvector)-1)
  return(zAvg)
}

findAnamoly <- function(fileInputPath,fileOutPutPath,featureToUse,windowSizeLength) {
library(stringr)
#WindowSize=7
WindowSize=windowSizeLength
#feature = 'degree'
#feature = 'clustCoeff' 
#feature = 'egoNet' 
feature=featureToUse
if(feature == 'egoNet') {
  library(sna)
  library(igraph)
# detach("package:igraph", unload=TRUE)
} else {
  library(igraph)
  library(sna)
  detach("package:sna", unload=TRUE)
}
  
#path = "F:\\Fall 2014\\GDM\\P4\\data\\reality_mining_voices"
#path = "F:\\Fall 2014\\GDM\\P4\\data\\enron"
path=fileInputPath
#get all the filenames present in the folder
filenames = list.files(path)
#Get all the numbers from the filenames
tmp <- sapply(filenames, function (k) strsplit(k, "[^0-9]"))
#Get only numbers from the tmp stripping other characters
tmp <- Reduce(union, tmp)
#convert them to numbers
tmp=as.numeric(tmp)
#sort numbers in increasing order
tmp=sort(tmp)
#Get the filename extension from the filenames eg:_enron_by_day.txt
filenameExt = "[_a-zA-Z.]{1,}"
filenameExt=str_extract(filenames[1], filenameExt)
numberOfVertices = 0
numberOfEdges = 0
data <- list() # creates a list
#create a list to hold the featuredata for each day
featureDataList <- list()
for (k in 1:length(filenames)){
  #Read data from the files in the order 0__enron_by_day.txt,1__enron_by_day.txt..etc
  data <- read.delim(sprintf("%s\\%s%s",path, tmp[k], filenameExt),header = FALSE,  sep = " ")
  #Get the number of edges for this day
  numberOfEdges = numberOfEdges + data[1,2]
  #Get the number of vertices in the entire graph
  numberOfVertices = data[1,1]
  #strip off the first line
  data = data[-1,]
  #cacluate the feature for the day and store it
  featureDataList[[k]]=getFeatureDataForADay(dataList = data,numberOfVertices = numberOfVertices,feature = feature) 
}
#get number of days
numberOfDays = length(tmp);
#get number of correlation matrices given by numberOfDays-windowSize+1
numberOfCorrelationMatrices = (numberOfDays-WindowSize)+1
#create matrix to hold data for each week consistibg of 7 rows and numberOfVertices as column
featureWeekData = matrix(0,ncol=numberOfVertices,nrow = 7)
currentDay = 1;
#Get data for first seven days
for (k in 1:WindowSize){
  featureWeekData[k,] = featureDataList[[k]]
  currentDay = currentDay +1
}

#caculcate the pearsons correlation matrix for the featureWeekData
correlationmatrix=getPearsonCor(featureDataMatrix=featureWeekData)
#We need matrix to hold principal eigen vector for each correlation matrix
principalEigenVector = matrix(0,nrow=numberOfVertices,ncol=numberOfCorrelationMatrices)
#Calculate eigen vectors for correlation data
EigenVectors<- eigen(correlationmatrix)
#Get the principal eigen vector
EigenVectors <- EigenVectors$vectors
#Store it for the first column
principalEigenVector[,1] = as.matrix(EigenVectors[,1])
currentCorelationMat = 2;
#caculcate the feature set for each day
for (k in 8:numberOfDays){
  #strip the first day from the previous week
  featureWeekData = featureWeekData[-1,]
  #calculate the feature data for current day
  featureDataTemp= featureDataList[[k]]
  #append the current day feature to the last six days
  featureWeekData = rbind(featureWeekData,featureDataTemp)
  #calculate the correlation matrix
  correlationmatrix = getPearsonCor(featureDataMatrix=featureWeekData)
  #Get principal eigen vector
  EigenVectors<- eigen(correlationmatrix)
  EigenVectors <- EigenVectors$vectors
  #save the correlation matrix's eigen vector
  principalEigenVector[,currentCorelationMat] = as.matrix(EigenVectors[,1])
  currentCorelationMat = currentCorelationMat +1;
}
#matrix to hold z vector
z = matrix(0,nrow=1,ncol=numberOfCorrelationMatrices)
#we use window W' as 5.hence start from the sixth principal eigen vector
for (k in 6:numberOfCorrelationMatrices){
#get last 5 eigen vectors
i=k-5
j=k-1
#calcuate the average of last five principal eigen vectors
averageEigenVector= rowMeans(principalEigenVector[,seq(i:j)])
averageEigenVector = as.matrix(averageEigenVector)
#get the currentEigen vector
currentEigenMatrix = principalEigenVector[,k]
currentEigenMatrix = as.matrix(currentEigenMatrix)
#Find the dot product of average eigen vector and current vector and caluclate x vector
z[,k] = 1-crossprod(averageEigenVector,currentEigenMatrix)
}
#Remove fist 5 values since it is zero and invalid
z = z[-c(1,2,3,4,5)]
#GEt moving range average of z scores
MRA = getmovingRangeAverage(z)
#get median of z scores
zMedian = median(z)
#Get the upper threshold
threshold = zMedian + (3 * MRA)
#Get z scor index which are more than threshold
anamolyIndex = which(z>=threshold)
#Get z score data which are more than threshold
anamolyData = z[which(z>=threshold)]
#Add 5 to index since days are forwarded by 5 to calculate z index
anamolyIndex = anamolyIndex + 5-1 
#combine results
anamolyData = cbind(anamolyIndex,anamolyData)
#sort the data in descending order based on z score data
anamolyData=anamolyData[order(anamolyData[,2],decreasing = TRUE),]
totalNumberOfAnamolies= length(anamolyIndex)

#Plot the z score against threshold
plotMain = sprintf("Anamoly Detection with Feature as %s",feature)
plot(z,main=plotMain,ylab = "Z Score",xlab = "Staring day for week")
abline(h=threshold,)

#Store the results in the output folder filename is generated as eg:output_egoNet_1415098449
#filename format output_<feature>_timestamp.txt
#outputFolder = "F:\\Fall 2014\\GDM\\P4\\output"
outputFolder = fileOutPutPath;
timeStamp = as.integer(Sys.time())
fileName = sprintf("output_day_%s_%s.txt",feature,timeStamp)
fullFilePath = sprintf("%s\\%s",outputFolder, fileName)

#Open file conenction
fileConn<-file(fullFilePath)
if(totalNumberOfAnamolies <= 10) {
  dataToOutput =c(totalNumberOfAnamolies,anamolyData[,1]);
}else if(totalNumberOfAnamolies<100) {
  dataToOutput =c(totalNumberOfAnamolies,anamolyData[seq(1:10),1]);
} else {
  totalNumberToOutPut = floor(totalNumberOfAnamolies/10)
  dataToOutput =c(totalNumberOfAnamolies,anamolyData[seq(1:totalNumberToOutPut),1]);
}
#Store the data
dataToOutput = as.character(dataToOutput)
writeLines(dataToOutput, fileConn)
close(fileConn)
}