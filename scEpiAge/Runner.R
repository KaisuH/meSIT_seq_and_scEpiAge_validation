###### Settings ###### 
setwd("~/public/4Students/KaisuHiltunen/scEpiAge")
options(bitmapType='cairo')
library(ggplot2)
##Tissue, options are: liver, lung, blood.
tissue="blood"

##Folder with COV files.
folderInput = "../euler_exports/S183/grcm38/"
#folderInput = "Gravina_2016/SingleCell/"

## output file with aging information.
#outputFileName="./preditionGravinaBulk.txt"
outputFileName= paste("./predS183_",tissue,".txt",sep="")
print(outputFileName)
##If there are known age calculate confidence on the difference observed.
nSimulations = 5

#this does not have to be changed. I modified the code so that it retrieves the correct rows
ageInfo = read.delim("./allSampleInfo.txt")
#ageInfo = read.delim("./GravinaBulkInfo.txt")
#ageInfo = read.delim("./GravinaScInfo.txt")
#ageInfo = read.delim("./")


## ouput file name with extended information.
outputFileNameExtended <- paste("./extended_", sub("^\\./", "", outputFileName), sep = "")
# outputFileNameExtended="./preditionGravinaBulk.extended.txt"
#outputFileNameExtended="./preditionGravinaSc.extended.txt"

plot=T

###################### 


## Actual code form here:

calcExtendedStats = F
if(!is.na(outputFileNameExtended) & !is.null(ageInfo)){
  calcExtendedStats=T
} else {
  print("Not given age information or output filename for the extended analysis")
}

###### Source ########
##source prediction functions.
source("./PredictionFunctions/functions.R")
###################### 

##### Load data ######
if(tissue=="liver"){
  expectedMethMatrix <- read.delim("./ExpectedMethylationMatrices/ExpectedMethMat_Liver.tsv",as.is=T,row.names=1,check.names = F)
  backupInformation <- read.delim("./SiteInformation/clockSites_Liver_BackUp.txt",as.is=T)
} else if(tissue=="lung"){
  expectedMethMatrix <- read.delim("./ExpectedMethylationMatrices/ExpectedMethMat_Lung.tsv",as.is=T,row.names=1,check.names = F)
  backupInformation <- read.delim("./SiteInformation/clockSites_Lung_BackUp.txt",as.is=T)
} else if(tissue=="blood"){
  expectedMethMatrix <- read.delim("./ExpectedMethylationMatrices/ExpectedMethMat_Blood.tsv",as.is=T,row.names=1,check.names = F)
  backupInformation <- read.delim("./SiteInformation/clockSites_Blood_BackUp.txt",as.is=T)
} else{
  print("No valid tissue selected")
  stop();
}

##Predict.
inputMethMatrix = readCovFiles(folderInput,backupInformation);

if(calcExtendedStats){
  predictionOutVersusExpected = predictAgesAndCalculateExpectedGivenAge(inputMethMatrix, backupInformation, expectedMethMatrix, ageInfo, nSimulations,plot);
  write.table(predictionOutVersusExpected,outputFileNameExtended,sep="\t",quote=F)
} else {
  predictionOut = predictAges(inputMethMatrix, backupInformation, expectedMethMatrix);
  write.table(predictionOut,outputFileName,sep="\t",quote=F)  
}

data <- read.delim(outputFileNameExtended)
#tissue = "blood"
#data <- read.delim("extended_predS183_blood.txt")

#tissue = "lung"
#data <- read.delim("extended_predS183_lung.txt")
#print(outputFileNameExtended)
# Extract the columns for actualAge and predictedAge
actualAge <- data$actualAge
predictedAge <- data$predictedAge

#print(rownames(data))
# Assuming that the sample names are in a column named "sampleName"
# You may need to adjust this based on your actual column name
tissueType <- sub(".*_(\\w+)_r\\d+_.*", "\\1", rownames(data))
data$TissueType <-tissueType


# Plot actual age vs. predicted age using ggplot
ggplot(data, aes(x = actualAge, y = predictedAge, color = TissueType)) +
  geom_point() +
  labs(title = paste("Actual vs Predicted Age with",tissue, "clock")) +
  scale_color_manual(values = c("BAT" = "red", "Blood" = "blue", "Liver" = "green", "scAT" = "purple"))
  #+theme_light()


ageDeviation <- data$ageDeviation
# Plot ageDeviation vs. actual age using ggplot

ggplot(data, aes(x = actualAge, y = ageDeviation, color = TissueType)) +
  geom_point() +
  labs(title = paste("Predicition deviation vs. actual age",tissue, "clock")) +
  scale_color_manual(values = c("BAT" = "red", "Blood" = "blue", "Liver" = "green", "scAT" = "purple"))

# View the data 
View(dropSites)
print(data)