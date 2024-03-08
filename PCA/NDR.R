#data needs to be rows samples and columns methylation sites
#the files in the methylation_coverage folder contain data that 
#looks the following
#each row: chr start_position end_position percentage_methylation 
#methylated_cytosines #unmethylated_cytosines

setwd("/home/khiltunen/public/4Students/KaisuHiltunen/PCA")
# Load necessary libraries
library(tidyverse) #contains ggplot2 dplyr and tidyR for plotting
library(stats)
library(impute) #for data imputation
library(readxl) #for reading excel data
library(dplyr)
library(viridis) #for creating colors in the plots (does not work yet)
library(tools) #for extracting the names


#Replace InputFolder with actual data
inputFolder = "../euler_exports/S183/grcm38/"
sampleFiles = list.files(inputFolder,full.names = T,recursive = T)
sampleFiles = sampleFiles[grep(sampleFiles,pattern = ".cov$|.cov.gz$")]


#read the mouse clock sites and create a dataframe with the sites
msClock <- read.table("msClockSites_mm10.bed", header = FALSE, 
                      sep = "\t",stringsAsFactors=FALSE, quote="")
msClock <- data.frame(paste(msClock[,1],msClock[,2],sep=":"))
df_full<-data.frame(matrix(ncol=1, nrow=nrow(msClock)))

#rename the rows
rownames(df_full) <- msClock[,1]
#insert a column called chr for merging the data from sample files 
df_full[,1] <- msClock[,1]
#rename the column as chr
colnames(df_full) <- c("chr")

#Create a vector for the sample file column names
colNames<-c("chr", "start_position", "end_position", "percentage_methylation",
            "methylated_cytosines", "unmethylated_cytosines")



#loop through all teh files in the folder and merge the rows with the 
#matching chromosome location
for(f in sampleFiles){
  #read data from one file
  scInfo <- read.delim(f,as.is=T,header=F, col.names=colNames,
                       colClasses = c("character","double","double",
                                      "double","double","double"))
  #modify first column to contain chromosome and start position
  scInfo[,1] = data.frame(paste(scInfo[,1],scInfo[,2],sep=":"))
  #add the sample into the matrix but only keep the clock sites
  df_full<-merge(df_full,scInfo[,c("chr","percentage_methylation")],
                 by="chr",all.x=TRUE)
  #remove the file path
  sampleName <- tools::file_path_sans_ext(basename(f))
  #edit the name to look like the rownames in age info
  colnames(df_full)[colnames(df_full) == "percentage_methylation"] <- sampleName
}


### Function for filtering all rows that contain more than 20 % NaN values
# input matrix contains data from multiple files
filter<-function(df, coverage, n_unique){
  # Count the number of non-NaN values in each row
  non_na_counts <- rowSums(!is.na(df));
  # Identify rows with at least 80% coverage
  selected_rows <- which(non_na_counts >= coverage * ncol(df));
  #create anew df with only these rows
  filtered_df <- df[selected_rows, ];
  # Transpose the dataframe to have samples as rows and genomic regions as columns 
  tdf <- t(filtered_df);
  
  #remove low variance columns
  low_variance_cols <- which(apply(tdf, 2, function(col) var(col, na.rm = TRUE) < 10))
  selected_columns <- !colnames(tdf) %in% low_variance_cols
  # Remove near-zero variance columns
  tdf <- tdf[,selected_columns]
  # Remove columns with less than n_unique unique values
  unique_counts <- lengths(apply(X=tdf, MARGIN=2, FUN=unique))
  tdf <- tdf[,which(unique_counts>n_unique)]
  return(tdf)
}


###A function for extracting the PCA scores###
perform_pca<-function(my_tdf){
  ### Preparing data for PCA ###
  #save rownames but drop the first one which is the chr
  rowNames<-rownames(my_tdf)[-1];
  rowNames
  # extract shortnames
  pattern <- "^(\\d+_\\w+_\\w+)_GRCm.*"
  sampleNames <- sub(pattern, "\\1", rowNames);
  sampleNames
  #transform to matrix
  matrixM <- as.matrix(my_tdf[,]);
  #now we can add the positions into colnames
  colnames(matrixM)=matrixM[1,];
  #and remove the redundant chr row
  matrixM=matrixM[-1,];
  
  # Convert all columns to numeric for analysis
  matrixM<-apply(matrixM, 2, as.numeric);
  # Impute missing values
  imputed_data <- impute.knn(matrixM);
  imputedMatrix <- imputed_data[["data"]];
  # Standardize the numeric columns
  stMatrix <- scale(imputedMatrix);
  ### Perform PCA ###
  rownames(stMatrix)<-sampleNames;
  pca_result <- prcomp(stMatrix, scale. = TRUE)
  # Access the principal components
  pc_scores <- as.data.frame(pca_result$x[, 1:2])
  
  ### Gathering all meta information ###
  #add back the rownames for plotting purposes
  pc_scores$ShortName<-sampleNames
  
  #extract tissue txpe
  tissueType <- sub(".*_(\\w+)_r\\d+_.*", "\\1", rownames(pc_scores))
  pc_scores$TissueType <-tissueType
  #read metadata
  Meta_data_file <- read_excel("Meta_data_file.xlsx")
  
  # Find rows where ShortNames match the names of samples
  matching_rows <- grepl(paste(sampleNames,collapse="|"), Meta_data_file$ShortName)
  ages<-Meta_data_file[matching_rows,c("ShortName","Age")]
  #add age
  pc_scores<-merge(pc_scores,ages[,c("ShortName","Age")],by="ShortName",all=TRUE)
  
  # Extract the first two letters of each shortname to distinguish replicates
  pc_scores$AnimalId <- substr(pc_scores$ShortName, 1, 2)
  n_sites=ncol(matrixM)
  return(list(pc_scores,n_sites))
}

### Plotting function 1 ###
plot1 <- function(pc_scores_list){
  pc_scores<-pc_scores_list[[1]]
  n_sites<-pc_scores_list[2]
  ggplot() +
    geom_point(data=pc_scores, size = 3, aes(x = PC1, y = PC2, color = pc_scores$Age, shape = TissueType)) +
    #geom_text(pc_scores,aes(label = AnimalId))+
    labs(title = paste("PCA of methylation percentage with",n_sites, "sites"),x = "Principal Component 1", y = "Principal Component 2",
         color = "Age", shape = "Tissue Type") +
    scale_color_gradient(low = "#000000", high = "#A4C8F6") +  # Use a color gradient for age
    scale_shape_manual(values = c("BAT" = 19, "Blood" = 17, "Liver" = 15, "scAT" = 8)) +
    theme_minimal() + # Use a minimal theme for a clean look
    theme(legend.position = "right")+
    guides(size = FALSE)
}

### Plotting function 2 ###
#this plot is meant for individual tissues
plot2 <- function(pc_scores_list){
  pc_scores <- pc_scores_list[[1]]
  n_sites <- pc_scores_list[2]
  plot2 <- ggplot(pc_scores, aes(x = PC1, y = PC2, color= Age, group = AnimalId)) +#shape = AnimalId)) +
    geom_point() +
    
    geom_line(aes(linecolor=pc_scores$AnimalId)) +
    labs(title = paste("PCA of methylation percentage with", n_sites, "sites from", pc_scores$TissueType[1], "samples"),
         x = "Principal Component 1", y = "Principal Component 2",
         color = "Age") +
    scale_color_gradient(low = "#000000", high = "#A4C8F6") +
    #scale_shape_manual(values = shapes)  +
    #scale_color_manual(values = animal_colors)+
    theme_minimal() +
    theme(legend.position = "right")
  
  print(plot2)
  
}
### All data ###
tdf <- filter(df_full,0.9,3);

pc_scores_all <-perform_pca(tdf);
plot1(pc_scores_all)
#plot2(pc_scores_all)


### Liver data ###
## filter from df_full just the liver data
#extract tissue txpe
tissueType <- sub(".*_(\\w+)_r\\d+_.*", "\\1", colnames(df_full));
#keep also the first columns which has chr information
liver_columns <- c(1,which(tissueType == "Liver"));
tdf_liver <- filter(df_full[,liver_columns],0.9,4);

pc_scores_liver <- perform_pca(tdf_liver);
plot1(pc_scores_liver)
plot2(pc_scores_liver)



### BAT data ###
#extract tissue txpe
BAT_columns <- c(1,which(tissueType == "BAT"));
tdf_BAT <- filter(df_full[,BAT_columns],0.9,4);
#Do PCA and plot
pc_scores_BAT <- perform_pca(tdf_BAT);
plot1(pc_scores_BAT)
plot2(pc_scores_BAT)

### scAT data ###
scAT_columns <- c(1,which(tissueType == "scAT"));
tdf_scAT <- filter(df_full[,scAT_columns],0.9,4);
#Do PCA and plot
pc_scores_scAT <- perform_pca(tdf_scAT);
plot1(pc_scores_scAT)
plot2(pc_scores_scAT)

### Blood data ###
blood_columns <- c(1,which(tissueType == "Blood"));
tdf_blood <- filter(df_full[,blood_columns],0.9,4);
#Do PCA and plot
pc_scores_blood <- perform_pca(tdf_blood);
plot1(pc_scores_blood)
plot2(pc_scores_blood)


# Visualize the results graveyard of functions
#maybe with autoplot
plot(pc_scores$PC1, pc_scores$PC2, main = "PCA of Methylation Percentage Data", 
     xlab = "Principal Component 1", ylab = "Principal Component 2")
# Plot actual age vs. predicted age using ggplot
ggplot(pc_scores, aes(x = PC1, y = PC2, color = TissueType)) +
  geom_point() +
  labs(title = paste("PCA of methylation percentage with",ncol(matrixM), "sites")) +
  scale_color_manual(values = c("BAT" = "red", "Blood" = "blue", "Liver" = "green", "scAT" = "purple")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
#scatter plot for all the tissues
ggplot() +
  geom_point(data=pc_scores, aes(x = PC1, y = PC2, color = pc_scores$Age, shape = TissueType, size = 100)) +
  labs(title = paste("PCA of methylation percentage with",ncol(matrixM), "sites"),x = "Principal Component 1", y = "Principal Component 2",
       color = "Age", shape = "Tissue Type") +
  scale_color_gradient(low = "#000000", high = "#A4C8F6") +  # Use a color gradient for age
  scale_shape_manual(values = c("BAT" = 16, "Blood" = 17, "Liver" = 0, "scAT" = 8)) +
  theme_minimal() + # Use a minimal theme for a clean look
  theme(legend.position = "right")




