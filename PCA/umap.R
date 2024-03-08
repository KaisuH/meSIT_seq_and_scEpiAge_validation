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
library(tibble) #for replacing rownames

metaDataFile<-"~/R/PCA/Meta_data_file.xlsx"
#Replace InputFolder with actual data
inputFolder = "~/R/euler_exports/S183/grcm38"
sampleFiles = list.files(inputFolder,full.names = T,recursive = T)
sampleFiles = sampleFiles[grep(sampleFiles,pattern = ".cov$|.cov.gz$")]


#read the mouse clock sites and create a dataframe with the sites
msClock <- read.table("~/R/PCA/msClockSites_mm10.bed", header = FALSE, 
                      sep = "\t",stringsAsFactors=FALSE, quote="")
msClock <- msClock %>%
  mutate(chr = paste(V1, V2, sep = ":")) %>%
  select(chr)

# Create a data frame with a single column named "chr"
df_full <- data.frame(chr = msClock$chr)

#rename the rows
rownames(df_full) <- msClock[,1]

#Create a vector for the sample file column names
colNames<-c("chr", "start_position", "end_position", "percentage_methylation",
            "methylated_cytosines", "unmethylated_cytosines")

#loop through all teh files in the folder and merge the rows with the 
#matching chromosome location
for(f in sampleFiles){
  #read data from one file
  sample_data <- read.delim(f,as.is=T,header=F, col.names=colNames,
                       colClasses = c("character","double","double",
                                      "double","double","double"))
  # Modify first column to contain chromosome and start position
  sample_data <- sample_data %>%
    mutate(chr = paste(chr, start_position, sep = ":")) %>%
    select(chr, percentage_methylation)
  # Merge with df_full, keeping only clock sites
  df_full <- left_join(df_full, sample_data, by = "chr")
  #remove the file path
  sampleName <- tools::file_path_sans_ext(basename(f))
  #edit the name to look like the rownames in age info
  df_full <- df_full %>%
    rename_with(~ sampleName, matches("percentage_methylation"))
}



### Function for filtering all rows that contain more than 20 % NaN values
# input matrix contains data from multiple files
filter_data<-function(df, min_coverage, n_unique){
  # Count the number of non-NaN values in each row.
  # Identify rows with at least 80% min_coverage
  # & create anew df with only these rows.
  # Transpose the dataframe to have samples as rows and genomic regions as columns 
  tdf <- df %>%
    filter(rowSums(!is.na(.[])) >= min_coverage * ncol(.)) %>%
    t() %>%
    as.data.frame()
  #remove low variance columns
  tdf <- tdf %>%
    as.data.frame() %>%
    select_if(~ var(.[-1], na.rm = TRUE) >=10) %>%
    select_if(~ length(unique(.)) > n_unique) %>%
    as.data.frame()
  return(tdf)
}

###A function for extracting the PCA scores###
impute_and_standardise <- function(my_tdf){
  ### Preparing data for PCA ###
  #save rownames but drop the first one which is the chr
  rowNames<-rownames(my_tdf)[-1];
   # extract shortnames
  pattern <- "^(\\d+_\\w+_\\w+)_GRCm.*"
  shortNames <- sub(pattern, "\\1", rowNames);
  
  #prepare the data matrix
  matrixM <- my_tdf %>%
    slice(-1) %>% #drop first row with the chr info
    rownames_to_column(var = "rowNames") %>% #Convert row names to a column named "rowNames"
    select(-rowNames) %>% #Exclude the "rowNames" column
    mutate_all(as.numeric) %>% # Convert all columns to numeric for analysis
    as.matrix() %>%# Transform to matrix
    impute.knn() # Impute missing values
  
  
  # Scale the imputed data
  stMatrix <- matrixM$data %>%
    scale() # Standardize the numeric columns
  rownames(stMatrix)=shortNames
  
  return(stMatrix)
}

perform_pca <- function(stMatrix){
  # This function takes as input the standardised data matrix and 
  # performs the pca analysis
  shortNames <- rownames(stMatrix)
  ### Perform PCA ###
  pca_result <- prcomp(stMatrix, scale. = TRUE)
  # Access the principal components
  pc_scores <- as.data.frame(pca_result$x[, 1:2])
  
  #Gather metadata
  Meta_data_file <- read_excel(metaDataFile)
  ages<-Meta_data_file[,c("ShortName","Age")]
  #extract tissue types into a vector
  tissueType <- sub(".*_(\\w+)_r\\d+_.*", "\\1", rownames(pc_scores))
  
  ### Gathering all meta information ###
  # Add back the row names for plotting purposes
  pc_scores <- pc_scores %>%
    mutate(ShortName = shortNames)%>% #add back the rownames for plotting purposes
    mutate(TissueType = tissueType) %>% #add the tissuetypes
    left_join(ages %>% select(ShortName, Age), by = "ShortName") %>% # Add age and animal Id
    mutate(AnimalId = substr(ShortName, 1, 2))
  
  n_sites=nrow(stMatrix)
  return(list(pc_scores,n_sites))
  
}

perform_umap <- function(stMatrix){
  # This function takes as input the standardised data matrix and 
  # performs the umap analysis
  shortNames <- rownames(stMatrix)
  ### Perform PCA ###
  umap_res <- umap(stMatrix)
  # Access the principal components
  #pc_scores <- as.data.frame(pca_result$x[, 1:2])
  umap_data <- as.data.frame(umap_res$layout)
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  #Gather metadata
  Meta_data_file <- read_excel(metaDataFile)
  ages<-Meta_data_file[,c("ShortName","Age")]
  #extract tissue types into a vector
  tissueType <- sub(".*_(\\w+)_r\\d+_.*", "\\1", rownames(umap_data))
  
  ### Gathering all meta information ###
  # Add back the row names for plotting purposes
  umap_data <- umap_data %>%
    mutate(ShortName = shortNames)%>% #add back the rownames for plotting purposes
    mutate(TissueType = tissueType) %>% #add the tissuetypes
    left_join(ages %>% select(ShortName, Age), by = "ShortName") %>% # Add age and animal Id
    mutate(AnimalId = substr(ShortName, 1, 2))
  
  n_sites=ncol(stMatrix)
  umap_data
  
  umap_data
  return(list(umap_data,n_sites))
  
}
### Plotting function 1 ###
plot1 <- function(pc_scores_list){
  pc_scores<-pc_scores_list[[1]]
  n_sites<-pc_scores_list[2]
  ggplot() +
    geom_point(data=pc_scores, size = 3, aes(x = PC1, y = PC2, color = Age, shape = TissueType)) +
    #geom_text(pc_scores,aes(label = AnimalId))+
    labs(title = paste("PCA of methylation percentage with",n_sites, "sites"),x = "Principal Component 1", y = "Principal Component 2",
         color = "Age", shape = "Tissue Type") +
    scale_color_gradient(low = "#000000", high = "#A4C8F6") +  # Use a color gradient for age
    scale_shape_manual(values = c("BAT" = 19, "Blood" = 17, "Liver" = 15, "scAT" = 8)) +
    theme_minimal() + # Use a minimal theme for a clean look
    theme(legend.position = "right")+
    guides(size = none)
}

### Plotting function 2 ###
#this plot is meant for individual tissues
plot2 <- function(pc_scores_list){
  pc_scores <- pc_scores_list[[1]]
  n_sites <- pc_scores_list[2]
  plot2 <- ggplot(pc_scores, aes(x = PC1, y = PC2, color= Age, group = AnimalId)) +#shape = AnimalId)) +
    geom_point() +
    
    geom_line(aes(linecolour=AnimalId)) +
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

plot_umap <- function(umap_data){
  umap_scores<-umap_data[[1]]
  n_sites<-umap_data[2]
  #ggplot(umap_data, aes(x = V1, y = V2))+
  #  geom_point()+
  #  theme_minimal()
  ggplot() +
    geom_point(data=umap_scores, size = 3, aes(x = UMAP1, y = UMAP2, color = Age, shape = TissueType)) +
    #geom_text(umap_scores,aes(label = AnimalId))+
    labs(title = paste("UMAP of methylation percentage with",n_sites, "sites"),x = "UMAP 1", y = "UMAP 2",
         color = "Age", shape = "Tissue Type") +
    scale_color_gradient(low = "#000000", high = "#A4C8F6") +  # Use a color gradient for age
    scale_shape_manual(values = c("BAT" = 19, "Blood" = 17, "Liver" = 15, "scAT" = 8)) +
    theme_minimal() + # Use a minimal theme for a clean look
    theme(legend.position = "right")+
    guides(size = none)
}
plot_umap_individual <- function(umap_data){
  umap_scores <- umap_data[[1]]
  n_sites <- umap_data[2]
  plot2 <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2, color= Age, group = AnimalId)) +#shape = AnimalId)) +
    geom_point() +
    
    geom_line(aes(linecolour=AnimalId)) +
    labs(title = paste("UMAP of methylation percentage with", n_sites, "sites from", umap_scores$TissueType[1], "samples"),
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
pc_scores_liver <- df_full %>%
  filter_data(0.9,3)%>%
  impute_and_standardise()%>%
  perform_pca();
plot1(pc_scores_all)
#plot2(pc_scores_all)
##umap
umap_all <- perform_umap(stM)
head(umap_all)
plot_umap(umap_all)
head(umap_all)

###Individual tissue analysis###
analyse_tissue <- function(tissue, df){
  tissue_columns <- c(1,which(tissueType == tissue));

  pc_scores_tissue <- df[,tissue_columns] %>%
    filter_data(0.9,5)%>%
    impute_and_standardise()%>%
    perform_pca();
  plot1(pc_scores_tissue)
  plot2(pc_scores_tissue)
  umap_liver <- df[,tissue_columns] %>%
    filter_data(0.9,5)%>%
    impute_and_standardise()%>%
    perform_umap();
  plot_umap_individual(umap_tissue)
}
### Liver data ###
#analyse_tissue("Liver",df_full)

## filter from df_full just the liver data
#extract tissue txpe
tissueType <- sub(".*_(\\w+)_r\\d+_.*", "\\1", colnames(df_full));
tissueType

#keep also the first columns which has chr information
liver_columns <- c(1,which(tissueType == "Liver"));
liver_columns
pc_scores_liver <- df_full[,liver_columns] %>%
  filter_data(0.9,5)%>%
  impute_and_standardise()%>%
  perform_pca();
plot1(pc_scores_liver)
plot2(pc_scores_liver)
umap_liver <- df_full[,liver_columns] %>%
  filter_data(0.9,5)%>%
  impute_and_standardise()%>%
  perform_umap();
plot_umap_individual(umap_liver)

### BAT data ###
#extract tissue txpe
BAT_columns <- c(1,which(tissueType == "BAT"));
tdf_BAT <- filter_data(df_full[,BAT_columns],0.9,4);
#Do PCA and plot
pc_scores_BAT <- perform_pca(tdf_BAT);
plot1(pc_scores_BAT)
plot2(pc_scores_BAT)

### scAT data ###
scAT_columns <- c(1,which(tissueType == "scAT"));
tdf_scAT <- filter_data(df_full[,scAT_columns],0.9,4);
#Do PCA and plot
pc_scores_scAT <- perform_pca(tdf_scAT);
plot1(pc_scores_scAT)
plot2(pc_scores_scAT)

### Blood data ###
blood_columns <- c(1,which(tissueType == "Blood"));
tdf_blood <- filter_data(df_full[,blood_columns],0.9,4);
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


