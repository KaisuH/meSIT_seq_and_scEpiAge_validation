library(Rsamtools)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tools) # For handling strings
library(ggridges)#used for ridge density plotting which was not used in the end
library(car)#variance analysis
setwd("C:/Users/khiltunen/Documents/R/coverage_analysis")




### All data ###
# Replace with the actual path to your folder
my_folder_path <- "~/R/coverage_analysis/S183_bam_coverage/"
S183<-create_data_frame(my_folder_path)
S184<-create_data_frame("~/R/coverage_analysis/S184_bam_coverage/")
S185<-create_data_frame("~/R/coverage_analysis/S185_bam_coverage/")
create_data_frame<-function(folder_path){
  # Get a list of file names in the folder
  file_names <- list.files(path = folder_path, pattern = ".txt", full.names = TRUE)
  
  #read the mouse clock sites and create a dataframe with the sites
  msClock <- read.table("~/R/PCA/msClockSites_mm10.bed", header = FALSE, 
                        sep = "\t",stringsAsFactors=FALSE, quote="")
  msClock <- msClock %>%
    mutate(chr = paste(V1, V2, sep = ":")) %>%
    dplyr::select(chr)
  
  # Create a data frame with a single column named "chr"
  df_full <- data.frame(chr = msClock$chr)
  
  #rename the rows
  rownames(df_full) <- msClock[,1]
  colNames <-  c("chr", "start_position", "end_position", "read", "coverage")
  df_full
  for(f in file_names){
    #read data from one file
    sample_data <- read.delim(f,as.is=T,header=F, col.names=colNames,
                         colClasses = c("character","double","double",
                                        "double","double"))
    sample_data <- sample_data %>%
      # Modify first column to contain chromosome and start position
      mutate(chr = paste(chr, start_position, sep = ":")) %>%
      dplyr::select(chr, coverage) %>%
      # Since these files contain coverage of both read 1 and read 2 
      # combine them into one number presenting the total coverage
      group_by(chr)%>%
      summarize(coverage=sum(coverage))
    # Merge with df_full, keeping only clock sites
    df_full <- right_join(df_full, sample_data, by = "chr")
    
    #remove the file path
    sampleName <- tools::file_path_sans_ext(basename(f))
    # extract shortnames
    pattern <- "^(\\d+_\\w+_\\w+)_GRCm.*"
    shortName <- sub(pattern, "\\1", sampleName);
    # Rename columns
    colnames(df_full)[colnames(df_full) == "coverage"] <- shortName
  }
  rownames(df_full)<-df_full[,"chr"]
  df_full <- df_full[,-1]
  return(df_full)
}
###heatmap stuff###
heatmap_matrix <- df_full %>%
  select_if(is.numeric) %>%
  as.matrix() 

heatmap(heatmap_matrix, scale="column", Colv = NA, Dendogram = FALSE)

pheatmap(heatmap_matrix, treeheight_col=0, treeheight_row=0)
ggplot(head(coverage_df,n=100), aes(x = Start, y = Coverage))+
  geom_point()+
  labs(x = "Genomic Position", y = "Coverage")+
  ggtitle("Sequencing Data Coverage Plot")

sampleID <- unique(sub("^(.*?_.*?)(_.*)$", "\\1", colnames(df)))



#plot the histogram of the coverage values in each libary
# mean average covereage of the libraries
means <- c(mean(as.matrix(S183)), mean(as.matrix(S184)),mean(as.matrix(S185)))

rmean_S183<-rowMeans(S183)
rmean_S184<-rowMeans(S184)
rmean_S185<-rowMeans(S185)
# mean coverages per row of all three libraries

df <- data.frame(
  Site = rep(row.names(S183), 3),  # Assuming your sites are numeric or can be represented as such
  Mean_coverage = c(rmean_S183, rmean_S184, rmean_S185),
  Avg_density = means,
  Enrichment = rep(c("Custom", "xGen", "none"), each = length(rmean_S183))
)

df_1 <- data.frame(
  Site = rep(row.names(S183), 1),  # Assuming your sites are numeric or can be represented as such
  Mean_coverage = c(rmean_S185),
  Avg_density = means[3],
  Enrichment = rep(c("none"), each = length(rmean_S183))
)
df_2 <- data.frame(
  Site = rep(row.names(S183), 2),  # Assuming your sites are numeric or can be represented as such
  Mean_coverage = c(rmean_S183, rmean_S184),
  Avg_density = means[1:2],
  Enrichment = rep(c("Custom", "Xgen"), each = length(rmean_S183))
)

# Normalize coverage values by library-specific read depth
#Normalized Coverage= Raw Coverage/Library-Specific Read Depth



df_2_normalized_coverage <- data.frame(
  Site = rep(row.names(S183), 2),  # Assuming your sites are numeric or can be represented as such
  Mean_coverage = c(rmean_S183/C_library_depth, rmean_S184/X_library_depth),
  Avg_density = c(mean(as.matrix(rmean_S183)/C_library_depth), 
                  mean(as.matrix(rmean_S184)/X_library_depth)),
  Enrichment = rep(c("Custom", "Xgen"), each = length(rmean_S183))
)


plot_histogram(df_2_normalized_coverage)
plot_boxplot(df_2_normalized_coverage)
plot_histogram <- function(df){
  ggplot(data = df, aes(x = Mean_coverage, fill = Enrichment)) +
    stat_bin(bins = 40, position = "identity", alpha = 0.5, aes(y = ..count../sum(..count..))) +
    scale_fill_manual(values = c("#2978a0", "#a13d63")) +
    labs(x = "Coverage", y = "Normalized Frequency", title = "Histogram of Clock Site Coverage Values", fill = "Enrichment method") +
    theme_minimal()+
    theme(
      legend.position = c(.98, .98),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6),
      plot.title = element_text(size = 15),  # Adjust title size
      axis.title.x = element_text(size = 14),  # Adjust x-axis label size
      axis.title.y = element_text(size = 14),  # Adjust y-axis label size
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),# Adjust axis text size
      legend.title = element_text(size = 13),  # Adjust legend title size
      legend.text = element_text(size = 12)  # Adjust legend text size
    )
}
#by groups
df$Enrichment <- factor(df$Enrichment, levels = c("none", "xGen", "Custom"))
#COLORS
#80ded9, #315659, #2978a0, #dd8f46, #a13d63
# RIGDE DENSITY
ggplot(df, aes(x = Mean_coverage, y = Enrichment)) +
  geom_density_ridges(aes(fill = Enrichment),linewidth=0.5, scale=4, alpha=0.6) +
  scale_fill_manual(values = c("#80ded9", "#a13d63","#2978a0")) +
  geom_vline(data = df, aes(xintercept = Avg_density), linetype = "dashed", size = 0.5) +  # Add average density line
  #geom_segment(aes(x = Avg_density, y = 0, xend = Avg_density, yend=Enrichment),linetype = "dashed",linewidth=0.5)+
  #geom_vline(xintercept = df$Avg_density, yintercept = df$Enrichment, linetype = "dashed", size = 0.5) +  # Add horizontal dashed line
  labs(group = "Enricment method", x = "Coverage", y = "Density", title = "Density Plot of Clock Site Coverage Values") +
  theme(legend.position = "none")+
  xlim(0, 37)+  # Set x-axis limits
  theme_minimal()+
  theme(
    plot.title = element_text(size = 15),  # Adjust title size
    axis.title.x = element_text(size = 14),  # Adjust x-axis label size
    axis.title.y = element_text(size = 14),  # Adjust y-axis label size
    axis.text.y = element_blank(),  # Adjust axis text size
    legend.title = element_text(size = 13),  # Adjust legend title size
    legend.text = element_text(size = 12)  # Adjust legend text size
  )
# HISTOGRAM
ggplot(data = df_2, aes(x = Mean_coverage, fill = Enrichment)) +
  geom_histogram(binwidth = 3, position = "identity", alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("#2978a0", "#a13d63")) +
  geom_vline(aes(xintercept = Avg_density), linetype = "dashed", size = 0.5, color = "red") +  # Add average density line
  geom_vline( aes(xintercept = means[1]), linetype = "dashed", color="#2978a0", size = 0.5) +  # Add average density line
  geom_vline( aes(xintercept = means[2]), linetype = "dashed", color="#a13d63", size = 0.5) +  # Add average density line
  labs(x = "Coverage", y = "Frequency", title = "Histogram of Clock Site Coverage Values", fill = "Enrichment method") +
  theme_minimal()+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    plot.title = element_text(size = 15),  # Adjust title size
    axis.title.x = element_text(size = 14),  # Adjust x-axis label size
    axis.title.y = element_text(size = 14),  # Adjust y-axis label size
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),# Adjust axis text size
    legend.title = element_text(size = 13),  # Adjust legend title size
    legend.text = element_text(size = 12)  
    )


# BOXPLOT
plot_boxplot <- function(df){
  ggplot(data = df, aes(x = Enrichment, y = Mean_coverage, fill = Enrichment)) +
    geom_boxplot(alpha = 0.5, width = 0.5) +
    scale_fill_manual(values = c("#80ded9", "#a13d63","#2978a0")) +
    geom_hline(yintercept = Avg_density, linetype = "dashed", size = 0.5, color = "red") +  # Add average density line
    labs(x = "Enrichment Method", y = "Coverage", title = "Boxplot of Clock Site Coverage Values") +
    theme_minimal()+
    guides(fill=FALSE)+
    theme(
      plot.title = element_text(size = 15),  # Adjust title size
      axis.title.x = element_text(size = 14),  # Adjust x-axis label size
      axis.title.y = element_text(size = 14),  # Adjust y-axis label size
      axis.text.y = element_text(size = 9),  # Adjust axis text size
      axis.text.x = element_text(size = 12),
      legend.title = element_text(size = 13),  # Adjust legend title size
      legend.text = element_text(size = 12)  # Adjust legend text size
    )
}
  
 
# NORMALISED HISTOGRAM
ggplot(data = df, aes(x = Mean_coverage, fill = Enrichment)) +
  stat_bin(binwidth = 3, position = "identity", alpha = 0.5, aes(y = ..count../sum(..count..)), color = "black") +
  scale_fill_manual(values = c("#2978a0", "#a13d63")) +
  #geom_vline(aes(xintercept = Avg_density), linetype = "dashed", size = 0.5, color = "red") +  # Add average density line
  geom_vline( aes(xintercept = means[1]), linetype = "dashed", color="#2978a0", size = 0.5) +  # Add average density line
  geom_vline( aes(xintercept = means[2]), linetype = "dashed", color="#a13d63", size = 0.5) +  # Add average density line
  labs(x = "Coverage", y = "Normalized Frequency", title = "Histogram of Clock Site Coverage Values", fill = "Enrichment method") +
  theme_minimal()+
  theme(
    legend.position = c(.98, .98),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    plot.title = element_text(size = 15),  # Adjust title size
    axis.title.x = element_text(size = 14),  # Adjust x-axis label size
    axis.title.y = element_text(size = 14),  # Adjust y-axis label size
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),# Adjust axis text size
    legend.title = element_text(size = 13),  # Adjust legend title size
    legend.text = element_text(size = 12)  # Adjust legend text size
  )

### Statistical Analysis ###
# Perform a t-test
t_test_result <- t.test(rmean_S183, rmean_S184)

# Print the result
print(t_test_result)

# Perform a t-test for normalised data
t_test_result <- t.test(rmean_S183/C_library_depth, rmean_S184/X_library_depth)

# Print the result
print(t_test_result)
### ANOVA ### 
# Create a dataframe
data <- data.frame(
    Values = c(rmean_S183, rmean_S184),
    Group = rep(c("S183", "S184"), each = length(rmean_S183))
)
# Perform ANOVA of libraries
aov_result <- aov(Values ~ Group, data = data)

# Display the results
print(aov_result)
summary(aov_result)

### Variance analysis of sample data with Levene test
C_df <- libs$C_df
C_samples <- 
var(rmean_S183)
var(rmean_S184)
sd(rmean_S183)
sd(rmean_S184)
levene_test_result <- leveneTest(Values ~ Group, data = data)
print(levene_test_result)


### Demultiplexing analysis ###
# Define the file location for demultiplexing statistics
demux_file <- "~/R/coverage_analysis/demultiplexing_statistics.txt"
# Read the file to a table
demux_df <- read.table(demux_file, header = TRUE, sep = "\t")

get_logsums <- function(df){
  # A function to calculate the respective number of reads 
  # for Lambda, NC, and samples in the data frame df.
  # The sum is converted to logarithmis scale for plotting.
  Lam <- log(sum(df[grep("^La_", df$sample_id), ]$templates))
  NC <- log(sum(df[grep("^NC_", df$sample_id), ]$templates))
  Samples <- log(sum(df[grep("^[01]", df$sample_id), ]$templates))
  Unmatched <- log(sum(df[grep("unmatched", df$sample_id), ]$templates))
  # Create a data frame
  result <- data.frame(
    category = factor(c("Lambda", "Samples", "NC", "Unmatched"), levels = c("Lambda", "Samples", "NC", "Unmatched")),
    value = c(Lam, Samples, NC, Unmatched)
  )
  return(result)
}

get_sums <- function(df){
  # A function to calculate the respective number of reads 
  # for Lambda, NC, and samples in the data frame df.
  Lam <- sum(df[grep("^La_", df$sample_id), ]$templates)
  NC <- sum(df[grep("^NC_", df$sample_id), ]$templates)
  Samples <- sum(df[grep("^[01]", df$sample_id), ]$templates)
  Unmatched <- sum(df[grep("unmatched", df$sample_id), ]$templates)
  # Create a data frame
  result <- data.frame(
    category = factor(c("Lambda", "Samples", "NC", "Unmatched"), levels = c("Lambda", "Samples", "NC", "Unmatched")),
    value = c(Lam, Samples, NC, Unmatched)
  )
  return(result)
}

get_library <- function(df){
  # A function to subset the data frame to include 
  # only samples ending with "_C", "_X", or "_N"
  
  # Extract samples ending with "_C" and the subsequent "unmatched" rows
  selected <- df[grep("_C$|unmatched", df$sample_id), ]
  # Subset the data frame up to the last occurrence
  first_index <- min(grep("_C$", selected$sample_id))
  last_index <- max(grep("_C$", selected$sample_id))+1
  C_df <- selected[1:last_index, ]
  
  # Repeat for X and N
  selected <- df[grep("_X$|unmatched", df$sample_id), ]
  first_index <- min(grep("_X$", selected$sample_id))
  last_index <- max(grep("_X$", selected$sample_id))+1
  X_df <- selected[first_index:last_index, ]
  
  selected <- df[grep("_N$|unmatched", df$sample_id), ]
  first_index <- min(grep("_N$", selected$sample_id))
  last_index <- max(grep("_N$", selected$sample_id))+1
  N_df <- selected[first_index:last_index, ]
  
  reslist <- list(C_df, X_df, N_df)
  names(reslist) <- c("C_df", "X_df", "N_df")
  return(reslist)
}
# Construct the tables for demultiplexing statistics of custom (C), xgen (X),
# and no enrichment (N) libraries
libs <- get_library(demux_df)

# Calculate logarithm of sum of the total number of templates
N_plot_data <- get_logsums(libs$N_df)
C_plot_data <- get_logsums(libs$C_df)
X_plot_data <- get_logsums(libs$X_df)

# Calculate logarithm of sum of the total number of templates
N_plot_data <- get_sums(libs$N_df)
C_plot_data <- get_sums(libs$C_df)
X_plot_data <- get_sums(libs$X_df)

# Compare template distributions of Lambda, Neg Control, and 
# Savmples. Within libraries. Plot a bar plot.
bar_plot(X_plot_data, "xGen")

bar_plot(C_plot_data, "Custom")

bar_plot(N_plot_data, "non-enriched")

#COLORS
#80ded9, #315659, #2978a0, #dd8f46, #a13d63

bar_plot <- function(data, library){
  # Function to create a bar plot for one library
  ggplot(data, aes(x = category, y = value, fill = category)) +
    geom_bar(stat = "identity", color = "black") +
    theme_minimal() +
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb","#a13d63"))+
    labs(title = paste0("Distribution of Reads in ",library, " library"), 
         fill = "Category", x = "Category", y = "Logarithm of Number of Reads")
}


compiled_bar_plot <- function(data){
  # Function to create a bar plot for all libraries
  ggplot(data, aes(x = category, y = value, fill = category)) +
    geom_bar(stat = "identity", color = "black") +
    theme_minimal() +
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) +  # Customize colors
    labs(title = "Distribution of Reads", fill = "Category", x = "Category", y = "Logarithm of Number of Reads")
}

sum(C_df$templates)
C_average_reads <- mean(C_df$templates)
X_average_reads <- mean(X_df$templates)
N_average_reads <- mean(N_df$templates)
N_library_depth <- sum(N_df$templates)
C_library_depth <- sum(C_df$templates)
X_library_depth <- sum(X_df$templates)

# These are the numbers that will normalise the coverage 
N_library_depth
C_library_depth
X_library_depth
# Print the result
print(average_reads)

#Analysing the coverage of the selected area of 11:12756587-120641815
# Read the coverage output file
coverage_data <- read.table("nonzero_coverage_chr11.bed", header = FALSE, col.names = c("chromosome", "start", "end", "position", "coverage"))
coverage_data <- coverage_data %>%
  mutate(position = position + 12756586) #make the position match the actual position in the chr
chr11_sites <- subset(msClock, V1 == 11)[1:2]
names(chr11_sites)[2] <- "position"
chr11_sites <- subset(chr11_sites, position >= 12756586)
chr11_sites$y <-rep(0,length(chr11_sites$V1))
chr11_sites$is_near_clock_site <- rep(TRUE,length(chr11_sites$V1))
# Add metadata indicating if coverage base is among the clock sites
# Create a logical column indicating whether the position is among clock sites
coverage_data <- coverage_data %>%
  mutate(is_clock_site = position %in% chr11_sites$position)
# create a logical column indicating if the position is near clock_site
coverage_data <- coverage_data %>%
  mutate(is_near_clock_site = sapply(position, function(pos) any(abs(chr11_sites$position - pos) <= 200)))
         

# Order data so that the true values are displayed on top
coverage_data <- coverage_data %>%
  arrange(is_clock_site)
library(scales)
# Plot the coverage
ggplot() +
  geom_col(data=coverage_data, aes(x = position, y = coverage, fill = factor(is_near_clock_site, levels = c(FALSE,TRUE))), width = 4) +
  labs(title = "Coverage plot of a section in Chr 11", x = "Position", y = "Coverage")+
  scale_fill_manual(values = c("#fc8d62", "#66c2a5")) +  # Customize colors
  # Add clock sites to the plot
  geom_point(data=chr11_sites,aes(x=chr11_sites$position, y=chr11_sites$y,shape = chr11_sites$is_near_clock_site, color=chr11_sites$is_near_clock_site))+
  scale_shape_manual(values=c(20))+
  scale_color_manual(values = c("#66c2a5")) +  # Customize colors
  theme_minimal()+
  guides(shape=FALSE)+
  labs(fill = "Within 200 bp of a clock site")+
  labs(color = "Clock site")+
  xlim(82000000,86000000)+
  scale_x_continuous(labels = scales::number_format(scale = 1e-6))

         