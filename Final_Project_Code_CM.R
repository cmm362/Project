##########################
##The purpose of this function is to take bacterial DNA sequencing data from fecal DNA isolates
##and using input provided by the user, determine the changes in relative abundance of each taxon
##detected in the fecal sample over time. The output is a csv file with the percent change in abundance as well as p values
##from paired t tests between the samples at different time points found under the "Results" folder.
##There is also a graph output in the "Figures" folder showing the percent change for each
##taxon that changed to a significant (p < 0.05) degree in either treatment group.
##Some potential user input is provided in parentheses next to each prompt.
##########################

#Define function to detemine whether packages are installed and install/load them if needed
pkgTest <- function(x)
{
  if (!require(x, character.only = TRUE))
  {
    install.packages(x, dep=TRUE)
    if(!require(x, character.only = TRUE)) stop("Package not found")
  }
}

#Run function on all necessary packages
pkgTest("ggplot2")
pkgTest("reshape")
pkgTest("dplyr")

##Take simple user input to find the correct files for the desired experiment, treatment groups and time points
ExperimentID <- readline(prompt="Experiment (CM002.2): ")
Group1 <- readline(prompt="Group 1 (Control): ")
Group2 <- readline(prompt="Group 2 (Treatment): ")
Group1Base = readline(prompt="Baseline timepoint for Group 1 (0): ") #Set baseline timepoint for group 1
Group2Base = readline(prompt="Baseline timepoint for Group 2 (0): ") #Set baseline timepoint for group 2
Group1End = readline(prompt="End timepoint for Group 1 (3, 7, or 14): ") #Set experimental timepoint for group 1
Group2End = readline(prompt="End timepoint for Group 2 (3, 7, or 14): ") #Set experimental timepoint for group 2

#Find sequencing data files for each taxon level and put them in a list
SeqDataFileNames <- grep("otu_table_filtered_L[0-9]+", list.files(path = "./"), value = TRUE)

#Extract the specific taxon levels that exist for the data set (for use with labeling graphs and columns)
TaxonLevelsAvailable <- sub("otu_table_filtered_", "", SeqDataFileNames)
TaxonLevelsAvailable <- sub(".csv", "", TaxonLevelsAvailable)
SeqDataFileCount <- length(SeqDataFileNames)

#Make the function loop to produce a csv and graph for each taxon level
for (TaxonLevel in 1:SeqDataFileCount) {
  TargetSeqDataFile <- paste("./", SeqDataFileNames[TaxonLevel], sep = "", collapse = "")

#Rename columns in data file as CM Experiment numbers
seqdata <- read.csv(TargetSeqDataFile, header = TRUE)
TaxonCount <- nrow(seqdata)
ColumnCount <- ncol(seqdata)
KLnames <- read.csv("./KL_Map.csv", header = TRUE)
KLnamesfrom <- as.vector(KLnames$Sequencing_Tag)
KLnamesto <- as.vector(KLnames$Sequencing_ID)
for (n in 1:500) {
  colnames(seqdata)[colnames(seqdata) == KLnamesfrom[n]] <- KLnamesto[n]
  n <- n + 1
}

# Rename columns in data file as Experiment_Treatment_Sample_Day
CMnames <- read.csv("./CM_Map.csv", header = TRUE)
CMnamesfrom <- as.vector(CMnames$Sequencing_ID)
CMnamesexp <- as.vector(CMnames$Experiment)
CMnamestreatment <- as.vector(CMnames$Treatment)
CMnamessample <- as.vector(CMnames$Sample)
CMnamesday <- as.vector(CMnames$Day)
for (a in 1:500) {
  colnames(seqdata)[colnames(seqdata) == CMnamesfrom[a]] <- paste(CMnamesexp[a], CMnamestreatment[a], CMnamesday[a], CMnamessample[a], sep = "_")
  a = a + 1
}

#Define vectors containing the column names for all samples related to the experiment
TargetExperiment <- grep(paste(ExperimentID, "_", sep = "", collapse = ""), colnames(seqdata), value = TRUE)

#Define vectors containing the column names for each treatment group
TargetGroup1 <- grep(paste("_", Group1, "_", sep = "", collapse = ""), TargetExperiment, value = TRUE)
TargetGroup2 <- grep(paste("_", Group2, "_", sep = "", collapse = ""), TargetExperiment, value = TRUE)

#Make vectors for each Experimental group for given time points in a consistent sample order (sort is needed to run paired t test)
TargetGroup1Base <- sort(grep(paste("_", Group1Base, "_", sep = "", collapse = ""), TargetGroup1, value = TRUE))
TargetGroup2Base <- sort(grep(paste("_", Group2Base, "_", sep = "", collapse = ""), TargetGroup2, value = TRUE))
TargetGroup1End <- sort(grep(paste("_", Group1End, "_", sep = "", collapse = ""), TargetGroup1, value = TRUE))
TargetGroup2End <- sort(grep(paste("_", Group2End, "_", sep = "", collapse = ""), TargetGroup2, value = TRUE))

#Find average difference between time points for group 1
DifferenceGroup1Results <- 1:TaxonCount
for (Taxon in 1:TaxonCount) {
  DifferenceGroup1Results[Taxon] <- 100 * mean(((as.numeric(seqdata[Taxon, TargetGroup1End])) - mean(as.numeric(seqdata[Taxon, TargetGroup1Base])))/mean(as.numeric(seqdata[Taxon, TargetGroup1Base])))
  Taxon <- Taxon + 1
}

#Find average difference between time points for group 2
DifferenceGroup2Results <- 1:TaxonCount
for (Taxon in 1:TaxonCount) {
  DifferenceGroup2Results[Taxon] <- 100 * ((mean(as.numeric(seqdata[Taxon, TargetGroup2End])) - mean(as.numeric(seqdata[Taxon, TargetGroup2Base])))/mean(as.numeric(seqdata[Taxon, TargetGroup2Base])))
  Taxon <- Taxon + 1
}

#Perform Paired T tests for each taxon (row) in Group 1
TTestGroup1Results <- 1:TaxonCount
for (Taxon in 1:TaxonCount) {
  TTestGroup1Results[Taxon] <- t.test(as.numeric(seqdata[Taxon, TargetGroup1Base]), 
                                      as.numeric(seqdata[Taxon, TargetGroup1End]), 
                                      paired = TRUE)$p.value
  Taxon <- Taxon + 1
}

#Perform Paired T tests for each taxon (row) in Group 2
TTestGroup2Results <- 1:TaxonCount
for (Taxon in 1:TaxonCount) {
  TTestGroup2Results[Taxon] <- t.test(as.numeric(seqdata[Taxon, TargetGroup2Base]), 
                                      as.numeric(seqdata[Taxon, TargetGroup2End]), 
                                      paired = TRUE)$p.value
  Taxon <- Taxon + 1
}

#Convert differences into matrix
DifferenceOutput <- matrix(c(as.numeric(DifferenceGroup1Results), as.numeric(DifferenceGroup2Results)), nrow = TaxonCount)
rownames(DifferenceOutput) <- as.character(seqdata[1:TaxonCount, 1])
colnames(DifferenceOutput) <- c(paste(Group1, "_days_", Group1Base, "_to_", Group1End, sep = "", collapse = ""),
                                 paste(Group2, "_days_", Group2Base, "_to_", Group2End, sep = "", collapse = ""))

#Convert p values from T test into matrix
PValueOutput <- matrix(c(as.numeric(TTestGroup1Results), as.numeric(TTestGroup2Results)), nrow = TaxonCount)
rownames(PValueOutput) <- as.character(seqdata[1:TaxonCount, 1])
colnames(PValueOutput) <- c(paste(Group1, "_days_", Group1Base, "_to_", Group1End, sep = "", collapse = ""),
                                paste(Group2, "_days_", Group2Base, "_to_", Group2End, sep = "", collapse = ""))

#Use melt function on each matrix and then combine the matrices so that it is in a format that ggplot can interpret
MeltedMatrixDifferences <- melt(DifferenceOutput)
MeltedMatrixPValues <- melt(PValueOutput)

FinishedMatrix <- matrix(nrow = TaxonCount * 2, ncol = 4)
FinishedMatrix[,1] <- as.character(MeltedMatrixDifferences[,1])
FinishedMatrix[,2] <- as.character(MeltedMatrixDifferences[,2])
FinishedMatrix[,3] <- MeltedMatrixDifferences[,3]
FinishedMatrix[,4] <- MeltedMatrixPValues[,3]
colnames(FinishedMatrix) <- c("Taxon", "Group", "Difference", "P_value")

##Convert results matrix into CSV and save with title reflecting the comparison that was done
system("mkdir Results") ## will print a warning if folder exists, but won't stop code
CSVOutputFileName <- paste("Results/", TaxonLevelsAvailable[TaxonLevel], "_", ExperimentID, "_", Group1, "_Days_", Group1Base, "_to_", Group1End, "_vs_", Group2, "_Days_", Group1Base, "_to_", Group1End, ".csv", sep = "", collapse = "")
write.csv(FinishedMatrix, file = CSVOutputFileName)

#Bar plot of p values for each taxon, omitting all taxa that had NaN as the result (because the taxon was not detected in any of the samples for the comparison)
#First, melt the df to a format that can be interpreted by ggplot with three columns: the taxon, the group it belongs to, and the p value from the paired t test
DF <- na.omit(read.csv(file = CSVOutputFileName, header = TRUE))
NewTitles <- sub(';c__', paste(';', '
      ', 'c__'), DF[,2])
NewTitles <- sub(';f__', paste(';', '
      ', 'f__'), NewTitles)
DF[,2] <- NewTitles
DF <- filter(DF, DF[,5] < 0.05)
DF[,5]

#Then make the bar plot, with the x axis as taxon, the y axis as p value, and the color as group
#It looks a little weird to have super wide bars wherever only one group had an omitted data point, but I actually like it better that way because it differentiates between NA and a very low p value, which would have extremely different consequences
PlotTitle <- paste(ExperimentID, " ", TaxonLevelsAvailable[TaxonLevel], sep = "", collapse = "")
DifferencesPlot <- ggplot(data = DF, aes(x = Taxon, y = Difference, fill = Group)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        legend.key.size = unit(0.5, "cm"),
        legend.title=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "gray"),
        axis.line.y = element_line(colour = "gray")) + 
  xlab("Taxon") + 
  ylab("% Difference over Time") + 
  ggtitle(PlotTitle)

##Overwrite original CSV with a version that's not usable by ggplot but looks better for analysis
EndMatrix <- matrix(nrow = TaxonCount, ncol = 5)
EndMatrix[,1] <- DifferenceOutput[,1]
EndMatrix[,2] <- PValueOutput[,1]
EndMatrix[,3] <- DifferenceOutput[,2]
EndMatrix[,4] <- PValueOutput[,2]
rownames(EndMatrix) <- as.character(seqdata[1:TaxonCount, 1])
colnames(EndMatrix) <- c("Taxon", 
                         paste("Difference_", Group1, "_days_", Group1Base, "_to_", Group1End, sep = "", collapse = ""), 
                         paste("P_Value_", Group1, "_days_", Group1Base, "_to_", Group1End, sep = "", collapse = ""), 
                         paste("Difference_", Group1, "_days_", Group1Base, "_to_", Group1End, sep = "", collapse = ""), 
                         paste("P_Value_", Group1, "_days_", Group1Base, "_to_", Group1End, sep = "", collapse = ""))
write.csv(EndMatrix, file = CSVOutputFileName)

##Export this plot as jpeg
system("mkdir Figures") ## will print a warning if folder exists, but won't stop code
GGPlotOutputFileName <- paste("Figures/", TaxonLevelsAvailable[TaxonLevel], "_", ExperimentID, "_", Group1, "_Days_", Group1Base, "_to_", Group1End, "_vs_", Group2, "_Days_", Group1Base, "_to_", Group1End, ".jpg", sep = "", collapse = "")
ggsave(GGPlotOutputFileName, DifferencesPlot)

TaxonLevel <- TaxonLevel + 1
}