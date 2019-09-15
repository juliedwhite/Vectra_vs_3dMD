#Analysis comparing 3dMDface to Vectra H1 using participant sparse landmark data
#Author: Julie D. White

##### Necessary packages and functions #####
library(tidyverse)
library(reshape2)
library(RRPP)
library(geomorph)
library(proxy)
library(plotly)
library(mixOmics)
library(ggpubr)
library(Morpho)
library(abind)
library(R.matlab)
library(nlme)

#Plot1Face script
#This script was written by Tomas Gonzalez-Zarzar and can be found at https://github.com/Arslan-Zaidi/Facial_masculinity_MHC/tree/master/Tutorial_visualizing_facial_heatmaps
source('Plot1Face.R')
RefScan <- read.table("FacialTemplate_20171201.obj", sep = "\t", header = F, colClasses = c("character", rep("numeric", 3)))
RefScan_Facets <- as.matrix(RefScan[which(RefScan$V1 == "f"), 2:4])
RefScan_Vertices <- as.matrix(RefScan[which(RefScan$V1 == "v"), 2:4])

#Remove RefScan since we don't need it anymore
rm(RefScan)

##### Read in the 19 landmark data #####
md.35ind.sparse <- readMat("RawImageData/3DMD/XYZ_35ind_3DMD_19LM.mat")
md.35ind.sparse <- md.35ind.sparse[[1]]

#Dimensions: 
#x = 35 (individuals)
#y = 57 (19 * 3 landmarks)
#z = 3 (replicates)
dimnames(md.35ind.sparse)[[1]] <- paste0("P", 1:35)
dimnames(md.35ind.sparse)[[2]] <- paste0(rep(c("x", "y", "z")), rep(1:19, each=3))
dimnames(md.35ind.sparse)[[3]] <- c("R1", "R2", "R3")

#Convert to geomorph format:
#x = number of landmarks (19)
#y = axes (3)
#z = individuals (105)

#First convert to 2d matrix that we can use array specs
#First replicate
r1 = md.35ind.sparse[,,1]
rownames(r1) <- paste0(rownames(r1), ".R1")

#Second replicate
r2 = md.35ind.sparse[,,2]
rownames(r2) <- paste0(rownames(r2), ".R2")

#Third replicate
r3 = md.35ind.sparse[,,3]
rownames(r3) <- paste0(rownames(r3), ".R3")

#bind back together
md.35ind.sparse2d <- rbind(r1, r2, r3)

#Convert to geomorph
md.35ind.sparse.geomorph <- arrayspecs(md.35ind.sparse2d, p=19, k=3, sep="")

#Remove clutter
rm(md.35ind.sparse, md.35ind.sparse2d, r1, r2, r3)

#Vectra
v.35ind.sparse <- readMat("RawImageData/Vectra/XYZ_35ind_Vectra_19LM.mat")
v.35ind.sparse <- v.35ind.sparse[[1]]

#Dimensions: 
#x = 35 (individuals)
#y = 57 (19 * 3 landmarks)
#z = 3 (replicates)
dimnames(v.35ind.sparse)[[1]] <- paste0("P", 1:35)
dimnames(v.35ind.sparse)[[2]] <- paste0(rep(c("x", "y", "z")), rep(1:19, each=3))
dimnames(v.35ind.sparse)[[3]] <- c("R1", "R2", "R3")

#Convert to geomorph format:
#x = number of landmarks (19)
#y = axes (3)
#z = individuals (105)

#First convert to 2d matrix that we can use array specs
#First replicate
r1 = v.35ind.sparse[,,1]
rownames(r1) <- paste0(rownames(r1), ".R1")

#Second replicate
r2 = v.35ind.sparse[,,2]
rownames(r2) <- paste0(rownames(r2), ".R2")

#Third replicate
r3 = v.35ind.sparse[,,3]
rownames(r3) <- paste0(rownames(r3), ".R3")

#bind back together
v.35ind.sparse2d <- rbind(r1, r2, r3)

#Convert to geomorph
v.35ind.sparse.geomorph <- arrayspecs(v.35ind.sparse2d, p=19, k=3, sep="")

#Remove clutter
rm(v.35ind.sparse, v.35ind.sparse2d, r1, r2, r3)

##### Align and average three replicate images #####
#3dMD
#Make a list containing the indices for each indiviudal. Should be a list of 35 with three items each. 
ind <- paste0("P", 1:35, "\\.")
ind.index <- lapply(ind, function(x) grep(x, dimnames(md.35ind.sparse.geomorph)[[3]]))

#3d array to store the GPA aligned coordinates
md.35ind.sparse.aligned <- NULL

#3d array to store the mshape, which is the centroid of the three replicates
md.35ind.sparse.aligned.centroid <- NULL

#List to store the names so we dont get lost
namelist <- NULL

#For each item in ind.index, pull out the three replicates for each image and convert to a 3d array. 
for (i in 1:length(ind.index)){
  #Pull out the three matrices that correspond to a single person
  tmp <- md.35ind.sparse.geomorph[,,ind.index[[i]]]
  #GPA align without scaling or reflecting
  tmp.gpa <- ProcGPA(tmp, scale = FALSE, reflection = FALSE)
  #Pull out the aligned coordinates and place them to a 3d array
  md.35ind.sparse.aligned <- abind(md.35ind.sparse.aligned, tmp.gpa$rotated, along = 3)
  #Pull out the dimnames and add them to a list of names
  namelist <- append(namelist, dimnames(tmp.gpa$rotated)[[3]])
  #Pull out the mean shape of the three replicate images and put into 3d array 
  md.35ind.sparse.aligned.centroid <- abind(md.35ind.sparse.aligned.centroid, tmp.gpa$mshape, along = 3)
}

#Assign names to the centroid of the replicates, using a pared down version of the namelist. The names in this list are in the same order as the data because they were constructed together. 
dimnames(md.35ind.sparse.aligned.centroid)[[3]] <- unique(unlist(strsplit(namelist, split = ".R[0-3]")))

#Remove the clutter
rm(ind, ind.index, namelist, tmp, tmp.gpa, i)

#Vectra
#Make a list containing the indices for each indiviudal. Should be a list of 35 with three items each. 
ind <- paste0("P", 1:35, "\\.")
ind.index <- lapply(ind, function(x) grep(x, dimnames(v.35ind.sparse.geomorph)[[3]]))

#3d array to store the GPA aligned coordinates
v.35ind.sparse.aligned <- NULL

#3d array to store the mshape, which is the centroid of the three replicates
v.35ind.sparse.aligned.centroid <- NULL

#List to store the names so we dont get lost
namelist <- NULL

#For each item in ind.index, pull out the three replicates for each image and convert to a 3d array. 
for (i in 1:length(ind.index)){
  #Pull out the three matrices that correspond to a single person
  tmp <- v.35ind.sparse.geomorph[,,ind.index[[i]]]
  #GPA align without scaling or reflecting
  tmp.gpa <- ProcGPA(tmp, scale = FALSE, reflection = FALSE)
  #Pull out the aligned coordinates and place them to a 3d array
  v.35ind.sparse.aligned <- abind(v.35ind.sparse.aligned, tmp.gpa$rotated, along = 3)
  #Pull out the dimnames and add them to a list of names
  namelist <- append(namelist, dimnames(tmp.gpa$rotated)[[3]])
  #Pull out the mean shape of the three replicate images and put into 3d array 
  v.35ind.sparse.aligned.centroid <- abind(v.35ind.sparse.aligned.centroid, tmp.gpa$mshape, along = 3)
}

#Assign names to the centroid of the replicates, using a pared down version of the namelist. The names in this list are in the same order as the data because they were constructed together. 
dimnames(v.35ind.sparse.aligned.centroid)[[3]] <- unique(unlist(strsplit(namelist, split = ".R[0-3]")))

#Remove the clutter
rm(ind, ind.index, namelist, tmp, tmp.gpa, i)

##### Camera error #####
#Align the two images per person
#Combine the replicate centroid arrays into a single array
b.35ind.sparse.replicate.centroid <- abind(md.35ind.sparse.aligned.centroid, v.35ind.sparse.aligned.centroid, along = 3)

#Alter the dimnames so that we know which camera each individual came from 
dimnames(b.35ind.sparse.replicate.centroid)[[3]] <- paste(dimnames(b.35ind.sparse.replicate.centroid)[[3]], rep(c("MD", "V"), each = 35), sep = ".")

#Make a list containing the indices for each indiviudal. Should be a list of 35 with two items each. 
ind <- paste0("P", 1:35, "\\.")
ind.index <- lapply(ind, function(x) grep(x, dimnames(b.35ind.sparse.replicate.centroid)[[3]]))

#3d array to store the GPA aligned coordinates
b.35ind.sparse.replicate.centroid.aligned <- NULL

#For each item in ind.index, pull out the two images ang align
for (i in 1:length(ind.index)){
  #Pull out the two matrices that correspond to a single person
  tmp <- b.35ind.sparse.replicate.centroid[,,ind.index[[i]]]
  #GPA align without scaling or reflecting
  tmp.gpa <- ProcGPA(tmp, scale = FALSE, reflection = FALSE)
  #Pull out the aligned coordinates and place them to a 3d array
  b.35ind.sparse.replicate.centroid.aligned <- abind(b.35ind.sparse.replicate.centroid.aligned, tmp.gpa$rotated, along = 3)
}

#Remove the clutter
rm(ind, ind.index, tmp, tmp.gpa, i)

#Calculate the euclidean distance between the aligned coordinates from each camera
#Make a list containing the indices for each indiviudal. Should be a list of 35 with two items each. 
ind <- paste0("P", 1:35, "\\.")
ind.index <- lapply(ind, function(x) grep(x, dimnames(b.35ind.sparse.replicate.centroid.aligned)[[3]]))

#Empty matrix to store the euclidean distances
b.35ind.sparse.camerror.euclid <- matrix(nrow = 35, ncol = 19)

#For each of the 35 people, calculate the euclidean distance between their 3dMD and Vectra coordinates and place into matrix
for (i in 1:length(ind.index)){
  b.35ind.sparse.camerror.euclid[i,] <- diag(dist(b.35ind.sparse.replicate.centroid.aligned[,,ind.index[[i]][1]], b.35ind.sparse.replicate.centroid.aligned[,,ind.index[[i]][2]], method = "euclidean"))
}

#Convert result matrix to dataframe
b.35ind.sparse.camerror.euclid <- as.data.frame(b.35ind.sparse.camerror.euclid)

#Create a Individual column for grouping
b.35ind.sparse.camerror.euclid$Individual <- unique(unlist(strsplit(dimnames(b.35ind.sparse.replicate.centroid.aligned)[[3]], split = c(".MD", ".V"))))

#Remove clutter
rm(ind, ind.index, i)

#Plot per person
#Reorder the factor levels by number
b.35ind.sparse.camerror.euclid$Individual <- factor(b.35ind.sparse.camerror.euclid$Individual, levels = c(paste0("P", 1:35)))

#Plot
p <- ggplot(melt(b.35ind.sparse.camerror.euclid, id.vars = "Individual"), aes(x = Individual, y = value))+geom_boxplot(fill = "#70B77E", outlier.size = 0.25)+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")+labs(x="Individual", y="Euclidean distance (mm)")
ggsave(filename = "FiguresAndTables/b.35ind.sparse.camerror.euclid.pdf", plot = p, device = "pdf", width = 3.25, height = 4.5, units = "in")

#Plot per landmark
lmk <- data.frame(c(1:19), c("Glabella", "Nasion", "Pronasale", "Subnasale", "Labiale superius", "Labiale inferius", "Pogonion", "Endocanthion left", "Endocanthion right", "Exocanthion left", "Exocanthion right", "Alar curvature left", "Alar curvature right", "Subalare left", "Subalare right", "Crista philtri left", "crista philtri right", "Cheilion left", "Cheilion right"))
colnames(lmk) <- c("LM.Num", "LM.Name")
b.35ind.sparse.camerror.euclid.perlmk <- as.data.frame(b.35ind.sparse.camerror.euclid[-20])
colnames(b.35ind.sparse.camerror.euclid.perlmk) <- lmk$LM.Name

p <- ggplot(melt(b.35ind.sparse.camerror.euclid.perlmk), aes(x = variable, y = value))+geom_boxplot(color = "#70B77E", outlier.colour = "black", outlier.size = 0.25)+theme_bw(base_size = 8)+theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))+labs(x="Landmark", y="Euclidean distance (mm)")
ggsave(filename = "FiguresAndTables/b.35ind.sparse.camerror.euclid.perlmk.pdf", plot = p, device = "pdf", width = 3.25, height = 4.5, units = "in")

#Remove clutter
rm(lmk, p)

##### ANOVA #####
#Put the 3dMD and Vectra images together in a single array
b.35ind.sparse.array <- abind(md.35ind.sparse.geomorph, v.35ind.sparse.geomorph, along = 3)

#Name the items
dimnames(b.35ind.sparse.array)[[3]] <- c(paste0("MD.", dimnames(md.35ind.sparse.geomorph)[[3]]), paste0("V.", dimnames(v.35ind.sparse.geomorph)[[3]]))

#GPA align
b.35ind.sparse.gpa <- ProcGPA(b.35ind.sparse.array, scale = FALSE, reflection = FALSE)

#Get info from the dimnames
info <- unlist(strsplit(dimnames(b.35ind.sparse.array)[[3]], split = "\\."))

#Camera
Camera <- c(rep("MD", times = dim(md.35ind.sparse.geomorph)[3]), rep("V", times = dim(v.35ind.sparse.geomorph)[3]))

#Individual
Individual <- info[grep("P", info)]

#Replicate
Replicate <- info[grep("R[1-3]", info)]

#Bind together into dataframe
b.35ind.sparse.covar <- as.data.frame(cbind(Camera, Individual, Replicate))
rm(Camera, Individual, Replicate, info)

#Create rrpp data frame
b.35ind.sparse.rrpp <- rrpp.data.frame(coords = two.d.array(b.35ind.sparse.gpa$rotated), Camera = b.35ind.sparse.covar$Camera, Individual = b.35ind.sparse.covar$Individual, Replicate = b.35ind.sparse.covar$Replicate)

#RCBD model with subsampling
b.35ind.sparse.rrpp.III <- lm.rrpp(coords ~ Camera*Individual, data = b.35ind.sparse.rrpp, print.progress = FALSE, iter = 99, SS.type = "III")

#Write anova results to file
capture.output(anova(b.35ind.sparse.rrpp.III, error = c("Camera:Individual", "Camera:Individual", "Residuals")), file = "FiguresAndTables/b.35ind.sparse.camerror.anova.txt")
