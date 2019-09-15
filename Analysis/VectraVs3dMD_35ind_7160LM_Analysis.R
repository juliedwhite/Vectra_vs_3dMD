#Analysis comparing 3dMDface to Vectra H1 using participant dense quasi-landmark data
#Author: Julie D. White

#We are restricted in sharing the dense-landmark coordinates or raw image files for the 35 participants, so this code is just for reference. 

#Necessary packages
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

#Load the dense quasi-landmarked participant data from MeshMonk
md.35ind <- read_csv("RawImageData/XYZ_35ind_3DMD_7160LM.txt", col_names = paste0(c("x", "y", "z"), rep(1:7160, each=3)), col_types = cols(.default = col_double()))
v.35ind <- read_csv("RawImageData/XYZ_35ind_Vectra_7160LM.txt", col_names = paste0(c("x", "y", "z"), rep(1:7160, each=3)), col_types = cols(.default = col_double()))

#Add columns to the dataframe as grouping variables
md.35ind$Individual <- factor(c(rep(paste0("P", rep(1:35, each=3)), times=3)), levels = c(paste0("P", 1:35)))
md.35ind$Replicate <- factor(c(rep(c("R1", "R2", "R3"), length = 315)))
md.35ind$Mapping <- factor(c(rep(c("M1", "M2", "M3"), each = 105)))

v.35ind$Individual <- factor(c(rep(paste0("P", rep(1:35, each=3)), times=3)), levels = c(paste0("P", 1:35)))
v.35ind$Replicate <- factor(c(rep(c("R1", "R2", "R3"), length = 315)))
v.35ind$Mapping <- factor(c(rep(c("M1", "M2", "M3"), each = 105)))

#Within-camera registration error
#Euclidean distance between mapping centroid and each mapping
#3dMD
#Calculate the centroid of each replicate image (across 3 mappings)
md.35ind.mapping.centroid <- md.35ind %>%
  group_by(Individual, Replicate) %>%
  dplyr::select(-Mapping) %>%
  summarise_all(.funs = mean) %>%
  {. ->> tmp} %>%
  ungroup() %>%
  dplyr::select(-c(Individual, Replicate)) %>%
  arrayspecs(., p=7160, k=3)

#Name the third dimension with Individual and Replicate info
dimnames(md.35ind.mapping.centroid)[[3]] <- paste(tmp$Individual, tmp$Replicate, sep=".")

#Remove clutter
rm(tmp)

#Reorder the mappings by individual and then replicate and then turn into 3d array
md.35ind.array <- md.35ind %>%
  group_by(Individual, Replicate) %>%
  arrange(.by_group = TRUE) %>%
  {. ->> tmp} %>%
  ungroup () %>%
  dplyr::select(-c(Individual, Replicate, Mapping)) %>%
  arrayspecs(., p=7160, k=3)

#Name with Individual.Replicate.Mapping
dimnames(md.35ind.array)[[3]] <- paste(tmp$Individual, tmp$Replicate, tmp$Mapping, sep = ".")

#Remove clutter
rm(tmp)

#Calculate the euclidean distance between M1, M2, M3 and M.Centroid
#Create a matrix to store the euclidean distances
md.35ind.regerror.euclid <- matrix(nrow=315, ncol=7160)

#Create an index with numbers we'll need to extract the array elements
indx <- seq(1,315,3)

#Calculate the distance between the mappings and the centroid of the mappings, for each landmark
for (i in 1:length(indx)){
  #Distance between P#.R1 and P#.R1.M1
  md.35ind.regerror.euclid[indx[i],] <- diag(dist(md.35ind.mapping.centroid[,,i], md.35ind.array[,,indx[i]], method = "euclidean"))
  
  #Distance between P#.R1 and P#.R1.M2
  md.35ind.regerror.euclid[indx[i]+1,] <- diag(dist(md.35ind.mapping.centroid[,,i], md.35ind.array[,,indx[i]+1], method = "euclidean"))
  
  #Distance between P#.R1 and P#.R1.M3
  md.35ind.regerror.euclid[indx[i]+2,] <- diag(dist(md.35ind.mapping.centroid[,,i], md.35ind.array[,,indx[i]+2], method = "euclidean"))
}

#Convert the euclidean distances to a dataframe
md.35ind.regerror.euclid <- as.data.frame(md.35ind.regerror.euclid)

#Create a list of the naming information so that we can add it to the distance dataframe. We're relying on everything being in the same order to place these
names <- unlist(strsplit(dimnames(md.35ind.array)[[3]], split = "\\."))
md.35ind.regerror.euclid$Individual <- names[grep(pattern="P", x=names)]
md.35ind.regerror.euclid$Replicate <- names[grep(pattern="R", x=names)]
md.35ind.regerror.euclid$Mapping <- names[grep(pattern="M", x=names)]

#Remove clutter
rm(indx, i, names)

#Calculate the average euclidean distance per replicate
md.35ind.regerror.euclid.mean <- md.35ind.regerror.euclid %>%
  group_by(Individual, Replicate) %>%
  dplyr::select(-Mapping) %>%
  summarise_all(.funs = mean)

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(md.35ind.regerror.euclid.mean, select=-c(Individual,Replicate))))

#Vectra
#Calculate the centroid of each replicate image (across 3 mappings)
v.35ind.mapping.centroid <- v.35ind %>%
  group_by(Individual, Replicate) %>%
  dplyr::select(-Mapping) %>%
  summarise_all(.funs = mean) %>%
  {. ->> tmp} %>%
  ungroup() %>%
  dplyr::select(-c(Individual, Replicate)) %>%
  arrayspecs(., p=7160, k=3)

#Name the third dimension with Individual and Replicate info
dimnames(v.35ind.mapping.centroid)[[3]] <- paste(tmp$Individual, tmp$Replicate, sep=".")

#Remove clutter
rm(tmp)

#Reorder the mappings by individual and then replicate and then turn into 3d array
v.35ind.array <- v.35ind %>%
  group_by(Individual, Replicate) %>%
  arrange(.by_group = TRUE) %>%
  {. ->> tmp} %>%
  ungroup () %>%
  dplyr::select(-c(Individual, Replicate, Mapping)) %>%
  arrayspecs(., p=7160, k=3)

#Name with Individual.Replicate.Mapping
dimnames(v.35ind.array)[[3]] <- paste(tmp$Individual, tmp$Replicate, tmp$Mapping, sep = ".")

#Remove clutter
rm(tmp)

#Calculate the euclidean distance between M1, M2, M3 and M.Centroid
#Create a matrix to store the euclidean distances
v.35ind.regerror.euclid <- matrix(nrow=315, ncol=7160)

#Create an index with numbers we'll need to extract the array elements
indx <- seq(1,315,3)

#Calculate the distance between the mappings and the centroid of the mappings, for each landmark
for (i in 1:length(indx)){
  #Distance between P#.R1 and P#.R1.M1
  v.35ind.regerror.euclid[indx[i],] <- diag(dist(v.35ind.mapping.centroid[,,i], v.35ind.array[,,indx[i]], method = "euclidean"))
  
  #Distance between P#.R1 and P#.R1.M2
  v.35ind.regerror.euclid[indx[i]+1,] <- diag(dist(v.35ind.mapping.centroid[,,i], v.35ind.array[,,indx[i]+1], method = "euclidean"))
  
  #Distance between P#.R1 and P#.R1.M3
  v.35ind.regerror.euclid[indx[i]+2,] <- diag(dist(v.35ind.mapping.centroid[,,i], v.35ind.array[,,indx[i]+2], method = "euclidean"))
}

#Convert the euclidean distances to a dataframe
v.35ind.regerror.euclid <- as.data.frame(v.35ind.regerror.euclid)

#Create a list of the naming information so that we can add it to the distance dataframe. We're relying on everything being in the same order to place these
names <- unlist(strsplit(dimnames(v.35ind.array)[[3]], split = "\\."))
v.35ind.regerror.euclid$Individual <- names[grep(pattern="P", x=names)]
v.35ind.regerror.euclid$Replicate <- names[grep(pattern="R", x=names)]
v.35ind.regerror.euclid$Mapping <- names[grep(pattern="M", x=names)]

#Remove clutter
rm(indx, i, names)

#Calculate the average euclidean distance per replicate
v.35ind.regerror.euclid.mean <- v.35ind.regerror.euclid %>%
  group_by(Individual, Replicate) %>%
  dplyr::select(-Mapping) %>%
  summarise_all(.funs = mean)

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(v.35ind.regerror.euclid.mean, select=-c(Individual,Replicate))))

#Compare 3dMD and Vectra
#Stack the two dataframes of sd statistics on top of each other
b.35ind.regerror.euclid.mean <- rbind(md.35ind.regerror.euclid.mean, v.35ind.regerror.euclid.mean)

#Add a column telling which camera each came from
b.35ind.regerror.euclid.mean$Camera <- c(rep("3dMD", times = nrow(md.35ind.regerror.euclid.mean)), rep("Vectra", times = nrow(v.35ind.regerror.euclid.mean)))

#Reorder the factor levels by number
b.35ind.regerror.euclid.mean$Individual <- factor(b.35ind.regerror.euclid.mean$Individual, levels = c(paste0("P", 1:35)))

#Plot with outliers
ggplot(melt(b.35ind.regerror.euclid.mean, id.vars = c("Individual", "Replicate", "Camera")), aes(x = Individual, y=value, fill=Camera))+ylab("Euclidean distance (mm)")+geom_boxplot(outlier.size = 0.25)+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+scale_fill_manual(values = c("#0D646B", "#C75E24"))+scale_color_manual(values = c("#0D646B", "#C75E24"))

#Plot per camera instead of per person, and perform Wilcoxon test on mean comparison
tmp <- b.35ind.regerror.euclid.mean %>% 
  group_by(Individual, Replicate, Camera) %>%
  ungroup() %>%
  mutate(Avg = rowMeans(.[,-c(which(colnames(b.35ind.regerror.euclid.mean) %in% c("Individual", "Replicate", "Camera")))])) %>%
  dplyr::select(Individual, Replicate, Camera, Avg)

ggplot(tmp, aes(x = Camera, y=Avg, fill = Camera))+geom_point(aes(x = Camera, y = Avg), size = 0.25)+geom_boxplot(outlier.alpha = 0, coef = 0, alpha = 0.75)+scale_fill_manual(values = c("#0D646B", "#C75E24"))+ylab("Euclidean distance (mm)")+xlab("Camera")+theme_bw()+theme(legend.position = "none", axis.title = element_text())+stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 2)

#Remove clutter
rm(tmp)

#Within-camera "biological" error
#Here we'll calculate the differences between the three different replicate images of the same person. 
#We calculated the centroid of each replicate image by averaging together the coordinates of the three mappings already. Now we want to pull out the three replicate images for each person and GPA align those, then extract the aligned coordinates. 
#3dMD
#Make a list containing the indices for each indiviudal. Should be a list of 35 with three items each. 
ind <- paste0("P", 1:35, "\\.")
ind.index <- lapply(ind, function(x) grep(x, dimnames(md.35ind.mapping.centroid)[[3]]))

#3d array to store the GPA aligned coordinates
md.35ind.mapping.centroid.aligned <- NULL

#3d array to store the mshape, which is the centroid
md.35ind.replicate.centroid <- NULL

#List to store the names so we dont get lost
namelist <- NULL

#For each item in ind.index, pull out the three replicates for each image and convert to a 3d array. 
for (i in 1:length(ind.index)){
  #Pull out the three matrices that correspond to a single person
  tmp <- md.35ind.mapping.centroid[,,ind.index[[i]]]
  #GPA align without scaling or reflecting
  tmp.gpa <- ProcGPA(tmp, scale = FALSE, reflection = FALSE)
  #Pull out the aligned coordinates and place them to a 3d array
  md.35ind.mapping.centroid.aligned <- abind(md.35ind.mapping.centroid.aligned, tmp.gpa$rotated, along = 3)
  #Pull out the dimnames and add them to a list of names
  namelist <- append(namelist, dimnames(tmp.gpa$rotated)[[3]])
  #Pull out the mean shape of the three replicate images and put into 3d array 
  md.35ind.replicate.centroid <- abind(md.35ind.replicate.centroid, tmp.gpa$mshape, along = 3)
}

#Assign names to the centroid of the replicates, using a pared down version of the namelist. The names in this list are in the same order as the data because they were constructed together. 
dimnames(md.35ind.replicate.centroid)[[3]] <- unique(unlist(strsplit(namelist, split = ".R[0-3]")))

#Remove the clutter
rm(ind, ind.index, namelist, tmp, tmp.gpa, i)

#Vectra
#Make a list containing the indices for each indiviudal. Should be a list of 35 with three items each. 
ind <- paste0("P", 1:35, "\\.")
ind.index <- lapply(ind, function(x) grep(x, dimnames(v.35ind.mapping.centroid)[[3]]))

#3d array to store the GPA aligned coordinates
v.35ind.mapping.centroid.aligned <- NULL

#3d array to store the mshape, which is the centroid of 
v.35ind.replicate.centroid <- NULL

#List to store the names so we dont get lost
namelist <- NULL

#For each item in the ind.index list, pull out the three replicates for each image and convert to a 3d array. 
for (i in 1:length(ind.index)){
  #Pull out the three matrices that correspond to a single person
  tmp <- v.35ind.mapping.centroid[,,ind.index[[i]]]
  #GPA align without scaling or reflecting
  tmp.gpa <- ProcGPA(tmp, scale = FALSE, reflection = FALSE)
  #Pull out the aligned coordinates and place them to a 3d array
  v.35ind.mapping.centroid.aligned <- abind(v.35ind.mapping.centroid.aligned, tmp.gpa$rotated, along = 3)
  #Pull out the dimnames and add them to a list of names
  namelist <- append(namelist, dimnames(tmp.gpa$rotated)[[3]])
  #Pull out the mean shape of the three replicate images and put into 3d array 
  v.35ind.replicate.centroid <- abind(v.35ind.replicate.centroid, tmp.gpa$mshape, along = 3)
}

#Assign names to the centroid of the replicates, using a pared down version of the namelist. The names in this list are in the same order as the data because they were constructed together. 
dimnames(v.35ind.replicate.centroid)[[3]] <- unique(unlist(strsplit(namelist, split = ".R[0-3]")))

#Remove the clutter
rm(ind, ind.index, namelist, tmp, tmp.gpa, i)

#Euclidean distance between each replicate image and the centroid of the replicates
#3dMD
#Create a matrix to store the euclidean distances
md.35ind.bioerror.euclid <- matrix(nrow=105, ncol=7160)

#Create an index with numbers we'll need to extract the list elements
indx <- seq(1,105,3)

#Calculate the distance between the replicates and the centroid of the replicates, for each landmark
for (i in 1:length(indx)){
  #Distance between P#.R1 and P#
  md.35ind.bioerror.euclid[indx[i],] <- diag(dist(md.35ind.replicate.centroid[,,i], md.35ind.mapping.centroid.aligned[,,indx[i]], method = "euclidean"))
  
  #Distance between P#.R2 and P#
  md.35ind.bioerror.euclid[indx[i]+1,] <- diag(dist(md.35ind.replicate.centroid[,,i], md.35ind.mapping.centroid.aligned[,,indx[i]+1], method = "euclidean"))
  
  #Distance between P#.R3 and P#
  md.35ind.bioerror.euclid[indx[i]+2,] <- diag(dist(md.35ind.replicate.centroid[,,i], md.35ind.mapping.centroid.aligned[,,indx[i]+2], method = "euclidean"))
}

#Convert the euclidean distances to a dataframe
md.35ind.bioerror.euclid <- as.data.frame(md.35ind.bioerror.euclid)

#Create a list of the naming information so that we can add it to the distance dataframe. We're relying on everything being in the same order to place these
names <- unlist(strsplit(dimnames(md.35ind.mapping.centroid.aligned)[[3]], split = "\\."))
md.35ind.bioerror.euclid$Individual <- names[grep(pattern="P", x=names)]
md.35ind.bioerror.euclid$Replicate <- names[grep(pattern="R", x=names)]

#Remove clutter
rm(indx, i, names)

#Calculate the average euclidean distance per replicate
md.35ind.bioerror.euclid.mean <- md.35ind.bioerror.euclid %>%
  group_by(Individual) %>%
  dplyr::select(-Replicate) %>%
  summarise_all(.funs = mean)

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(md.35ind.bioerror.euclid.mean, select=-Individual)))

#Vectra
#Create a matrix to store the euclidean distances
v.35ind.bioerror.euclid <- matrix(nrow=105, ncol=7160)

#Create an index with numbers we'll need to extract the list elements
indx <- seq(1,105,3)

#Calculate the distance between the replicates and the centroid of the replicates, for each landmark
for (i in 1:length(indx)){
  #Distance between P#.R1 and P#
  v.35ind.bioerror.euclid[indx[i],] <- diag(dist(v.35ind.replicate.centroid[,,i], v.35ind.mapping.centroid.aligned[,,indx[i]], method = "euclidean"))
  
  #Distance between P#.R2 and P#
  v.35ind.bioerror.euclid[indx[i]+1,] <- diag(dist(v.35ind.replicate.centroid[,,i], v.35ind.mapping.centroid.aligned[,,indx[i]+1], method = "euclidean"))
  
  #Distance between P#.R3 and P#
  v.35ind.bioerror.euclid[indx[i]+2,] <- diag(dist(v.35ind.replicate.centroid[,,i], v.35ind.mapping.centroid.aligned[,,indx[i]+2], method = "euclidean"))
}

#Convert the euclidean distances to a dataframe
v.35ind.bioerror.euclid <- as.data.frame(v.35ind.bioerror.euclid)

#Create a list of the naming information so that we can add it to the distance dataframe. We're relying on everything being in the same order to place these
names <- unlist(strsplit(dimnames(v.35ind.mapping.centroid.aligned)[[3]], split = "\\."))
v.35ind.bioerror.euclid$Individual <- names[grep(pattern="P", x=names)]
v.35ind.bioerror.euclid$Replicate <- names[grep(pattern="R", x=names)]

#Remove clutter
rm(indx, i, names)

#Calculate the average euclidean distance per replicate
v.35ind.bioerror.euclid.mean <- v.35ind.bioerror.euclid %>%
  group_by(Individual) %>%
  dplyr::select(-Replicate) %>%
  summarise_all(.funs = mean)

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(v.35ind.bioerror.euclid.mean, select=-Individual)))

#Compare 3dMD and Vectra
#Stack the two dataframes of sd statistics on top of each other
b.35ind.bioerror.euclid.mean <- rbind(md.35ind.bioerror.euclid.mean, v.35ind.bioerror.euclid.mean)

#Add a column telling which camera each came from
b.35ind.bioerror.euclid.mean$Camera <- c(rep("3dMD", times = nrow(md.35ind.bioerror.euclid.mean)), rep("Vectra", times = nrow(v.35ind.bioerror.euclid.mean)))

#Reorder the factor levels by number
b.35ind.bioerror.euclid.mean$Individual <- factor(b.35ind.bioerror.euclid.mean$Individual, levels = c(paste0("P", 1:35)))

#Plot with outliers
ggplot(melt(b.35ind.bioerror.euclid.mean, id.vars = c("Individual", "Camera")), aes(x = Individual, y=value, fill=Camera))+geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ylab("Euclidean distance (mm)")+scale_fill_manual(values = c("#0D646B", "#C75E24"))

#Plot per camera instead of per person
tmp <- b.35ind.bioerror.euclid.mean %>% 
  group_by(Individual, Camera) %>%
  ungroup() %>%
  mutate(Avg = rowMeans(.[,-c(which(colnames(b.35ind.bioerror.euclid.mean) %in% c("Individual", "Camera")))])) %>%
  dplyr::select(Individual, Camera, Avg)

ggplot(melt(tmp, id.vars = c("Individual", "Camera")), aes(x = Camera, y=value, fill = Camera))+geom_point(colour = "black", size = 0.9)+geom_line(aes(group = Individual), colour = "black", alpha = 0.25)+geom_boxplot(outlier.shape = NA, alpha = 0.75, coef = 0)+geom_boxplot(aes(color = Camera), fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0)+scale_fill_manual(values = c("#0D646B", "#C75E24"))+scale_color_manual(values = c("#0D646B", "#C75E24"))+ylab("Euclidean distance (mm)")+xlab("Camera")+theme_bw()+theme(legend.position = "none", axis.title = element_text())

#Remove clutter
rm(tmp)

##### Camera error #####
#Align the two images from each camera
#We calculated the centroid of each person by averaging together the coordinates of the three images already. Now we want to pull out the two images from each camera for each person and GPA align those, then extract the aligned coordinates. 

#Combine the replicate centroid arrays into a single array
b.35ind.replicate.centroid <- abind(md.35ind.replicate.centroid, v.35ind.replicate.centroid, along = 3)

#Alter the dimnames so that we know which camera each individual came from 
dimnames(b.35ind.replicate.centroid)[[3]] <- paste(dimnames(b.35ind.replicate.centroid)[[3]], rep(c("MD", "V"), each = 35), sep = ".")

#Make a list containing the indices for each indiviudal. Should be a list of 35 with two items each. 
ind <- paste0("P", 1:35, "\\.")
ind.index <- lapply(ind, function(x) grep(x, dimnames(b.35ind.replicate.centroid)[[3]]))

#3d array to store the GPA aligned coordinates
b.35ind.replicate.centroid.aligned <- NULL

#For each item in ind.index, pull out the two images ang align
for (i in 1:length(ind.index)){
  #Pull out the two matrices that correspond to a single person
  tmp <- b.35ind.replicate.centroid[,,ind.index[[i]]]
  #GPA align without scaling or reflecting
  tmp.gpa <- ProcGPA(tmp, scale = FALSE, reflection = FALSE)
  #Pull out the aligned coordinates and place them to a 3d array
  b.35ind.replicate.centroid.aligned <- abind(b.35ind.replicate.centroid.aligned, tmp.gpa$rotated, along = 3)
}

#Remove the clutter
rm(ind, ind.index, tmp, tmp.gpa, i)

#Calculate the euclidean distance between the aligned coordinates from each camera
#Make a list containing the indices for each indiviudal. Should be a list of 35 with two items each. 
ind <- paste0("P", 1:35, "\\.")
ind.index <- lapply(ind, function(x) grep(x, dimnames(b.35ind.replicate.centroid.aligned)[[3]]))

#Empty matrix to store the euclidean distances
b.35ind.camerror.euclid <- matrix(nrow = 35, ncol = 7160)

#For each of the 35 people, calculate the euclidean distance between their 3dMD and Vectra coordinates and place into matrix
for (i in 1:length(ind.index)){
  b.35ind.camerror.euclid[i,] <- diag(dist(b.35ind.replicate.centroid.aligned[,,ind.index[[i]][1]], b.35ind.replicate.centroid.aligned[,,ind.index[[i]][2]], method = "euclidean"))
}

#Convert result matrix to dataframe
b.35ind.camerror.euclid <- as.data.frame(b.35ind.camerror.euclid)

#Create a Individual column for grouping
b.35ind.camerror.euclid$Individual <- unique(unlist(strsplit(dimnames(b.35ind.replicate.centroid.aligned)[[3]], split = c(".MD", ".V"))))

#Remove clutter
rm(ind, ind.index, i)

#Plot per person
#Reorder the factor levels by number
b.35ind.camerror.euclid$Individual <- factor(b.35ind.camerror.euclid$Individual, levels = c(paste0("P", 1:35)))

#Plot
p <- ggplot(melt(b.35ind.camerror.euclid, id.vars = "Individual"), aes(x = Individual, y = value))+geom_boxplot(outlier.size = 0.25, fill = "#70B77E")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")+labs(x="Individual", y="Euclidean distance (mm)")
ggsave(filename = "FiguresAndTables/b.35ind.camerror.euclid.pdf", plot = p, device = "pdf", width = 6.5, height = 6, units = "in")

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(b.35ind.camerror.euclid, select=-Individual)))

#Normal displacement
#We used Matlab functions to calculate normal displacement, so this part of the code is to write PLY files for reading into Matlab
#Ply header information, the same in all files.
Ply.Header <- "ply\nformat ascii 1.0\nelement vertex 7160\nproperty float x\nproperty float y\nproperty float z\nelement face 14050\nproperty list uchar int vertex_index\nend_header"
Faces <- cbind(matrix(3, nrow = 14050, ncol = 1), RefScan_Facets)

#For each person, export a ply file with their 3d image information
#Create a list of filenames with the folder and name where we want the ply file stored, in the same order as the data are in the aligned array.
filenames <- paste0("NormalDisplacement/", dimnames(b.35ind.replicate.centroid.aligned)[[3]], ".ply")

for (i in 1:length(filenames)){
  Vertices <- b.35ind.replicate.centroid.aligned[,,i]
  writeLines(Ply.Header, filenames[i])
  write.table(Vertices, filenames[i], sep = " ", col.names = F, row.names = F, quote = F, append = T)
  write.table(Faces, filenames[i], sep = " ", col.names = F, row.names = F, quote = F, append = T)
}

#Read in the normal displacement statistics for each person, average, and plot
#List of files to be read in
files <- paste0("NormalDisplacement", paste0("P", 1:35), "_3dMDvsVectraStats.txt")

#Create an empty 3d matrix where we'll store the full distance information
b.35ind.camerror.diststats <- array(NA, c(7160, 1, 35))

#Read in files
for (i in 1:length(files)){
  b.35ind.camerror.diststats[,,i] <- as.matrix(read.table(files[i], header = T, sep = ",", colClasses = "numeric"))
}

tmp <- read.table(files[i], header = T, sep = ",", colClasses = "numeric")

#Add dimnames
dimnames(b.35ind.camerror.diststats)[[2]] <- colnames(tmp)
dimnames(b.35ind.camerror.diststats)[[3]] <- paste0("P", 1:35)

#Remove clutter
rm(files, i, tmp)

#Plot average normal displacement on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = rowMeans(b.35ind.camerror.diststats[,"NormalDistances",]), title = "Normal displacement from 3dMD to Vectra\n Red: Vectra front, Blue = 3dMD front")

##### ANOVA #####
#Put the 3dMD and Vectra images together in a single array
b.35ind.array <- abind(md.35ind.array, v.35ind.array, along = 3)

#Name the items
dimnames(b.35ind.array)[[3]] <- c(paste0("MD.", dimnames(md.35ind.array)[[3]]), paste0("V.", dimnames(v.35ind.array)[[3]]))

#GPA align
b.35ind.gpa <- ProcGPA(b.35ind.array, scale = FALSE, reflection = FALSE)

#Get info from the dimnames
info <- unlist(strsplit(dimnames(b.35ind.array)[[3]], split = "\\."))

#Camera
Camera <- c(rep("MD", times = dim(md.35ind.array)[3]), rep("V", times = dim(v.35ind.array)[3]))

#Individual
Individual <- info[grep("P", info)]

#Replicate
Replicate <- info[grep("R[1-3]", info)]

#Mapping
Mapping <- info[grep("M[1-3]", info)]

#Bind together into dataframe
b.35ind.covar <- as.data.frame(cbind(Camera, Individual, Replicate, Mapping))
rm(Camera, Individual, Replicate, Mapping, info)

#Create rrpp data frame
b.35ind.rrpp <- rrpp.data.frame(coords = two.d.array(b.35ind.gpa$rotated), Camera = b.35ind.covar$Camera, Individual = b.35ind.covar$Individual, Replicate = b.35ind.covar$Replicate)

#Nested RCBD design ANOVA
#Using Type III SS because the interaction between Camera:Individual is significant
b.35ind.rrpp.III <- lm.rrpp(coords ~ Camera+Individual+Camera:Individual+Camera:Individual:Replicate, data = b.35ind.rrpp, print.progress = FALSE, iter = 99, SS.type = "III")

#Treat individual and replicate as random effects
anova(b.35ind.rrpp.III, error = c("Camera:Individual", "Camera:Individual", "Camera:Individual:Replicate", "Residuals"))

