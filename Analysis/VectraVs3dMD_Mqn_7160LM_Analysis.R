#Analysis comparing 3dMDface to Vectra H1 using mannequin dense quasi-landmark data
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

##### Load the dense quasi-landmark mannequin data from MeshMonk #####
md.mqn <- read_csv("RawImageData/3DMD/Mannequin_WithPowder_Mapped/XYZ_Mqn_3DMD_7160LM.txt", col_names = paste0(c("x", "y", "z"), rep(1:7160, each=3)), col_types = cols(.default = col_double()))
v.mqn <- read_csv("RawImageData/Vectra/Mannequin_WithPowder_Mapped/XYZ_Mqn_Vectra_7160LM.txt", col_names = paste0(c("x", "y", "z"), rep(1:7160, each=3)), col_types = cols(.default = col_double()))

#Add columns to the dataframe as grouping variables
md.mqn$Replicate <- factor(c(rep(c("R1", "R2", "R3"), times = 3)))
md.mqn$Mapping <- factor(c(rep(c("M1", "M2", "M3"), each = 3)))

v.mqn$Replicate <- factor(c(rep(c("R1", "R2", "R3"), times = 3)))
v.mqn$Mapping <- factor(c(rep(c("M1", "M2", "M3"), each = 3)))

##### Within-camera registration error #####
#Euclidean distance between mapping centroid and each mapping
#3dMD
#Calculate the centroid of each replicate image (across 3 mappings)
md.mqn.mapping.centroid <- md.mqn %>%
  group_by(Replicate) %>%
  dplyr::select(-Mapping) %>%
  summarise_all(.funs = mean) %>%
  {. ->> tmp} %>%
  ungroup() %>%
  dplyr::select(-Replicate) %>%
  arrayspecs(., p=7160, k=3)

#Name the items of the array with Individual and Replicate info
dimnames(md.mqn.mapping.centroid)[[3]] <- tmp$Replicate

#Remove clutter
rm(tmp)

#Reorder the mappings by replicate and then turn into 3d array
md.mqn.array <- md.mqn %>%
  group_by(Replicate) %>%
  arrange(.by_group = TRUE) %>%
  {. ->> tmp} %>%
  ungroup () %>%
  dplyr::select(-c(Replicate, Mapping)) %>%
  arrayspecs(., p=7160, k=3)

#Name with Individual.Replicate.Mapping
dimnames(md.mqn.array)[[3]] <- paste(tmp$Replicate, tmp$Mapping, sep = ".")

#Remove clutter
rm(tmp)

#Calculate the euclidean distance between M1, M2, M3 and M.Centroid
#Create a matrix to store the euclidean distances
md.mqn.regerror.euclid <- matrix(nrow=9, ncol=7160)

#Create an index with numbers we'll need to extract the array elements
indx <- seq(1,9,3)

#Calculate the distance between the mappings and the centroid of the mappings, for each landmark
for (i in 1:length(indx)){
  #Distance between P#.R1 and P#.R1.M1
  md.mqn.regerror.euclid[indx[i],] <- diag(dist(md.mqn.mapping.centroid[,,i], md.mqn.array[,,indx[i]], method = "euclidean"))
  
  #Distance between P#.R1 and P#.R1.M2
  md.mqn.regerror.euclid[indx[i]+1,] <- diag(dist(md.mqn.mapping.centroid[,,i], md.mqn.array[,,indx[i]+1], method = "euclidean"))
  
  #Distance between P#.R1 and P#.R1.M3
  md.mqn.regerror.euclid[indx[i]+2,] <- diag(dist(md.mqn.mapping.centroid[,,i], md.mqn.array[,,indx[i]+2], method = "euclidean"))
}

#Convert the euclidean distances to a dataframe
md.mqn.regerror.euclid <- as.data.frame(md.mqn.regerror.euclid)

#Create a list of the naming information so that we can add it to the distance dataframe. We're relying on everything being in the same order to place these
names <- unlist(strsplit(dimnames(md.mqn.array)[[3]], split = "\\."))
md.mqn.regerror.euclid$Replicate <- names[grep(pattern="R", x=names)]
md.mqn.regerror.euclid$Mapping <- names[grep(pattern="M", x=names)]

#Remove clutter
rm(indx, i, names)

#Calculate the average euclidean distance per replicate
md.mqn.regerror.euclid.mean <- md.mqn.regerror.euclid %>%
  group_by(Replicate) %>%
  dplyr::select(-Mapping) %>%
  summarise_all(.funs = mean)

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(md.mqn.regerror.euclid.mean, select=-c(Replicate))))

#Vectra
#Calculate the centroid of each replicate image (across 3 mappings)
v.mqn.mapping.centroid <- v.mqn %>%
  group_by(Replicate) %>%
  dplyr::select(-Mapping) %>%
  summarise_all(.funs = mean) %>%
  {. ->> tmp} %>%
  ungroup() %>%
  dplyr::select(-Replicate) %>%
  arrayspecs(., p=7160, k=3)

#Name the items of the array with Individual and Replicate info
dimnames(v.mqn.mapping.centroid)[[3]] <- tmp$Replicate

#Remove clutter
rm(tmp)

#Reorder the mappings by replicate and then turn into 3d array
v.mqn.array <- v.mqn %>%
  group_by(Replicate) %>%
  arrange(.by_group = TRUE) %>%
  {. ->> tmp} %>%
  ungroup () %>%
  dplyr::select(-c(Replicate, Mapping)) %>%
  arrayspecs(., p=7160, k=3)

#Name with Individual.Replicate.Mapping
dimnames(v.mqn.array)[[3]] <- paste(tmp$Replicate, tmp$Mapping, sep = ".")

#Remove clutter
rm(tmp)

#Calculate the euclidean distance between M1, M2, M3 and M.Centroid
#Create a matrix to store the euclidean distances
v.mqn.regerror.euclid <- matrix(nrow=9, ncol=7160)

#Create an index with numbers we'll need to extract the array elements
indx <- seq(1,9,3)

#Calculate the distance between the mappings and the centroid of the mappings, for each landmark
for (i in 1:length(indx)){
  #Distance between P#.R1 and P#.R1.M1
  v.mqn.regerror.euclid[indx[i],] <- diag(dist(v.mqn.mapping.centroid[,,i], v.mqn.array[,,indx[i]], method = "euclidean"))
  
  #Distance between P#.R1 and P#.R1.M2
  v.mqn.regerror.euclid[indx[i]+1,] <- diag(dist(v.mqn.mapping.centroid[,,i], v.mqn.array[,,indx[i]+1], method = "euclidean"))
  
  #Distance between P#.R1 and P#.R1.M3
  v.mqn.regerror.euclid[indx[i]+2,] <- diag(dist(v.mqn.mapping.centroid[,,i], v.mqn.array[,,indx[i]+2], method = "euclidean"))
}

#Convert the euclidean distances to a dataframe
v.mqn.regerror.euclid <- as.data.frame(v.mqn.regerror.euclid)

#Create a list of the naming information so that we can add it to the distance dataframe. We're relying on everything being in the same order to place these
names <- unlist(strsplit(dimnames(v.mqn.array)[[3]], split = "\\."))
v.mqn.regerror.euclid$Replicate <- names[grep(pattern="R", x=names)]
v.mqn.regerror.euclid$Mapping <- names[grep(pattern="M", x=names)]

#Remove clutter
rm(indx, i, names)

#Calculate the average euclidean distance per replicate
v.mqn.regerror.euclid.mean <- v.mqn.regerror.euclid %>%
  group_by(Replicate) %>%
  dplyr::select(-Mapping) %>%
  summarise_all(.funs = mean)

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(v.mqn.regerror.euclid.mean, select=-c(Replicate))))

#Compare 3dMD and Vectra
#Stack the two dataframes of sd statistics on top of each other
b.mqn.regerror.euclid.mean <- rbind(md.mqn.regerror.euclid.mean, v.mqn.regerror.euclid.mean)

#Add a column telling which camera each came from
b.mqn.regerror.euclid.mean$Camera <- c(rep("3dMD", times = nrow(md.mqn.regerror.euclid.mean)), rep("Vectra", times = nrow(v.mqn.regerror.euclid.mean)))

#Plot
p <- ggplot(melt(b.mqn.regerror.euclid.mean, id.vars = c("Replicate", "Camera")), aes(x = Replicate, y=value, fill=Camera))+scale_fill_manual(values = c("#0D646B", "#C75E24"))+geom_boxplot(outlier.size = 0.25)+theme_bw()+ylab("Euclidean distance (mm)")+theme(legend.position = "bottom")
ggsave(filename = "FiguresAndTables/b.mqn.regerror.euclid.mean.pdf", plot = p, device = "pdf", width = 6, height = 6, units = "in")

##### Within-camera technological error #####
#Here we'll calculate the differences between the three different replicate images of the mannequin. 
#We calculated the centroid of each replicate image by averaging together the coordinates of the three mappings already. Now we want to GPA align the three replicate mannequin images, then extract the aligned coordinates. 
#3dMD
#GPA align the images
md.mqn.mapping.centroid.gpa <- ProcGPA(md.mqn.mapping.centroid, scale = FALSE, reflection = FALSE)

#Extract the aligned coordinates
md.mqn.mapping.centroid.aligned <- md.mqn.mapping.centroid.gpa$rotated

#Extract the mean shape
md.mqn.replicate.centroid <- md.mqn.mapping.centroid.gpa$mshape

#Remove the clutter
rm(md.mqn.mapping.centroid.gpa)

#Vectra
#GPA align the images
v.mqn.mapping.centroid.gpa <- ProcGPA(v.mqn.mapping.centroid, scale = FALSE, reflection = FALSE)

#Extract the aligned coordinates
v.mqn.mapping.centroid.aligned <- v.mqn.mapping.centroid.gpa$rotated

#Extract the mean shape
v.mqn.replicate.centroid <- v.mqn.mapping.centroid.gpa$mshape

#Remove the clutter
rm(v.mqn.mapping.centroid.gpa)

#Calculate the euclidean distance between R1, R2, R3 and Replicate.Centroid
#3dMD
#Create a matrix to store the euclidean distances
md.mqn.techerror.euclid <- matrix(nrow=3, ncol=7160)

#Distance between R1 and centroid
md.mqn.techerror.euclid[1,] <- diag(dist(md.mqn.replicate.centroid, md.mqn.mapping.centroid.aligned[,,1], method = "euclidean"))

#Distance between R2 and centroid
md.mqn.techerror.euclid[2,] <- diag(dist(md.mqn.replicate.centroid, md.mqn.mapping.centroid.aligned[,,2], method = "euclidean"))

md.mqn.techerror.euclid[3,] <- diag(dist(md.mqn.replicate.centroid, md.mqn.mapping.centroid.aligned[,,3], method = "euclidean"))

#Convert the euclidean distances to a dataframe
md.mqn.techerror.euclid <- as.data.frame(md.mqn.techerror.euclid)

#Add names
md.mqn.techerror.euclid$Replicate <- dimnames(md.mqn.mapping.centroid.aligned)[[3]]

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(md.mqn.techerror.euclid, select=-Replicate)))

#Vectra
#Create a matrix to store the euclidean distances
v.mqn.techerror.euclid <- matrix(nrow=3, ncol=7160)

#Distance between R1 and centroid
v.mqn.techerror.euclid[1,] <- diag(dist(v.mqn.replicate.centroid, v.mqn.mapping.centroid.aligned[,,1], method = "euclidean"))

#Distance between R2 and centroid
v.mqn.techerror.euclid[2,] <- diag(dist(v.mqn.replicate.centroid, v.mqn.mapping.centroid.aligned[,,2], method = "euclidean"))

v.mqn.techerror.euclid[3,] <- diag(dist(v.mqn.replicate.centroid, v.mqn.mapping.centroid.aligned[,,3], method = "euclidean"))

#Convert the euclidean distances to a dataframe
v.mqn.techerror.euclid <- as.data.frame(v.mqn.techerror.euclid)

#Add names
v.mqn.techerror.euclid$Replicate <- dimnames(v.mqn.mapping.centroid.aligned)[[3]]

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = colMeans(subset(v.mqn.techerror.euclid, select=-Replicate)))

#Compare 3dMD and Vectra
#Stack the two dataframes of sd statistics on top of each other
b.mqn.techerror.euclid <- rbind(md.mqn.techerror.euclid, v.mqn.techerror.euclid)

#Add a column telling which camera each came from
b.mqn.techerror.euclid$Camera <- c(rep("3dMD", times = nrow(md.mqn.techerror.euclid)), rep("Vectra", times = nrow(v.mqn.techerror.euclid)))

#Plot
p <- ggplot(melt(b.mqn.techerror.euclid, id.vars = c("Replicate", "Camera")), aes(x = Replicate, y=value, fill=Camera))+scale_fill_manual(values = c("#0D646B", "#C75E24"))+geom_boxplot(outlier.size = 0.25)+theme_bw()+ylab("Euclidean distance (mm)")+theme(legend.position = "bottom")
ggsave(filename = "FiguresAndTables/b.mqn.techerror.euclid.mean.pdf", plot = p, device = "pdf", width = 6.5, height = 5, units = "in")

##### Camera error #####
#Align the two images from each camera
#Put the 3dMD and Vectra average image together into an array
b.mqn.replicate.centroid <- abind(md.mqn.replicate.centroid, v.mqn.replicate.centroid, along = 3)

#Assign names so that we know set of coordinates came from which camera
dimnames(b.mqn.replicate.centroid)[[3]] <- c("3dMD", "Vectra")

#GPA align 
b.mqn.replicate.centroid.gpa <- ProcGPA(b.mqn.replicate.centroid, scale = FALSE, reflection = FALSE)

#Pull out the aligned coordinates into a 3d array
b.mqn.replicate.centroid.aligned <- b.mqn.replicate.centroid.gpa$rotated

#Remove clutter
rm(b.mqn.replicate.centroid.gpa)

#Calculate the euclidean distance between the aligned coordinates from each camera
b.mqn.camerror.euclid <- diag(dist(b.mqn.replicate.centroid.aligned[,,1], b.mqn.replicate.centroid.aligned[,,2], method = "euclidean"))

#Plot on a face
Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = b.mqn.camerror.euclid)

#Normal displacement
#We used Matlab functions to calculate normal displacement, so this part of the code is to write PLY files for reading into Matlab
#Ply header information, the same in all files.
Ply.Header <- "ply\nformat ascii 1.0\nelement vertex 7160\nproperty float x\nproperty float y\nproperty float z\nelement face 14050\nproperty list uchar int vertex_index\nend_header"
Faces <- cbind(matrix(3, nrow = 14050, ncol = 1), RefScan_Facets)

#Make the faces zero based instead of 1 based
FacesZero <- cbind(Faces[,1], apply(Faces[,-1], 2, function(x) x-1))

#3dMD
writeLines(Ply.Header, "NormalDisplacement/Mqn.MD.Powder.Zero.ply")
write.table(b.mqn.replicate.centroid.aligned[,,1], "NormalDisplacement/Mqn.MD.Powder.Zero.ply", sep = " ", col.names = F, row.names = F, quote = F, append = T)
write.table(FacesZero, "NormalDisplacement/Mqn.MD.Powder.Zero.ply", sep = " ", col.names = F, row.names = F, quote = F, append = T)

#Vectra
writeLines(Ply.Header, "NormalDisplacement/Mqn.V.Powder.Zero.ply")
write.table(b.mqn.replicate.centroid.aligned[,,2], "NormalDisplacement/Mqn.V.Powder.Zero.ply", sep = " ", col.names = F, row.names = F, quote = F, append = T)
write.table(FacesZero, "NormalDisplacement/Mqn.V.Powder.Zero.ply", sep = " ", col.names = F, row.names = F, quote = F, append = T)

#Read in information from Matlab and plot
#b.mqn.camerror.diststats <- read.table("NormalDisplacement/Mqn_Powder_3dMDvsVectraStats.txt", header = T, sep = ",", colClasses = "numeric")

#Plot on a face
#Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = b.mqn.camerror.diststats$NormalDistances, title = "Normal displacement from 3dMD to Vectra\n Red: Vectra front, Blue: 3dMD front")

##### ANOVA #####
#Put the 3dMD and Vectra images together in a single array
b.mqn.array <- abind(md.mqn.array, v.mqn.array, along = 3)

#Name the items
dimnames(b.mqn.array)[[3]] <- c(paste0("MD.", dimnames(md.mqn.array)[[3]]), paste0("V.", dimnames(v.mqn.array)[[3]]))

#GPA align
b.mqn.gpa <- ProcGPA(b.mqn.array, scale = FALSE, reflection = FALSE)

#Get info from the dimnames
info <- unlist(strsplit(dimnames(b.mqn.array)[[3]], split = "\\."))

#Camera
Camera <- c(rep("MD", times = dim(md.mqn.array)[3]), rep("V", times = dim(v.mqn.array)[3]))

#Replicate
#We don't want the replicates to be grouped together (i.e. R1.Vectra and R1.3dMD) so going to give them unique names
Replicate <- rep(c("R1", "R2", "R3", "R4", "R5", "R6"), each = 3)

#Bind together into dataframe
b.mqn.covar <- as.data.frame(cbind(Camera, Replicate))
rm(Camera, Replicate, info)

#Create rrpp data frame
b.mqn.rrpp <- rrpp.data.frame(coords = two.d.array(b.mqn.gpa$rotated), Camera = b.mqn.covar$Camera, Replicate = b.mqn.covar$Replicate)

#Mixed model nested anova 
b.mqn.rrpp.III <- lm.rrpp(coords ~ Camera*Replicate, data = b.mqn.rrpp, print.progress = FALSE, iter = 99, SS.type = "III")

#Write anova results to file, replicate is random effect
capture.output(anova(b.mqn.rrpp.III, error = c("Residuals", "Residuals")), file = "FiguresAndTables/b.mqn.camerror.anova.txt")

