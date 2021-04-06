library(shapes)
library(rgl)
# library(matlib)
# library(RiemBase)
# library(mvtnorm)
# library(matrixcalc) 
# library(tictoc)
# library(caret)
# library(e1071)
# library("kdensity")

#####################################################################################################
#####################################################################################################
#set working directory to file location

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#####################################################################################################
#####################################################################################################
# functions 

source("subFunctions/MathFunctions.R")

#####################################################################################################
#####################################################################################################
# read s-rep data 

# read sreps
srepsDataG1<- read.csv(file=paste("../files/hippocampusNeg.csv",sep = ""), header=TRUE, sep=",")
#View(srepsDataG1)
srepsDataG2<- read.csv(file=paste("../files/hippocampusPos.csv",sep = ""), header=TRUE, sep=",")
#View(srepsDataG2)
srepsDataG3<- read.csv(file=paste("../files/caudateNeg.csv",sep = ""), header=TRUE, sep=",")
#View(srepsDataG3)
srepsDataG4<- read.csv(file=paste("../files/caudatePos.csv",sep = ""), header=TRUE, sep=",")
#View(srepsDataG4)

# or
nSamplesG1<-max(srepsDataG1$srepNumber)
nSamplesG2<-max(srepsDataG2$srepNumber)
nSamplesG3<-max(srepsDataG3$srepNumber)
nSamplesG4<-max(srepsDataG4$srepNumber)

nSamplesG1
nSamplesG2
nSamplesG3
nSamplesG4

#####################################################################################################
#####################################################################################################
# extract number of srep GOPs

#No. samples and spokes
upSpoeksNumber<-max(srepsDataG1$SpokesNumber[which(srepsDataG1$srepNumber==1 & srepsDataG1$Spoke=='up')])
downSpoeksNumber<-max(srepsDataG1$SpokesNumber[which(srepsDataG1$srepNumber==1 & srepsDataG1$Spoke=='down')])
crestSpoksNumber<-max(srepsDataG1$SpokesNumber[which(srepsDataG1$srepNumber==1 & srepsDataG1$Spoke=='crest')])
nTotalRadii <- upSpoeksNumber + downSpoeksNumber + crestSpoksNumber
skelPointNo <- nTotalRadii-downSpoeksNumber
skelRange<-c(1:downSpoeksNumber,(2*downSpoeksNumber+1):nTotalRadii)

#####################################################################################################
#####################################################################################################
# skeletal and Boundary PDM of G1 and G2

source("subFunctions/readSrepsData.R")

tempG1<-readSrepsData(srepsData = srepsDataG1)
SkeletalPDMG1<-tempG1$SkeletalPDM
BoundaryPDMG1<-tempG1$BoundaryPDM
boundaryPlusSkeletal_G1<-tempG1$boundaryPlusSkeletal

tempG2<-readSrepsData(srepsData = srepsDataG2)
SkeletalPDMG2<-tempG2$SkeletalPDM
BoundaryPDMG2<-tempG2$BoundaryPDM
boundaryPlusSkeletal_G2<-tempG2$boundaryPlusSkeletal

tempG3<-readSrepsData(srepsData = srepsDataG3)
SkeletalPDMG3<-tempG3$SkeletalPDM
BoundaryPDMG3<-tempG3$BoundaryPDM
boundaryPlusSkeletal_G3<-tempG3$boundaryPlusSkeletal

tempG4<-readSrepsData(srepsData = srepsDataG4)
SkeletalPDMG4<-tempG4$SkeletalPDM
BoundaryPDMG4<-tempG4$BoundaryPDM
boundaryPlusSkeletal_G4<-tempG4$boundaryPlusSkeletal

#####################################################################################################
#####################################################################################################
# plot

sampleNo<-1
srep1<-rbind(SkeletalPDMG1[,,sampleNo],BoundaryPDMG1[,,sampleNo])
plot3d(SkeletalPDMG1[skelRange,,sampleNo],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
for (i in 1:nTotalRadii) {
  plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 2,col = "blue",expand = 10,box=FALSE,add = TRUE)
}
srep3<-rbind(SkeletalPDMG3[,,sampleNo],BoundaryPDMG3[,,sampleNo])
plot3d(SkeletalPDMG3[skelRange,,sampleNo],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
for (i in 1:nTotalRadii) {
  plot3d(srep3[c(i,(i+nTotalRadii)),],type="l",lwd = 2,col = "blue",expand = 10,box=FALSE,add = TRUE)
}

plot3d(BoundaryPDMG1[,,sampleNo],type="s", radius = 0.1,col = "green",expand = 10,box=FALSE,add = TRUE)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)

#####################################################################################################
#####################################################################################################
# check for relabeling

#plot test
for (k in c(1,2,3,4)) {
  skelG1_1<-SkeletalPDMG1[skelRange,,k]
  # open3d()
  referencePoint<-skelG1_1[6,]
  plot3d(skelG1_1,type="p", size=0.5,col = "red",expand = 10,box=FALSE,add = TRUE)
  for (i in 1:(dim(skelG1_1)[1]-1)) {
    plot3d(skelG1_1[c(i,i+1),],type="l",lwd = 0.5,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  plot3d(rbind(referencePoint,referencePoint),type="s",size = 0.8,col = "blue",expand = 10,box=FALSE,add = TRUE)
}

#####################################################################################################
#####################################################################################################
# relabe by fixFilipedSreps function
source("subFunctions/fixFlipIssue.R")

referenceSrep<-srepsDataG1[srepsDataG1$srepNumber==1,]

# NB!!!! define the grid size as "5_7" or "5_9"
gridSize = "5_9"

if(gridSize == "5_7"){
  centerPoint<-16 #intrinsic centroid
  p2<-13
  p3<-17
}else if(gridSize == "5_9"){
  centerPoint<-19 #intrinsic centroid
  p2<-16
  p3<-20
}
v0<-BoundaryPDMG1[centerPoint,,1]-SkeletalPDMG1[centerPoint,,1]
reference_u0<-v0/sqrt(sum(v0^2))
v1<-SkeletalPDMG1[centerPoint,,1]-SkeletalPDMG1[p2,,1]
reference_u1<-v1/sqrt(sum(v1^2))
v2<-SkeletalPDMG1[centerPoint,,1]-SkeletalPDMG1[p3,,1]
reference_u2<-v2/sqrt(sum(v2^2))

sorted_srepDataAdjustedG1<-fixFilipedSreps(srepsData = srepsDataG1,
                                           reference_u0 = reference_u0,
                                           reference_u1 = reference_u1,
                                           reference_u2 = reference_u2,
                                           centerPoint = centerPoint, p2 = p2, p3 = p3,
                                           gridSize = gridSize)
#View(sorted_srepDataAdjustedG1)
sorted_srepDataAdjustedG2<-fixFilipedSreps(srepsData = srepsDataG2,
                                           reference_u0 = reference_u0,
                                           reference_u1 = reference_u1,
                                           reference_u2 = reference_u2,
                                           centerPoint = centerPoint, p2 = p2, p3 = p3,
                                           gridSize = gridSize)
#View(sorted_srepDataAdjustedG2)

sorted_srepDataAdjustedG3<-fixFilipedSreps(srepsData = srepsDataG3,
                                           reference_u0 = reference_u0,
                                           reference_u1 = reference_u1,
                                           reference_u2 = reference_u2,
                                           centerPoint = centerPoint, p2 = p2, p3 = p3,
                                           gridSize = gridSize)

sorted_srepDataAdjustedG4<-fixFilipedSreps(srepsData = srepsDataG4,
                                           reference_u0 = reference_u0,
                                           reference_u1 = reference_u1,
                                           reference_u2 = reference_u2,
                                           centerPoint = centerPoint, p2 = p2, p3 = p3,
                                           gridSize = gridSize)


#####################################################################################################
#####################################################################################################
# read again skeletal and Boundary PDM of G1 and G2 after relabeling
source("subFunctions/readSrepsData.R")
tempG1<-readSrepsData(srepsData = sorted_srepDataAdjustedG1)
SkeletalPDMG1<-tempG1$SkeletalPDM
BoundaryPDMG1<-tempG1$BoundaryPDM
boundaryPlusSkeletal_G1<-tempG1$boundaryPlusSkeletal

tempG2<-readSrepsData(srepsData = sorted_srepDataAdjustedG2)
SkeletalPDMG2<-tempG2$SkeletalPDM
BoundaryPDMG2<-tempG2$BoundaryPDM
boundaryPlusSkeletal_G2<-tempG2$boundaryPlusSkeletal

tempG3<-readSrepsData(srepsData = sorted_srepDataAdjustedG3)
SkeletalPDMG3<-tempG3$SkeletalPDM
BoundaryPDMG3<-tempG3$BoundaryPDM
boundaryPlusSkeletal_G3<-tempG3$boundaryPlusSkeletal

tempG4<-readSrepsData(srepsData = sorted_srepDataAdjustedG4)
SkeletalPDMG4<-tempG4$SkeletalPDM
BoundaryPDMG4<-tempG4$BoundaryPDM
boundaryPlusSkeletal_G4<-tempG4$boundaryPlusSkeletal

#plot test
for (k in 1:4) {
  skelG1_1<-SkeletalPDMG1[skelRange,,k]
  # open3d()
  referencePoint<-skelG1_1[6,]
  plot3d(skelG1_1,type="p", size=0.5,col = "red",expand = 10,box=FALSE,add = TRUE)
  for (i in 1:(dim(skelG1_1)[1]-1)) {
    plot3d(skelG1_1[c(i,i+1),],type="l",lwd = 0.5,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  plot3d(rbind(referencePoint,referencePoint),type="s",size = 0.8,col = "blue",expand = 10,box=FALSE,add = TRUE)
}
# Important plot to see that up and down is OK
# You must see up down and skeletal points in 3 clean seperated colors
open3d()
for (k in 1:nSamplesG1) {
  plot3d(SkeletalPDMG1[1:upSpoeksNumber,,k],type="p",col = "yellow",expand = 10,box=FALSE,add = TRUE)
  plot3d(BoundaryPDMG1[1:upSpoeksNumber,,k],type="p",col = "orange",expand = 10,box=FALSE,add = TRUE)
  plot3d(BoundaryPDMG1[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),,k],type="p",col = "green",expand = 10,box=FALSE,add = TRUE)
  plot3d(SkeletalPDMG1[(2*upSpoeksNumber+1):nTotalRadii,,k],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
}

sampleNo<-1
srep1<-rbind(SkeletalPDMG1[,,sampleNo],BoundaryPDMG1[,,sampleNo])
open3d()
plot3d(SkeletalPDMG1[skelRange,,sampleNo],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
for (i in 1:nTotalRadii) {
  plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.5,col = "red",expand = 10,box=FALSE,add = TRUE)
}

#####################################################################################################
#####################################################################################################
# alignment

#Hippo
all_G1_G2<-abind(boundaryPlusSkeletal_G1,boundaryPlusSkeletal_G2)
procAllG1_G2<-procGPA(all_G1_G2, scale = F)   # scale = T for shape analysis
alignedAllG1_G2<-procAllG1_G2$rotated  #pooled group
alignedSkeletalPlusBoundaryG1<-alignedAllG1_G2[,,1:nSamplesG1]  #G1 group
alignedSkeletalPlusBoundaryG2<-alignedAllG1_G2[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)] #G2 group
skeletalAlignedG1<-alignedSkeletalPlusBoundaryG1[1:nTotalRadii,,]
boundaryAlignedG1<-alignedSkeletalPlusBoundaryG1[(nTotalRadii+1):(2*nTotalRadii),,]
skeletalAlignedG2<-alignedSkeletalPlusBoundaryG2[1:nTotalRadii,,]
boundaryAlignedG2<-alignedSkeletalPlusBoundaryG2[(nTotalRadii+1):(2*nTotalRadii),,]

#Caud
all_G3_G4<-abind(boundaryPlusSkeletal_G3,boundaryPlusSkeletal_G4)
procAllG3_G4<-procGPA(all_G3_G4, scale = F)   # scale = T for shape analysis
alignedAllG3_G4<-procAllG3_G4$rotated  #pooled group
alignedSkeletalPlusBoundaryG3<-alignedAllG3_G4[,,1:nSamplesG3]  #G1 group
alignedSkeletalPlusBoundaryG4<-alignedAllG3_G4[,,(nSamplesG3+1):(nSamplesG3+nSamplesG4)] #G2 group
skeletalAlignedG3<-alignedSkeletalPlusBoundaryG3[1:nTotalRadii,,]
boundaryAlignedG3<-alignedSkeletalPlusBoundaryG3[(nTotalRadii+1):(2*nTotalRadii),,]
skeletalAlignedG4<-alignedSkeletalPlusBoundaryG4[1:nTotalRadii,,]
boundaryAlignedG4<-alignedSkeletalPlusBoundaryG4[(nTotalRadii+1):(2*nTotalRadii),,]

open3d()
for (k in 1:nSamplesG1) {
  plot3d(skeletalAlignedG1[1:upSpoeksNumber,,k],type="p",col = "yellow",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryAlignedG1[1:upSpoeksNumber,,k],type="p",col = "orange",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryAlignedG1[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),,k],type="p",col = "green",expand = 10,box=FALSE,add = TRUE)
  plot3d(skeletalAlignedG1[(2*upSpoeksNumber+1):nTotalRadii,,k],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
}
for (k in 1:nSamplesG4) {
  plot3d(skeletalAlignedG4[1:upSpoeksNumber,,k],type="p",col = "yellow",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryAlignedG4[1:upSpoeksNumber,,k],type="p",col = "orange",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryAlignedG4[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),,k],type="p",col = "green",expand = 10,box=FALSE,add = TRUE)
  plot3d(skeletalAlignedG4[(2*upSpoeksNumber+1):nTotalRadii,,k],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
}
#####################################################################################################
#####################################################################################################
#find radii after scale and align
radiiG1<-array(NA, dim=c(nTotalRadii,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:nTotalRadii) {
    radiiG1[i,k]<-norm(alignedSkeletalPlusBoundaryG1[i,,k]-
                         alignedSkeletalPlusBoundaryG1[(i+nTotalRadii),,k],type = "2")
  }
}
radiiG2<-array(NA, dim=c(nTotalRadii,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:nTotalRadii) {
    radiiG2[i,k]<-norm(alignedSkeletalPlusBoundaryG2[i,,k]-
                         alignedSkeletalPlusBoundaryG2[(i+nTotalRadii),,k],type = "2")
  }
}
radiiG3<-array(NA, dim=c(nTotalRadii,nSamplesG3))
for (k in 1:nSamplesG3) {
  for (i in 1:nTotalRadii) {
    radiiG3[i,k]<-norm(alignedSkeletalPlusBoundaryG3[i,,k]-
                         alignedSkeletalPlusBoundaryG3[(i+nTotalRadii),,k],type = "2")
  }
}
radiiG4<-array(NA, dim=c(nTotalRadii,nSamplesG4))
for (k in 1:nSamplesG4) {
  for (i in 1:nTotalRadii) {
    radiiG4[i,k]<-norm(alignedSkeletalPlusBoundaryG4[i,,k]-
                         alignedSkeletalPlusBoundaryG4[(i+nTotalRadii),,k],type = "2")
  }
}

#####################################################################################################
#####################################################################################################
# find directions
vectorsG1<-boundaryAlignedG1-skeletalAlignedG1
DirectionsG1<-array(NA,dim = dim(vectorsG1))
for (i in 1:nSamplesG1) {
  DirectionsG1[,,i]<-vectorsG1[,,i]/sqrt(rowSums(vectorsG1[,,i]^2))
}
vectorsG2<-boundaryAlignedG2-skeletalAlignedG2
DirectionsG2<-array(NA,dim = dim(vectorsG2))
for (i in 1:nSamplesG2) {
  DirectionsG2[,,i]<-vectorsG2[,,i]/sqrt(rowSums(vectorsG2[,,i]^2))
}
vectorsG3<-boundaryAlignedG3-skeletalAlignedG3
DirectionsG3<-array(NA,dim = dim(vectorsG3))
for (i in 1:nSamplesG3) {
  DirectionsG3[,,i]<-vectorsG3[,,i]/sqrt(rowSums(vectorsG3[,,i]^2))
}
vectorsG4<-boundaryAlignedG4-skeletalAlignedG4
DirectionsG4<-array(NA,dim = dim(vectorsG4))
for (i in 1:nSamplesG4) {
  DirectionsG4[,,i]<-vectorsG4[,,i]/sqrt(rowSums(vectorsG4[,,i]^2))
}

open3d()
spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
plot3d(t(DirectionsG3[1,,]),type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
plot3d(t(DirectionsG4[1,,]),type="p",col = "red",expand = 10,box=FALSE,add = TRUE)

#####################################################################################################
#####################################################################################################
# Check a sample just to make sure the correctness of the obtained radii and directions after alignment
# we can skip this part 

sampleNo<-6 #change it to check other samples

sumDis<-0
for (i in 1:nTotalRadii) {
  temp1<-skeletalAlignedG1[i,,sampleNo]+DirectionsG1[i,,sampleNo]*radiiG1[i,sampleNo]
  temp2<-boundaryAlignedG1[i,,sampleNo]
  sumDis<-sumDis+norm(temp1-temp2,type = "2")
}
if(sumDis<1e-5){print("Perfect match! Radii and directions are acceptable.")}else{"Imperfect match!!! Do not proceed."}

# plot 
# here we should see only one color!!!
plot3d(skeletalAlignedG1[,,sampleNo]+DirectionsG1[,,sampleNo]*radiiG1[,sampleNo],
       type="p",col = "red" ,expand = 10,box=TRUE,add = TRUE)
plot3d(boundaryAlignedG1[,,sampleNo],type="p",col = "blue" ,expand = 10,box=TRUE,add = TRUE)

#####################################################################################################
#####################################################################################################
#################################### CPNS - CPNG  ###################################################
#####################################################################################################
#####################################################################################################
# CPNS and CPNG

# library(DWD)
# library(caret)
# library(e1071)
# library("kdensity")

# plot function
plotEignmodes <- function(matrixZ,mainTitle) {
  #plot eigenmodes contribution
  s<-prcomp(t(matrixZ))
  pScores<-s$x        #The coordinates of the observations on the principal components
  eigenVectors<-s$rotation
  eigenValues<-(s$sdev)^2
  #Plot the graph for PCA proportions and save as an image
  eigen_prop_G1<- eigenValues/sum(eigenValues)*100
  xx=barplot(eigen_prop_G1[1:20],xlab = 'No. of Eigenmodes=20',
             ylab = 'Percentage contribution', col = "blue",
             main= mainTitle,border = par("fg"),ylim =c(0,110))
  y<-round(eigen_prop_G1[1:20], digits = 2)
  text(x=xx,y=y,labels=as.character(y),pos = 3 , col = "blue", cex = 0.9)
  box()
  cum_prop_G1 <- cumsum(eigenValues)/sum(eigenValues)*100
  lines(cum_prop_G1)
}

# First calculate pooled Zcomp (for classification data must be in a common feature space)
skelG1<-skeletalAlignedG1[skelRange,,]
skelG2<-skeletalAlignedG2[skelRange,,]
skelG3<-skeletalAlignedG3[skelRange,,]
skelG4<-skeletalAlignedG4[skelRange,,]

plot3d(skelG3[,,1],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)

for (i in 1:nSamplesG1) {
  plot3d(skelG1[,,i],type="p", size=0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
}
#####################################################################################################
#####################################################################################################
# Zcomp Hippocampus

# G1 PDMs into the unit sphere and save scale factors
scaleOfSkeletal_G1<-c()
preshapeMatrix_G1<-c()
for (i in 1:nSamplesG1) {
  centeredTemp<-skelG1[,,i] - rep(colMeans(skelG1[,,i]), each = nrow(skelG1[,,i]))
  scaleOfSkeletal_G1<-c(scaleOfSkeletal_G1,centroid.size(centeredTemp))
  preshapeMatrix_G1<-cbind(preshapeMatrix_G1,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G1^2) #test that all shapes are on hypersphere
# G2 PDMs into the unit sphere and save scale factors
scaleOfSkeletal_G2<-c()
preshapeMatrix_G2<-c()
for (i in 1:nSamplesG2) {
  centeredTemp<-skelG2[,,i] - rep(colMeans(skelG2[,,i]), each = nrow(skelG2[,,i]))
  scaleOfSkeletal_G2<-c(scaleOfSkeletal_G2,centroid.size(centeredTemp))
  preshapeMatrix_G2<-cbind(preshapeMatrix_G2,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G2^2) #test that all shapes are on hypersphere

# PNG for pooled skeletal positions
preshapeMatrix_Hippo_pooled<-cbind(preshapeMatrix_G1,preshapeMatrix_G2)

source('mypns.R')
source('mygetSubSphere.R')
source('mysphereFit.R')
# comment out to recalculate PNS for pooled PDMs. Takes a lot of time!!!
####
pns4pdm_Hippo_pooled<-mypns(x = preshapeMatrix_Hippo_pooled, sphere.type = "small")
####
# save(pns4pdm_Hippo_pooled,file = "../files/pns4pdm_Hippo_pooled.Rdata")
####
# load calculated PNG of pooled PDMs. 
# load("../files/pns4pdm_Hippo_pooled.Rdata")

# plot eigenmodes of ZShape_pooled
ZShape_pooled_Hippo<-pns4pdm_Hippo_pooled$resmat
plotEignmodes(ZShape_pooled_Hippo,mainTitle = "Pooled PDMs")

# pooled scale
pooledScale_Hippo<-c(scaleOfSkeletal_G1,scaleOfSkeletal_G2)

#Size of pooled PDMs
sizePDMPooled_Hippo <-pooledScale_Hippo  #gamma(i)
meanSizePDMPooled_Hippo <-exp(mean(log(sizePDMPooled_Hippo)))  #gamma_bar = GM of gamma(i)
normalizedSizePDMPooled_Hippo <- log( sizePDMPooled_Hippo / meanSizePDMPooled_Hippo ) #gamma(i)* = log(gamma(i)/gamma_bar): normalized size variables

#Radii for pooled CPNS
radii_pooled_Hippo<-cbind(radiiG1,radiiG2)
logR_pooled_Hippo <-log(radii_pooled_Hippo)
meanLogR_pooled_Hippo <-rowMeans(logR_pooled_Hippo)
meanRadii_pooled_Hippo <- exp( meanLogR_pooled_Hippo )  #r_bar_i: mean of i-th radius along all samples

#r_bar(i, :): mean of i-th radius for all samples
#used as a scaling factor for radii and spoke directions
rScaleFactors_pooled_Hippo<-array(meanRadii_pooled_Hippo,dim=c(nTotalRadii,(nSamplesG1+nSamplesG2)))

# directions scale factors
uScaleFactors_pooled_Hippo <-  array(NA, dim=c(2*nTotalRadii, (nSamplesG1+nSamplesG2)))
for (i in c(1:nTotalRadii)){
  uScaleFactors_pooled_Hippo[2*(i-1) + 1,]<-rScaleFactors_pooled_Hippo[i,]
  uScaleFactors_pooled_Hippo[2*(i-1) + 2,]<-rScaleFactors_pooled_Hippo[i,]
}

#shifted log r-i's
RStar_pooled_Hippo <-logR_pooled_Hippo - array( meanLogR_pooled_Hippo, dim=c(nTotalRadii,(nSamplesG1+nSamplesG2)))     #R*_i = logR - logR_bar_i

# PNG for Pooled directions
ZSpoke_pooled_Hippo<-c()
DirectionsPooled_Hippo<-abind(DirectionsG1,DirectionsG2)
for (i in c(1 : nTotalRadii)){
  PNS_Directions<- pns(DirectionsPooled_Hippo[i,,], sphere.type = "small")
  residualsTemp<-PNS_Directions$resmat
  ZSpoke_pooled_Hippo<-rbind(ZSpoke_pooled_Hippo,residualsTemp)
}

# CPNG Pooled space 
# CPNG: Construct composite Z matrix and perform PCA
ZComp_pooled_Hippo <- rbind(meanSizePDMPooled_Hippo * ZShape_pooled_Hippo,
                            meanSizePDMPooled_Hippo * normalizedSizePDMPooled_Hippo,
                            rScaleFactors_pooled_Hippo * RStar_pooled_Hippo,
                            uScaleFactors_pooled_Hippo * ZSpoke_pooled_Hippo)

dim(ZComp_pooled_Hippo) # dimention of the feature space
plotEignmodes(ZComp_pooled_Hippo,mainTitle = "Pooled CPNG")

# # extract ZComp features of G1 and G2
# ZComp_G1<-ZComp_pooled_Hippo[,1:nSamplesG1]
# ZComp_G2<-ZComp_pooled_Hippo[,(nSamplesG1+1):(nSamplesG1+nSamplesG2)]

#####################################################################################################
#####################################################################################################
# Zcomp Caudate

# G3 PDMs into the unit sphere and save scale factors
scaleOfSkeletal_G3<-c()
preshapeMatrix_G3<-c()
for (i in 1:nSamplesG3) {
  centeredTemp<-skelG3[,,i] - rep(colMeans(skelG3[,,i]), each = nrow(skelG3[,,i]))
  scaleOfSkeletal_G3<-c(scaleOfSkeletal_G3,centroid.size(centeredTemp))
  preshapeMatrix_G3<-cbind(preshapeMatrix_G3,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G3^2) #test that all shapes are on hypersphere
# G4 PDMs into the unit sphere and save scale factors
scaleOfSkeletal_G4<-c()
preshapeMatrix_G4<-c()
for (i in 1:nSamplesG4) {
  centeredTemp<-skelG4[,,i] - rep(colMeans(skelG4[,,i]), each = nrow(skelG4[,,i]))
  scaleOfSkeletal_G4<-c(scaleOfSkeletal_G4,centroid.size(centeredTemp))
  preshapeMatrix_G4<-cbind(preshapeMatrix_G4,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G4^2) #test that all shapes are on hypersphere

# PNG for pooled skeletal positions
preshapeMatrix_Caud_pooled<-cbind(preshapeMatrix_G3,preshapeMatrix_G4)

source('mypns.R')
source('mygetSubSphere.R')
source('mysphereFit.R')
# comment out to recalculate PNS for pooled PDMs. Takes a lot of time!!!
####
pns4pdm_Caud_pooled<-mypns(x = preshapeMatrix_Caud_pooled, sphere.type = "small")
####
# save(pns4pdm_Caud_pooled,file = "../files/pns4pdm_Caud_pooled.Rdata")
####
# load calculated PNG of pooled PDMs. 
# load("../files/pns4pdm_Caud_pooled_pooled.Rdata")

# plot eigenmodes of ZShape_pooled
ZShape_pooled_Caud<-pns4pdm_Caud_pooled$resmat
plotEignmodes(ZShape_pooled_Caud,mainTitle = "Pooled PDMs")

# pooled scale
pooledScale_Caud<-c(scaleOfSkeletal_G3,scaleOfSkeletal_G4)

#Size of pooled PDMs
sizePDMPooled_Caud <-pooledScale_Caud  #gamma(i)
meanSizePDMPooled_Caud <-exp(mean(log(sizePDMPooled_Caud)))  #gamma_bar = GM of gamma(i)
normalizedSizePDMPooled_Caud <- log( sizePDMPooled_Caud / meanSizePDMPooled_Caud ) #gamma(i)* = log(gamma(i)/gamma_bar): normalized size variables

#Radii for pooled CPNS
radii_pooled_Caud<-cbind(radiiG3,radiiG4)
logR_pooled_Caud <-log(radii_pooled_Caud)
meanLogR_pooled_Caud <-rowMeans(logR_pooled_Caud)
meanRadii_pooled_Caud <- exp( meanLogR_pooled_Caud )  #r_bar_i: mean of i-th radius along all samples

#r_bar(i, :): mean of i-th radius for all samples
#used as a scaling factor for radii and spoke directions
rScaleFactors_pooled_Caud<-array(meanRadii_pooled_Caud,dim=c(nTotalRadii,(nSamplesG3+nSamplesG4)))

# directions scale factors
uScaleFactors_pooled_Caud <-  array(NA, dim=c(2*nTotalRadii, (nSamplesG3+nSamplesG4)))
for (i in c(1:nTotalRadii)){
  uScaleFactors_pooled_Caud[2*(i-1) + 1,]<-rScaleFactors_pooled_Caud[i,]
  uScaleFactors_pooled_Caud[2*(i-1) + 2,]<-rScaleFactors_pooled_Caud[i,]
}

#shifted log r-i's
RStar_pooled_Caud <-logR_pooled_Caud - array( meanLogR_pooled_Caud, dim=c(nTotalRadii,(nSamplesG3+nSamplesG4)))     #R*_i = logR - logR_bar_i

# PNG for Pooled directions
ZSpoke_pooled_Caud<-c()
DirectionsPooled_Caud<-abind(DirectionsG3,DirectionsG4)
for (i in c(1 : nTotalRadii)){
  PNS_Directions<- pns(DirectionsPooled_Caud[i,,], sphere.type = "small")
  residualsTemp<-PNS_Directions$resmat
  ZSpoke_pooled_Caud<-rbind(ZSpoke_pooled_Caud,residualsTemp)
}

# CPNG Pooled space 
# CPNG: Construct composite Z matrix and perform PCA
ZComp_pooled_Caud <- rbind(meanSizePDMPooled_Caud * ZShape_pooled_Caud,
                           meanSizePDMPooled_Caud * normalizedSizePDMPooled_Caud,
                           rScaleFactors_pooled_Caud * RStar_pooled_Caud,
                           uScaleFactors_pooled_Caud * ZSpoke_pooled_Caud)

dim(ZComp_pooled_Caud) # dimention of the feature space
plotEignmodes(ZComp_pooled_Caud,mainTitle = "Pooled CPNG")

# # extract ZComp features of G3 and G4
# ZComp_G3<-ZComp_pooled_Caud[,1:nSamplesG3]
# ZComp_G4<-ZComp_pooled_Caud[,(nSamplesG3+1):(nSamplesG3+nSamplesG4)]


#####################################################################################################
#####################################################################################################
#normalize 
library(bestNormalize)

bigNum<-10000 #boxcox only accept positive values
ZComp_pooled_Hippo_Normalized<-array(NA,dim = dim(ZComp_pooled_Hippo))
for (i in 1:dim(ZComp_pooled_Hippo)[1]) {
  ZComp_pooled_Hippo_Normalized[i,]<-boxcox(ZComp_pooled_Hippo[i,]+bigNum)$x.t
}
ZComp_pooled_Caud_Normalized<-array(NA,dim = dim(ZComp_pooled_Caud))
for (i in 1:dim(ZComp_pooled_Caud)[1]) {
  ZComp_pooled_Caud_Normalized[i,]<-boxcox(ZComp_pooled_Caud[i,]+bigNum)$x.t
}

#plotEignmodes(ZComp_pooled_Caud_Normalized,mainTitle = "Pooled CPNG")
#####################################################################################################
#####################################################################################################
#classification

# diProperm
library(diproperm)
library(DWDLargeR)
library(r.jive)
# devtools::install_github("idc9/r_jive")
library(ajive)

# # example for DiProPerm
# X <- Matrix::t(mushrooms$X) #NB! 8124 objects with 112 features 
# y <- mushrooms$y            #NB! y labels 8124 objects
# dpp <- DiProPerm(X=X,y=y,B = 10)
# diproperm::plotdpp(dpp)

# example for AJIVE
# sample a toy dataset with true joint rank of 1
# blocks <- sample_toy_data(n=700, dx=170, dy=170)
# data_blocks_heatmap(blocks, show_color_bar=FALSE)
# initial_signal_ranks <- c(2, 3) # set by looking at scree plots
# jive_results <- ajive(blocks, initial_signal_ranks)
# # estimated joint rank
# jive_results$joint_rank
# decomposition_heatmaps(blocks, jive_results)
# jive_results$joint_scores
# X=jive_results$block_decomps[[2]][['joint']][['full']]
# y=c(rep(-1,10),rep(1,30))
# dpp <- DiProPerm(X=X,y=y)


blocks_Hippo_Caud<-list("Hippo"=ZComp_pooled_Hippo_Normalized, "Caud"=ZComp_pooled_Caud_Normalized)
data_blocks_heatmap(blocks_Hippo_Caud, show_color_bar=FALSE)
data_blocks_heatmap(blocks_Hippo_Caud, show_color_bar=FALSE)
# #jive #Doesn't work properly
# result<-jive(blocks_Hippo_Caud)
# result$rankJ
# showHeatmaps(result)
#ajive

initial_signal_ranks <- c(50, 50) # set by looking at scree plots
jive_results_Hippo_Caud <- ajive(blocks_Hippo_Caud, initial_signal_ranks)
jive_results_Hippo_Caud$joint_rank
decomposition_heatmaps(blocks_Hippo_Caud, jive_results_Hippo_Caud)
hippo_Joint<-jive_results_Hippo_Caud$block_decomps[[1]][['joint']][['full']]
dim(hippo_Joint)
# hippo_neg<-hippo_Joint[1:nSamplesG1,]
# hippo_pos<-hippo_Joint[(nSamplesG1+1):(nSamplesG1+nSamplesG2),]
X<-t(hippo_Joint)                          #NB! 117 objects with 694 features (each row one object!)
y<-c(rep(-1,nSamplesG1),rep(1,nSamplesG2))   #NB! y labels 117 objects
dpp1 <- DiProPerm(X=X, y=y, cores = 1, classifier = "md")
dpp1$pvalue
diproperm::plotdpp(dpp1)
dpp2 <- DiProPerm(X=X, y=y, cores = 1, classifier = "dwd")
dpp2$pvalue
diproperm::plotdpp(dpp2)

######################################
# cut! we need to cut zero rows from Zcomp matrix.

blocks_Hippo_Caud<-list("Hippo"=t(ZComp_pooled_Hippo_Normalized), "Caud"=t(ZComp_pooled_Caud_Normalized))
data_blocks_heatmap(blocks_Hippo_Caud, show_color_bar=TRUE)
n1<-dim(ZShape_pooled_Hippo)[1]+1
n2<-dim(ZComp_pooled_Hippo_Normalized)[1]
blocks_Hippo_Caud<-list("Hippo"=t(ZComp_pooled_Hippo_Normalized[c(1:20,n1:n2),]),
                        "Caud"=t(ZComp_pooled_Caud_Normalized[c(1:20,n1:n2),]))
data_blocks_heatmap(blocks_Hippo_Caud, show_color_bar=FALSE)
# #jive
# result<-jive(blocks_Hippo_Caud)
# result$rankJ
# showHeatmaps(result)
#ajive
r_init = 50
initial_signal_ranks <- c(r_init, r_init) # set by looking at scree plots
jive_results_Hippo_Caud <- ajive(blocks_Hippo_Caud, initial_signal_ranks, joint_rank = 45)
jive_results_Hippo_Caud$joint_rank
decomposition_heatmaps(blocks_Hippo_Caud, jive_results_Hippo_Caud)
hippo_Joint<-jive_results_Hippo_Caud$block_decomps[[1]][['joint']][['full']]
hippo_Joint <- jive_results_Hippo_Caud$joint_scores
dim(hippo_Joint)
# hippo_neg<-hippo_Joint[1:nSamplesG1,]
# hippo_pos<-hippo_Joint[(nSamplesG1+1):(nSamplesG1+nSamplesG2),]
#X<-t(hippo_Joint)                          #NB! 117 objects with 694 features (each row one object!)
X<-hippo_Joint
y<-c(rep(-1,nSamplesG1),rep(1,nSamplesG2))   #NB! y labels 117 objects
dpp1 <- DiProPerm(X=X, y=y, cores = 1, classifier = "md")
dpp1$pvalue
diproperm::plotdpp(dpp1)
dpp2 <- DiProPerm(X=X, y=y, cores = 1, classifier = "dwd")
dpp2$pvalue
diproperm::plotdpp(dpp2)
# dev.set()
# dev.list()
# dev.close(which=mydev)
