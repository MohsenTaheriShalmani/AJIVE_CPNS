library(mvtnorm)
library(shapes)
library(e1071)
library(caret)
library(ajive)

################################################
# simulation 
# individual has a better performance !!!
rm(list = ls())

group1<-gorm.dat
changeMatrix<-t(matrix(c(50,20,0,0,0,0,0,0,0,0,10,0,0,0,0,0),nrow = 2))
group2<-array(NA, dim = dim(group1))
for (i in 1:dim(group1)[3]) {
  group2[,,i]<-group1[,,i]+changeMatrix
}
plotshapes(group1,color = "blue")
plotshapes(group2,color = "red")

preshapeMatrix_G1<-c()
centered_G1<-array(NA,dim = dim(group1))
for (i in 1:dim(group1)[3]) {
  centeredTemp<-group1[,,i] - rep(colMeans(group1[,,i]), each = nrow(group1[,,i]))
  centered_G1[,,i]<-centeredTemp
  preshapeMatrix_G1<-cbind(preshapeMatrix_G1,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G1^2) #test that all shapes are on hypersphere
preshapeMatrix_G2<-c()
centered_G2<-array(NA,dim = dim(group2))
for (i in 1:dim(group2)[3]) {
  centeredTemp<-group2[,,i] - rep(colMeans(group2[,,i]), each = nrow(group2[,,i]))
  centered_G2[,,i]<-centeredTemp
  preshapeMatrix_G2<-cbind(preshapeMatrix_G2,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G2^2) #test that all shapes are on hypersphere

#simulated data
procTest<-procGPA(group1)
# plotshapes(procTest$rotated,color = "blue")
# par(new=T)
# plotshapes(centered_G1,color = "red")
group4simulation<-procTest$rotated
preshapeMatrix_G3<-c()
centered_G3<-array(NA,dim = dim(group1))
for (i in 1:dim(group1)[3]) {
  centeredTemp<-group4simulation[,,i] - rep(colMeans(group4simulation[,,i]), each = nrow(group4simulation[,,i]))
  centered_G3[,,i]<-centeredTemp
  preshapeMatrix_G3<-cbind(preshapeMatrix_G3,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G1^2) #test that all shapes are on hypersphere

preshapeMatrix_G4<-c()
centered_G4<-array(NA,dim = dim(group2))
for (i in 1:dim(group2)[3]) {
  # centeredTemp<-group2[,,i] - rep(colMeans(group2[,,i]), each = nrow(group2[,,i]))
  # centeredTemp<-prcomp(group2[,,i])$x #PCA of group2
  centeredTemp<-cbind(prcomp(group2[,,i])$x[,1],-prcomp(group2[,,i])$x[,2]) #PCA of group2
  centered_G4[,,i]<-centeredTemp
  preshapeMatrix_G4<-cbind(preshapeMatrix_G4,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G4^2) #test that all shapes are on hypersphere

# #block 1
# plotshapes(centered_G1,color = "blue")
# par(new=T)
# plotshapes(centered_G3,color = "green")
# #block 2
# plotshapes(centered_G2,color = "red")
# par(new=T)
# plotshapes(centered_G4,color = "pink")

xlim<-ylim<-c(-180,180)
par(new=F)
for (i in 1:dim(centered_G1)[3]) {
  plot(centered_G1[,,i], pch=1, col="blue" ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y")
  par(new=T)
}
for (i in 1:dim(centered_G3)[3]) {
  plot(centered_G3[,,i], pch=4, col="red" ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y")
  par(new=T)
}
legend("bottomright", legend=c("Group A", "Group B"),pch = c(1,4),
       col=c("blue","red"), cex=0.5,pt.cex = 1, horiz=T)

xlim<-ylim<-c(-180,180)
par(new=F)
for (i in 1:dim(centered_G2)[3]) {
  plot(centered_G2[,,i], pch=1, col="blue" ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y")
  par(new=T)
}
for (i in 1:dim(centered_G4)[3]) {
  plot(centered_G4[,,i], pch=4, col="red" ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y")
  par(new=T)
}
legend("bottomright", legend=c("Group A", "Group B"),pch = c(1,4),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)



dim(preshapeMatrix_G1)

g1<-t(preshapeMatrix_G1)
g1_2<-t(preshapeMatrix_G3)
n_g1<-dim(g1)[1]
g2<-t(preshapeMatrix_G2)
g2_2<-t(preshapeMatrix_G4)
n_g2<-dim(g2)[1]

block1<-rbind(g1,g1_2)
dim(block1)
block1_dataFrame <- data.frame(x=rbind(g1,g1_2))
SVM_block1 <- train(block1_dataFrame, as.factor(c(rep(1,n_g1),rep(0,n_g1))),
                      "svmLinear",
                      tuneLength = 10,
                      trControl = trainControl(method = "cv"))
SVM_block1$results

block2<-rbind(g2,g2_2)
dim(block2)
block2_dataFrame <- data.frame(x=rbind(g2,g2_2))
SVM_block2 <- train(block2_dataFrame, as.factor(c(rep(1,n_g2),rep(0,n_g2))),
                    "svmLinear",
                    tuneLength = 10,
                    trControl = trainControl(method = "cv"))
SVM_block2$results

# example for AJIVE
# sample a toy dataset with true joint rank of 1
# blocks <- sample_toy_data(n=200, dx=100, dy=150)
# dim(blocks$`1`) #NB!!! in block 1 & block 2 number of rows (features) must be equal not columns!!!
# dim(blocks$`2`)
# data_blocks_heatmap(blocks, show_color_bar=FALSE)

#NB!!! in block 1 & block 2 number of rows (features) must be equal not columns!!!
blocks<-list("B1"=block1, "B2"=block2)
dim(blocks$B1)
dim(blocks$B2)
data_blocks_heatmap(blocks, show_color_bar=FALSE)
initial_signal_ranks <- c(2,3)
ajive_results <- ajive(blocks,initial_signal_ranks)
ajive_results$joint_rank
decomposition_heatmaps(blocks, ajive_results)
hippo_Joint<-ajive_results$block_decomps[[1]][['joint']][['full']]
dim(hippo_Joint)
caud_Joint<-ajive_results$block_decomps[[1]][['joint']][['full']]
dim(caud_Joint)

joint_dataFrame_hippo <- data.frame(x=hippo_Joint)
dim(joint_dataFrame_hippo)
SVM_joint_hippo <- train(joint_dataFrame_hippo, as.factor(c(rep(1,n_g1),rep(0,n_g1))),
                         "svmLinear",
                         tuneLength = 10,
                         trControl = trainControl(method = "cv"))
SVM_joint_hippo$results

joint_dataFrame_caud <- data.frame(x=caud_Joint)
dim(joint_dataFrame_caud)
SVM_joint_caud <- train(joint_dataFrame_caud, as.factor(c(rep(1,n_g2),rep(0,n_g2))),
                        "svmLinear",
                        tuneLength = 10,
                        trControl = trainControl(method = "cv"))
SVM_joint_caud$results

#hippo+caud
joint_dataFrame <- data.frame(x=rbind(hippo_Joint,caud_Joint))
dim(joint_dataFrame)
SVM_joint <- train(joint_dataFrame, as.factor(c(rep(1,n_g1),rep(0,n_g1),rep(1,n_g2),rep(0,n_g2))),
                   "svmLinear",
                   tuneLength = 10,
                   trControl = trainControl(method = "cv"))
SVM_joint$results

################################################
# another simulation
# joint has a better performance !!!
rm(list = ls())

group1<-gorm.dat
changeMatrix<-t(matrix(c(50,20,0,0,0,0,0,0,0,0,10,0,0,0,0,0),nrow = 2))
group2<-array(NA, dim = dim(group1))
for (i in 1:dim(group1)[3]) {
  theta<-pi/4
  rotationMatrixTemp<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow = 2,byrow = T)
  group2[,,i]<-(group1[,,i]+changeMatrix)%*%rotationMatrixTemp
}
# plotshapes(group1,color = "blue")
# plotshapes(group2,color = "red")

preshapeMatrix_G1<-c()
centered_G1<-array(NA,dim = dim(group1))
for (i in 1:dim(group1)[3]) {
  centeredTemp<-group1[,,i] - rep(colMeans(group1[,,i]), each = nrow(group1[,,i]))
  centered_G1[,,i]<-centeredTemp
  preshapeMatrix_G1<-cbind(preshapeMatrix_G1,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G1^2) #test that all shapes are on hypersphere
preshapeMatrix_G2<-c()
centered_G2<-array(NA,dim = dim(group2))
for (i in 1:dim(group2)[3]) {
  centeredTemp<-group2[,,i] - rep(colMeans(group2[,,i]), each = nrow(group2[,,i]))
  centered_G2[,,i]<-centeredTemp
  preshapeMatrix_G2<-cbind(preshapeMatrix_G2,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G2^2) #test that all shapes are on hypersphere

#simulated data
procTest<-procGPA(group1)
# plotshapes(procTest$rotated,color = "blue")
# par(new=T)
# plotshapes(centered_G1,color = "red")
group4simulation<-procTest$rotated
preshapeMatrix_G3<-c()
centered_G3<-array(NA,dim = dim(group1))
for (i in 1:dim(group1)[3]) {
  centeredTemp<-group4simulation[,,i] - rep(colMeans(group4simulation[,,i]), each = nrow(group4simulation[,,i]))
  centered_G3[,,i]<-centeredTemp
  preshapeMatrix_G3<-cbind(preshapeMatrix_G3,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G1^2) #test that all shapes are on hypersphere

procTest2<-procGPA(group2)
group4simulation2<-procTest2$rotated
preshapeMatrix_G4<-c()
centered_G4<-array(NA,dim = dim(group2))
for (i in 1:dim(group2)[3]) {
  centeredTemp<-group4simulation2[,,i] - rep(colMeans(group4simulation2[,,i]), each = nrow(group4simulation2[,,i]))
  centered_G4[,,i]<-centeredTemp
  preshapeMatrix_G4<-cbind(preshapeMatrix_G4,as.vector(t(centeredTemp))/centroid.size(centeredTemp))
}
colSums(preshapeMatrix_G4^2) #test that all shapes are on hypersphere


# plotshapes(centered_G1,color = "blue")
# par(new=T)
# plotshapes(centered_G3,color = "green")
# plotshapes(centered_G2,color = "red")
# par(new=T)
# plotshapes(centered_G4,color = "pink")

xlim<-ylim<-c(-180,180)
par(new=F)
for (i in 1:dim(centered_G1)[3]) {
  plot(centered_G1[,,i], pch=1, col="blue" ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y")
  par(new=T)
}
for (i in 1:dim(centered_G3)[3]) {
  plot(centered_G3[,,i], pch=4, col="red" ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y")
  par(new=T)
}
legend("bottomright", legend=c("Group A", "Group B"),pch = c(1,4),
       col=c("blue","red"), cex=0.5,pt.cex = 1, horiz=T)

xlim<-ylim<-c(-180,180)
par(new=F)
for (i in 1:dim(centered_G2)[3]) {
  plot(centered_G2[,,i], pch=1, col="blue" ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y")
  par(new=T)
}
for (i in 1:dim(centered_G4)[3]) {
  plot(centered_G4[,,i], pch=4, col="red" ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y")
  par(new=T)
}
legend("bottomright", legend=c("Group A", "Group B"),pch = c(1,4),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)


dim(preshapeMatrix_G1)

g1<-t(preshapeMatrix_G1)
g1_2<-t(preshapeMatrix_G3)
n_g1<-dim(g1)[1]
g2<-t(preshapeMatrix_G2)
g2_2<-t(preshapeMatrix_G4)
n_g2<-dim(g2)[1]

block1<-rbind(g1,g1_2)
dim(block1)
block1_dataFrame <- data.frame(x=rbind(g1,g1_2))
SVM_block1 <- train(block1_dataFrame, as.factor(c(rep(1,n_g1),rep(0,n_g1))),
                    "svmLinear",
                    tuneLength = 10,
                    trControl = trainControl(method = "cv"))
SVM_block1$results

block2<-rbind(g2,g2_2)
dim(block2)
block2_dataFrame <- data.frame(x=rbind(g2,g2_2))
SVM_block2 <- train(block2_dataFrame, as.factor(c(rep(1,n_g2),rep(0,n_g2))),
                    "svmLinear",
                    tuneLength = 10,
                    trControl = trainControl(method = "cv"))
SVM_block2$results


# example for AJIVE
# sample a toy dataset with true joint rank of 1
# blocks <- sample_toy_data(n=200, dx=100, dy=150)
# dim(blocks$`1`) #NB!!! in block 1 & block 2 number of rows (features) must be equal not columns!!!
# dim(blocks$`2`)
# data_blocks_heatmap(blocks, show_color_bar=FALSE)

#NB!!! in block 1 & block 2 number of rows (features) must be equal not columns!!!
blocks<-list("B1"=block1, "B2"=block2)
dim(blocks$B1)
dim(blocks$B2)
data_blocks_heatmap(blocks, show_color_bar=FALSE)
initial_signal_ranks <- c(2,3)
ajive_results <- ajive(blocks,initial_signal_ranks)
ajive_results$joint_rank
decomposition_heatmaps(blocks, ajive_results)
hippo_Joint<-ajive_results$block_decomps[[1]][['joint']][['full']]
dim(hippo_Joint)
caud_Joint<-ajive_results$block_decomps[[2]][['joint']][['full']]
dim(caud_Joint)

joint_dataFrame_hippo <- data.frame(x=hippo_Joint)
dim(joint_dataFrame_hippo)
SVM_joint_hippo <- train(joint_dataFrame_hippo, as.factor(c(rep(1,n_g1),rep(0,n_g1))),
                         "svmLinear",
                         tuneLength = 10,
                         trControl = trainControl(method = "cv"))
SVM_joint_hippo$results

joint_dataFrame_caud <- data.frame(x=caud_Joint)
dim(joint_dataFrame_caud)
SVM_joint_caud <- train(joint_dataFrame_caud, as.factor(c(rep(1,n_g2),rep(0,n_g2))),
                        "svmLinear",
                        tuneLength = 10,
                        trControl = trainControl(method = "cv"))
SVM_joint_caud$results

#hippo+caud
joint_dataFrame <- data.frame(x=rbind(hippo_Joint,caud_Joint))
dim(joint_dataFrame)
SVM_joint <- train(joint_dataFrame, as.factor(c(rep(1,n_g1),rep(0,n_g1),rep(1,n_g2),rep(0,n_g2))),
                   "svmLinear",
                   tuneLength = 10,
                   trControl = trainControl(method = "cv"))
SVM_joint$results
