library(shapes)
shape1<-matrix(c(2,2,4,2,2,4,4,4),ncol = 2,byrow = TRUE)
shape2<-matrix(c(5,4,9,4,5,8,9,8),ncol = 2,byrow = TRUE)
plotshapes(abind(shape1,shape2))

preshape1<-as.vector(shape1)/sqrt(sum(as.vector(shape1)^2))
preshape2<-as.vector(shape2)/sqrt(sum(as.vector(shape2)^2))
    
sum(preshape1^2)
sum(preshape2^2)

preshape1_back<-matrix(preshape1,ncol = 2,byrow = FALSE)
preshape2_back<-matrix(preshape2,ncol = 2,byrow = FALSE)

plotshapes(abind(preshape1_back,preshape2_back))


#fix issue
middleOfCentroids<-(colMeans(shape1)+colMeans(shape2))/2
middleOfCentroids

centeredByMiddleOfCentroid_1<-shape1-rep(middleOfCentroids, each = nrow(shape1))
centeredByMiddleOfCentroid_2<-shape2-rep(middleOfCentroids, each = nrow(shape2))
plotshapes(abind(centeredByMiddleOfCentroid_1,centeredByMiddleOfCentroid_2))

preshapeCentered1<-as.vector(centeredByMiddleOfCentroid_1)/sqrt(sum(as.vector(centeredByMiddleOfCentroid_1)^2))
preshapeCentered2<-as.vector(centeredByMiddleOfCentroid_2)/sqrt(sum(as.vector(centeredByMiddleOfCentroid_2)^2))
sum(preshapeCentered1^2)
sum(preshapeCentered2^2)

preshapeCentered1_back<-matrix(preshapeCentered1,ncol = 2,byrow = FALSE)
preshapeCentered2_back<-matrix(preshapeCentered2,ncol = 2,byrow = FALSE)
plotshapes(abind(preshapeCentered1_back,preshapeCentered2_back))

