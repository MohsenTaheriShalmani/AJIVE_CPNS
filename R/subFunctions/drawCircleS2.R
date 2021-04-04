# draw circle on unit sphere S2 by center of small circle and r
# converted code of Sungkyu Jung MATLAB drawCircleS2.m
drawCircleS2 <- function(center,theta) {
  # NB!!! theta is the angle from center
  if(theta==pi/2){
    t<-cbind(cos(seq(0,2*pi,length.out = 50)),sin(seq(0,2*pi,length.out = 50)),rep(0,50))
    sCirc<-t%*%rotMat(center,c(0,0,1))
  }else{
    t<-cbind(sin(theta)*cos(seq(0,2*pi,length.out = 50)),sin(theta)*sin(seq(0,2*pi,length.out = 50)),cos(theta)*rep(1,50))
    sCirc<-t%*%rotMat(center,c(0,0,1))
  }
  spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
  plot3d(sCirc,type="l",col = "black",expand = 10,box=TRUE,add = TRUE)
}
