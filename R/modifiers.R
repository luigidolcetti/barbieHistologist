#'@export
.modifier.scramble<-function(x){
  x[sample(seq_along(x),length(x))]
}

#'@export
.modifier.addNoise<-function(x){
  x+sample(x,1)
}

#'@export
.matrix.gauss<-function(x,
                        sigma=1){
  raster::focalWeight(x,sigma)
}


#'@export
bh_modifier<-function(raster,
                      wMatrix=matrix(1,ncol=3,nrow=3),
                      fun = .modifier.scramble){
  raster::focal(x = raster,
                w = wMatrix,
                fun = fun)
}
