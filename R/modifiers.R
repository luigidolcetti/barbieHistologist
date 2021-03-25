#'@export
.modifier.scramble<-function(x){
  x[sample(seq_along(x),1)]
}

#'@export
.modifier.addNoise<-function(x){
  v <- x[x!=0]
  if (length(v)) out<-sample(v,1) else out<-0
  return(out)
}

#'@export
.modifier.multDiv<-function(x,quantity=2){
  whatOp<-(sample(c('*','/'),1))
  out<-eval(parse(text = paste0('x',whatOp,'quantity')))
  mean(out)
}

#'@export
.matrix.gauss<-function(x,
                        sigma=1){
  raster::focalWeight(x,sigma,'Gauss')
}

#'@export
.matrix.circle<-function(x,
                        radius=3){
  raster::focalWeight(x,radius,'circle')
}

#'@export
bh_modifier<-function(raster,
                      wMatrix=matrix(1,3,3),
                      fun = .modifier.scramble,
                      fun.param = NULL){
  
  if (!is.null(fun.param)){
    formals(fun)<-append(alist(x=),fun.param)
  }
  
  raster::focal(x = raster,
                w = wMatrix,
                fun = fun,
                pad=TRUE,
                padValue=0)
}
