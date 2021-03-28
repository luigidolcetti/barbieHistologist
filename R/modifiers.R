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
bh_focal_modifier<-function(raster,
                      wMatrix=matrix(1,3,3),
                      fun = .modifier.scramble,
                      fun.param = NULL){
  
  if (!is.null(fun.param)){
    newFormals<-formals(fun)
    for (i in names(fun.param)) newFormals[[i]]<-fun.param[[i]]
    formals(fun)<-newFormals
    }
  
  raster::focal(x = raster,
                w = wMatrix,
                fun = fun,
                pad=TRUE,
                padValue=0)
}

#' @export
.perturbation.constant<-function(field = c(0,10,0,10),
                                 amount = 0.5,
                                 mean = 10,
                                 sd = 1){

  XY<-expand.grid(seq(field[1],field[2]),seq(field[3],field[4]))
  colnames(XY)<-c('x','y')
  Nrow<-nrow(XY)
  Ncases<-round(Nrow*amount)
  Wcases<-sample(1:Nrow,Ncases)
  Zmodif<-rnorm(Ncases,mean,sd)
  Zmodif[Zmodif<0]<-0
  Z<-rep(0,Nrow)
  Z[Wcases]<-Zmodif
  XYZ<-cbind(XY,Z)
  colnames(XYZ)<-c('x','y','z')
  return(XYZ)
}

#'@export
bh_field_perturbation<-function(raster,
                            fun = .perturbation.constant,
                            fun.param = NULL){
  
  if (!is.null(fun.param)){
    newFormals<-formals(fun)
    for (i in names(fun.param)) newFormals[[i]]<-fun.param[[i]]
    formals(fun)<-newFormals
  }
  
  XYZ<-fun()
  
  tiles<-raster::cellFromXY(raster,XYZ[,1:2])
  raster[tiles]<-raster[tiles]+XYZ[,3]
  return(raster)
}
