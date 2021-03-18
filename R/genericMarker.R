#' @export
.pattern.random<-function(xy=expand.grid(x=1:10,y=1:10),
                          mean = 10,
                          sd = 1){
  z = rnorm(n = nrow(xy),
            mean = mean,
            sd = sd)
  out<-cbind(xy,z)
  return(out)
}

#' @export
.pattern.concentric<-function(xy=expand.grid(x=1:10,y=1:10),
                              mean = 10,
                              sd = 1,
                              ...){
  
  xC<-round(mean(xy[,1]))
  yC<-round(mean(xy[,2]))
  xyD<-sqrt((xy[1]-xC)^2+(xy[2]-yC)^2)
  
  z<-qnorm(p = unname(unlist((xyD-min(xyD))/(max(xyD)-min(xyD)))),
           mean = mean,
           sd = sd)
  
  z[z==Inf]<-max(z[!is.infinite(z)])
  z[z==-Inf]<-min(z[!is.infinite(z)])
  
  out<-cbind(xy,z)
  return(out)
}

#' @export
.pattern.invConcentric<-function(xy=expand.grid(x=1:10,y=1:10),
                                 mean = 10,
                                 sd = 1,
                                 ...){
  
  xC<-round(mean(xy[,1]))
  yC<-round(mean(xy[,2]))
  xyD<-sqrt((xy[1]-xC)^2+(xy[2]-yC)^2)
  
  z<-qnorm(p = unname(unlist((xyD-min(xyD))/(max(xyD)-min(xyD)))),
           mean = mean,
           sd = sd)
  
  z[z==Inf]<-max(z[!is.infinite(z)])
  z[z==-Inf]<-min(z[!is.infinite(z)])
  
  z<-abs(z-max(z))
  
  out<-cbind(xy,z)
  return(out)
}

#' @export
.pattern.skin<-function(xy=expand.grid(x=1:10,y=1:10),
                        mean = 10,
                        sd = 1,
                        trsh = 0.2,
                        ...){
  
  levelXY<-cbind(xy,n=rep(0,nrow(xy)),l=rep(0,nrow(xy)))
  # colnames(levelXY)<-c('x','y','n','l')
  l=1
  shft<-matrix(c(-1,0,0,1,1,0,0,-1),ncol=2,byrow = T)
  while (any(levelXY[,'l']==0)){
    levelXY<-apply(levelXY,1,function(r){
      if (r['n']==0){
        cell<-apply(shft,1,function(sh){
          cc<-which(levelXY[,'x']==r['x']+sh[1] & 
                      levelXY[,'y']==r['y']+sh[2])
          
          if (length(cc)==0) {
            return(1)} else {
              levelXY[cc,'n']
            }
        })
        if (sum(cell)>0) {
          r['n']<-1
          r['l']<-l}
        return(r)
      } else {
        return(r)
      }
    })
    levelXY<-t(levelXY)
    l=l+1
    
  }
  newTrsh<-round((max(levelXY[,'l'])-min(levelXY[,'l']))*trsh)
  z<-rep(0,nrow(xy))
  w<-which(levelXY[,'l']<=newTrsh)
  z[w]<-rnorm(n = length(w),mean = mean,sd = sd)
  
  out<-cbind(xy,z)
  return(out)
}

#' @export
.pattern.skin2<-function(xy=expand.grid(x=1:10,y=1:10),
                        mean = 10,
                        sd = 1,
                        k = 1,
                        X0 = 0.75,
                        ...){
  
  levelXY<-cbind(xy,n=rep(0,nrow(xy)),l=rep(0,nrow(xy)))
  # colnames(levelXY)<-c('x','y','n','l')
  l=1
  shft<-matrix(c(-1,0,0,1,1,0,0,-1),ncol=2,byrow = T)
  while (any(levelXY[,'l']==0)){
    levelXY<-apply(levelXY,1,function(r){
      if (r['n']==0){
        cell<-apply(shft,1,function(sh){
          cc<-which(levelXY[,'x']==r['x']+sh[1] & 
                      levelXY[,'y']==r['y']+sh[2])
          
          if (length(cc)==0) {
            return(1)} else {
              levelXY[cc,'n']
            }
        })
        if (sum(cell)>0) {
          r['n']<-1
          r['l']<-l}
        return(r)
      } else {
        return(r)
      }
    })
    levelXY<-t(levelXY)
    l=l+1
    
  }
  X0<-(max(levelXY[,'l'])-min(levelXY[,'l']))*X0
  
  logisticLevel<-1/(1+exp(k*(levelXY[,'l']-X0)))
  z<-rnorm(n = length(xy),mean = mean,sd = sd)*logisticLevel
  
  out<-cbind(xy,z)
  return(out)
}



#' @export
genericMarker<-methods::setClass(Class = 'genericMarker',
                                 slots = c(marker = 'character',
                                           channel = 'character',
                                           Rname = 'character',
                                           Cname = 'character',
                                           pattern = 'call',
                                           compartment = 'character'),
                                 prototype = list(marker = 'a',
                                                  channel = 'ch1',
                                                  Rname = 'x.a.ch1',
                                                  Cname = 'ch1-a(ch1)',
                                                  pattern = call('c'),
                                                  compartment = 'cytoplasm'))

methods::setClassUnion(name = 'genericMarker_NULL',
                       members = c('genericMarker','NULL'))

#' @export
methods::setMethod('initialize',
                   'genericMarker',
                   function(.Object,
                            marker = 'a',
                            channel = 'ch1',
                            Rname = 'x.ch1.ch1',
                            Cname = 'ch1-ch1(ch1)',
                            pattern = call('c'),
                            compartment = 'cytoplasm',
                            ...){
                     
                     .Object <- methods::callNextMethod(.Object, ...)
                     
                     .Object@marker<-marker
                     .Object@channel<-channel
                     
                     Rname<-paste0('x.',channel,'.',channel)
                     Cname<-paste0(channel,'-',channel,'(',channel,')')
                     
                     .Object@Rname<-Rname
                     .Object@Cname<-Cname
                     
                     if (!any(compartment %in% c('cytoplasm','nucleus','organelle'))) stop ('compartment must be one among cytoplasm, nucleus or organelle')
                     .Object@compartment<-compartment
                     .Object@pattern<-pattern
                     .Object
                   })

#' Define a new marker
#' 
#' helper function to create a new marker.
#' @param markerName character, a name for the marker e.g. CD4.
#' @param channelName character, the channel (layer) in which the marker
#'   will be detected. A corresponding layer should be created in a rasterStack
#'   or tissue.
#' @param patternFunction character, an existing function, either one among
#'   those coming in the package e.g. .pattern.skin() ore a custom function.
#' @param patternModifier named list, a list of parameters (e.g. mean, sd, tresh)
#'   that specify new values. Parameters names must match those of the original
#'   function. xy parameter is excluded because this parameter is passed 
#'   automatically at the moment of taking a picture.
#' @param compartment character, one among cytoplasm, nucleus or organelle
#'   (the last one not yet implemented), that specify the compartment that
#'   express the marker.
#' @return An object of class genericMarker
#' @export
bh_defineMarker<-function(markerName = NULL,
                          channelName = NULL,
                          patternFunction = NULL,
                          patternModifier = NULL,
                          compartment = NULL){
  
  if (is.null(markerName)) stop('give a maker name')
  if (is.null(channelName)) stop('give a channel name')
  if (is.null(patternFunction)) stop('give a pattern function')
  if (!exists(patternFunction)) stop('cannot find this pattern function')
  if (!is.null(patternModifier)){
    if (!is.list(patternModifier)) stop('pattern modifiers must be wrapped in a named list')
    frm<-names(formals(patternFunction))
    if (length(frm)==0) stop('name the list')
    if (!all(names(patternModifier) %in% frm)) stop('check pattern modifier arguments')
  }
  if (is.null(compartment)) stop('give a compartment')
  
  newCall<-rlang::call2(patternFunction,!!!patternModifier)
  
  new('genericMarker',
      marker = markerName,
      channel = channelName,
      pattern = newCall,
      compartment = compartment)
}