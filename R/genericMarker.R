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
                            Rname = 'x.a.ch1',
                            Cname = 'ch1-a(ch1)',
                            pattern = call('c'),
                            compartment = 'cytoplasm',
                            ...){
                     
                     .Object <- methods::callNextMethod(.Object, ...)
                     
                     Rname<-paste0('x.',marker,'.',channel)
                     Cname<-paste0(channel,'-',marker,'(',channel,')')
                     
                     .Object@Rname<-Rname
                     .Object@Cname<-Cname
                     
                     
                     if (!any(compartment %in% c('cytoplasm','nucleus','organelle'))) stop ('compartment must be one among cytoplasm, nucleus or organelle')
                     .Object@compartment<-compartment
                     
                     
                     .Object@pattern<-pattern
                     
                     .Object
                   })
