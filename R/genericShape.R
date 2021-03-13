#' @export

genericShape<-methods::setClass(Class = 'genericShape',
                        slots = c(majorAxis = 'numeric',
                                  minorAxis = 'numeric',
                                  roundness = 'numeric',
                                  nArms = 'numeric',
                                  armExt = 'numeric',
                                  armElbow = 'numeric',
                                  armSwing = 'numeric',
                                  armTrig = 'numeric'),
                       prototype = list( majorAxis = c(10,1),
                                         minorAxis = c(10,1),
                                         roundness = c(0.8,0.1),
                                         nArms = c(5,1),
                                         armExt = c(1,1),
                                         armElbow = 3,
                                         armSwing = 3,
                                         armTrig = c(10,1)))

#' @export
methods::setGeneric(name = 'bh_create',
                    signature = 'x',
                    function(x){})

methods::setMethod(f = 'bh_create',
                   signature = 'genericShape',
                   definition = function(x){
                     browser()
                     MA <- rnorm(n = 1, x@majorAxis)
                     MI <- rnorm(n = 1, x@minorAxis)
                     RO <- rnorm(n = 1, x@roundness)
                     NR <- rpois(n = 1, x@nArms)
                     AE <- rnorm(n = NR, MA, x@armExt[2]) * x@armExt[1]
                     AW <- x@armElbow
                     AS <- x@armSwing
                     AT <- rnorm(n = NR, x@armTrig)
                     
                     angl<-sample(1:360,NR)*pi/180
                     
                     skel<-lapply(1:NR,function(i){
                       elbList<-seq(from = 0,
                                    to = AE[i],
                                    by = AE[i]/AW)
                       elb<-lapply(1:AW,function(ii){
                         newAngl <- angl[i]+sample(seq(from=-(AS*pi/180),
                                                       to = (AS*pi/180),
                                                       length.out = 10),1)
                       x<-cos(newAngl)*elbList[ii]
                       y<-sin(newAngl)*elbList[ii]
                     return(matrix(c(x,y),ncol=2,byrow = T))
                       })
                       out<-do.call(rbind,elb)
                       out<-rbind(matrix(c(0,0),ncol = 2,byrow = T),
                                  out)
                       return(out)
                     })
                     skel<-sf::st_multilinestring(x = skel,dim = 'XYZ')
                   })
