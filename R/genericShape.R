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
                                                  armTrig = c(1,1)))

#' @export
methods::setGeneric(name = 'bh_create',
                    signature = 'x',
                    function(x){})

methods::setMethod(f = 'bh_create',
                   signature = 'genericShape',
                   definition = function(x){
                     
                     MA <- rnorm(n = 1, x@majorAxis)
                     MI <- rnorm(n = 1, x@minorAxis)
                     RO <- rnorm(n = 1, x@roundness)
                     NR <- rpois(n = 1, x@nArms)
                     AE <- rnorm(n = NR, MA, x@armExt[2]) * x@armExt[1]
                     AW <- x@armElbow
                     AS <- x@armSwing
                     AT <- rnorm(n = NR, x@armTrig)
                     
                     angl<-sample(1:360,NR)*pi/180
                     
                     angl<-sort(angl)
                     
                     skel<-lapply(1:NR,function(i){
                       elbListX<-seq(from = MA,
                                    to = MA+AE[i],
                                    by = AE[i]/AW)
                       elbListY<-seq(from = MI,
                                     to = MI+AE[i],
                                     by = AE[i]/AW)
                       
                       elb<-lapply(1:AW,function(ii){
                         newAngl <- angl[i]+sample(seq(from=-(AS*pi/180),
                                                       to = (AS*pi/180),
                                                       length.out = 10),1)
                         x<-cos(newAngl)*elbListX[ii]
                         y<-sin(newAngl)*elbListY[ii]
                         m1<--(x/y)
                         q1<-y-(m1*x)
                         
                         if (AT[i]<0) DD<-(AW-ii)*(-AT[i]) else DD<-ii*AT[i]
                         x1<-x+sqrt(DD/(1+m1^2))
                         x2<-x+-sqrt(DD/(1+m1^2))
                         y1<-m1*x1+q1
                         y2<-m1*x2+q1
                         out<-list(stem=matrix(c(x,y),ncol=2,byrow = T),
                                   branch=matrix(c(x1,y1,x2,y2),ncol=2,byrow=T))
                         return(out)
                       })
                       
                       out_stem<-do.call(rbind,lapply(elb,function(x)unlist(x$stem,recursive = F)))
                       out_stem<-rbind(matrix(c(0,0),ncol = 2,byrow = T),
                                       out_stem)
                       out_branch<-lapply(elb,function(x)unlist(x$branch,recursive = F))[-1]
                       
                       out<-list(stem=out_stem,
                                 branch=out_branch)
                       return(out)
                     })
                     
                     stm<-sf::st_multilinestring(x = lapply(skel,function(x)unlist(x$stem,recursive = F)),dim = 'XYZ')
                     brnch<-lapply(skel,function(x){
                       sf::st_multilinestring(x = x$branch,dim = 'XYZ')})
                     
                     outLine<-lapply(skel,function(x){
                    
                       lx<-lapply(x$branch,function(y) y[c(T,F),])
                       lx <- do.call(rbind,lx)
                       rx<-lapply(x$branch,function(y) y[c(F,T),])
                       rx<-do.call(rbind,rev(rx))
                       if (sum(lx[,2])>0) {
                       out<-rbind(lx,rx)} else {
                         out<-rbind(rx[nrow(rx):1,],lx[nrow(lx):1,])}
                       # out<-out[order(out[,2],decreasing = T),]
                       })
                     outLine<-do.call(rbind,outLine)
                     outLine<-rbind(outLine,outLine[1,])
                     outLine<-sf::st_polygon(list(outLine))
                     outLine<-sf::st_buffer(outLine,MA/10)
                     outLine<-sf::st_buffer(outLine,-MA/10)
                     
                     # spikes<-round(18*RO)
                     # 
                     # if (spikes<3) spikes<-3
                     # 
                     # angl<-seq(from = 1,
                     #           to = 360,
                     #           length.out = spikes)*pi/180
                     # 
                     # bodyVertex<-lapply(angl,function(a){
                     #   x <- cos(a)*MA
                     #   y <- sin(a)*MI
                     #   matrix(c(x,y),ncol=2,byrow = T)
                     # })
                     # bodyVertex <- do.call(rbind,bodyVertex)
                     # bodyVertex<-rbind(bodyVertex,bodyVertex[1,])
                     # 
                     # rotation<-.rot(sample(1:360,1)*pi/180)
                     # browser()
                     # body<-sf::st_polygon(list(t(apply(bodyVertex,1,'*',rotation))))
                     # 
                     # finalOutLine<-sf::st_union(body,outLine)
                     # 
                     out<-list(stem = stm,
                               branch = brnch,
                               outline = outLine)
                     return(out)
                   })

.rot = function(a) matrix(c(cos(a), sin(a)), ncol = 2, nrow = 1)
                     
                     