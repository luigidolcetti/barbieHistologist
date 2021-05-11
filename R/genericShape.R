#' @export
genericShape<-methods::setClass(Class = 'genericShape',
                                slots = c(majorAxis = 'numeric',
                                          minorAxis = 'numeric',
                                          roundness = 'numeric',
                                          nArms = 'numeric',
                                          fixedArms = 'logical',
                                          orientation = 'numeric',
                                          armExt = 'numeric',
                                          armElbow = 'numeric',
                                          armSwing = 'numeric',
                                          armTrig = 'numeric'),
                                prototype = list( majorAxis = c(10,1),
                                                  minorAxis = c(10,1),
                                                  roundness = c(0.8,0.1),
                                                  nArms = 5,
                                                  fixedArms = F,
                                                  orientation = NULL,
                                                  armExt = c(1,1),
                                                  armElbow = 3,
                                                  armSwing = 3,
                                                  armTrig = c(1,1)))

#'@export
genericShape_NULL<-methods::setClassUnion(name = 'genericShape_NULL',
                                          members = c('genericShape','NULL'))

#' @export
methods::setMethod('initialize',
                   'genericShape',
                   function(.Object,
                            majorAxis = c(10,1),
                            minorAxis = c(10,1),
                            roundness = c(0.8,0.1),
                            nArms = 5,
                            fixedArms = F,
                            orientation = NULL,
                            armExt = c(1,1),
                            armElbow = 3,
                            armSwing = 3,
                            armTrig = c(1,1),
                            ...){
                     .Object <- methods::callNextMethod(.Object, ...)
                     if (!is.null(majorAxis)){
                       if (length(majorAxis)!=2) stop('length of one argument is wrong')
                       .Object@majorAxis <- majorAxis}
                     if (!is.null(minorAxis)){
                       if (length(minorAxis)!=2) stop('length of one argument is wrong')
                       .Object@minorAxis <- minorAxis}
                     if (!is.null(roundness)){
                       if (length(roundness)!=2) stop('length of one argument is wrong')
                       if (roundness[1]<0 | roundness[1]>1) stop('roundness must be <1 and >0')
                       .Object@roundness <- roundness}
                     if (!is.null(nArms)){
                       if (length(nArms)!=1) stop('length of one argument is wrong')
                       .Object@nArms <- nArms}
                     if (!is.null(fixedArms)){
                       if (!is.logical(fixedArms)) stop('fixed arms must be T/F')
                       .Object@fixedArms<-fixedArms}
                     if (!is.null(orientation)){
                       if (!is.numeric(orientation)) stop('orientation must be either NULL or an angle expressed in degrees')
                       .Object@orientation<-orientation}
                     if (!is.null(armExt)){
                       if (length(armExt)!=2) stop('length of one argument is wrong')
                       .Object@armExt <- armExt}
                     if (!is.null(armElbow)){
                       if (length(armElbow)!=1) stop('length of one argument is wrong')
                       if (armElbow<2) stop('armElbow must be >=2')
                       .Object@armElbow <- armElbow}
                     if (!is.null(armSwing)){
                       if (length(armSwing)!=1) stop('length of one argument is wrong')
                       .Object@armSwing <- armSwing}
                     if (!is.null(armTrig)){
                       if (length(armTrig)!=2) stop('length of one argument is wrong')
                       .Object@armTrig <- armTrig}
                     .Object
                   })

#' @export
simpleShape<-methods::setClass(Class = 'simpleShape',
                                contains = 'list')


#' @export
methods::setMethod('initialize',
                   'simpleShape',
                   function(.Object,
                            centroid = sf::st_sfc(NULL),
                            stem = sf::st_sfc(NULL),
                            branch = sf::st_sfc(NULL),
                            outline = sf::st_sfc(NULL),
                            ...){
                     .Object <- methods::callNextMethod(.Object, ...)
                     .Object <- list(centroid = centroid,
                                     stem = stem,
                                     branch = branch,
                                     outline = outline)
                     attr(.Object,'class')<-'simpleShape'
                     attr(attr(.Object,'class'),'package')<-'barbieHistologist'
                     return(.Object)
                   })





#' Object instantiation
#' 
#' The method bh_create cast a real shape or cell from a prototype.
#' @param x object of class genericShape or cellPrototype
#' @param lox numeric.
#' @param loy numeric, only for cell creation lox and loy are the coordinate
#'   where to drop the cell
#' @return a list containing an istantiation of a genericShape or an 
#'   an object of class cell.
#' @export
methods::setGeneric(name = 'bh_create',
                    signature = 'x',
                    function(x,...){})

methods::setMethod(f = 'bh_create',
                   signature = 'genericShape',
                   definition = function(x,...){
                     
                     MA <- rnorm(n = 1, mean = x@majorAxis[1],sd = x@majorAxis[2])
                     MI <- rnorm(n = 1, mean = x@minorAxis[1],sd = x@minorAxis[2])
                     RO <- rnorm(n = 1, mean = x@roundness[1],sd = x@roundness[2])
                     
                     if (x@fixedArms) NR<-x@nArms else NR <- rpois(n = 1, lambda = x@nArms)
                     
                     if (NR<2) NR<-2
                     if (NR%%2==1) NR<-NR+1
                     # AE <- rnorm(n = NR, mean = MA, sd = x@armExt[2]) * x@armExt[1]
                     AE <- rnorm(n = NR, mean = x@armExt[1], sd = x@armExt[2]) *MA 
                     
                     AW <- x@armElbow
                     AS <- x@armSwing
                     AT <- rnorm(n = NR, mean = x@armTrig[1],sd = x@armTrig[2])
                     
                     if (!is.null(x@orientation)) ang1<-x@orientation %% 360 else ang1<-sample(1:360,1)
                     
                     ang2<-ang1+180
                     angSpread<-90*RO
                     
                     anglLeft<-sample(((-angSpread):(angSpread))+ang1,round(NR/2))*pi/180
                     anglRight<-sample(((-angSpread):(angSpread))+ang2,round(NR/2))*pi/180
                     angl<-sort(c(anglLeft,anglRight))
                     
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
                         m1<--1/(y/x)
                         q1<-y-(m1*x)
                         
                         if (AT[i]<0) DD<-(AW-ii)*(-AT[i]) else DD<-ii*AT[i]
                         
                         smpAng<-atan2(sin(newAngl),cos(newAngl)) + pi
                         if ((smpAng>=0 & smpAng<pi)) drc<-(-1) else drc<-1
                         
                         if (!is.infinite(m1)){
                           x1<-x+drc*sqrt(DD/(1+m1^2))
                           x2<-x-drc*sqrt(DD/(1+m1^2))
                           y1<-m1*x1+q1
                           y2<-m1*x2+q1
                         } else {
                           x1 <- x
                           x2 <- x
                           y1 <- drc*DD
                           y2 <- -drc*DD
                         }
                         
                         out<-list(stem=matrix(c(x,y),ncol=2,byrow = T),
                                   branch=matrix(c(x1,y1,x,y,x2,y2),ncol=2,byrow=T))
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
                       
                       lx<-lapply(x$branch,function(y) y[c(T,F,F),])
                       lx <- do.call(rbind,lx)
                       rx<-lapply(x$branch,function(y) y[c(F,F,T),])
                       rx<-do.call(rbind,rev(rx))
                       out<-rbind(lx,rx)
                       # if (sum(lx[,2])>0) {
                       #   out<-rbind(lx,rx)} else {
                       #     out<-rbind(rx[nrow(rx):1,],lx[nrow(lx):1,])}
                     })
                     outLine<-do.call(rbind,outLine)
                     outLine<-rbind(outLine,outLine[1,])
                     outLine<-sf::st_polygon(list(outLine))
                     outLine<-sf::st_buffer(outLine,0)
                     # outLine<-sf::st_buffer(outLine,MA/10)
                     # outLine<-sf::st_buffer(outLine,-MA/10)
                     if (!sf::st_is_valid(outLine)){
                       outLine<-sf::st_make_valid(outLine)}
                     
                     # out<-list(centroid = sf::st_sfc(sf::st_point(c(0,0),dim = 'XYZ')),
                     #           stem = sf::st_sfc(stm),
                     #           branch = sf::st_sfc(brnch),
                     #           outline = sf::st_sfc(outLine))
                     
                     out<-new('simpleShape',centroid = sf::st_sfc(sf::st_point(c(0,0))),
                               stem = sf::st_sfc(stm),
                               branch = sf::st_sfc(brnch),
                               outline = sf::st_sfc(outLine))
                     
                     
                     return(out)
                   })




#' Define a new shape
#' 
#' helper function to create a new shape prototype.
#' @param majorAxis numeric, mean and sd pair.
#' @param minorAxis numeric, mean and sd pair.
#' @param nArms numeric, number of arms to be produced. This value is passed to
#'   rpois to obtain the actual number of arms.
#' @param armExt numeric, mean and sd pair. Proportion of the major axis.
#'   E.g. if major axis is 10 AU and arm extent is 1.5 the total lenght of 
#'   the arm will be approx 15 AU with the branching starting after 10 AU.
#' @param armElow numeric, number of branching points.
#' @param armSwing numeric, arm swinging extent in degree.
#' @param armTrig numeric, mean and sd pair. Determine the behaviour of branches.
#'   A positive value determine the branches to get longer towards the edge,
#'   a negative value makes branches to get smaller.
#' @return An object of class genericShape
#' @export
bh_defineShape<-function(majorAxis = NULL,
                         minorAxis = NULL,
                         roundness = NULL,
                         nArms = NULL,
                         fixedArms = NULL,
                         orientation = NULL,
                         armExt = NULL,
                         armElbow = NULL,
                         armSwing = NULL,
                         armTrig = NULL){
  
  if (is.null(majorAxis)) stop('define major axis')
  if (length(majorAxis)!=2) stop('major axis need a mean and sd')
  if (is.null(minorAxis)) stop('define minor axis')
  if (length(minorAxis)!=2) stop('minor axis need a mean and sd')
  if (is.null(roundness)) stop('define roundness')
  if (length(roundness)!=2) stop('roundness need a mean and sd')
  if (is.null(nArms)) stop('define number of arms')
  if (length(nArms)!=1) stop('number of arms must be a single value')
  if (is.null(fixedArms)) stop('fixed arms needs T/F')
  if (is.null(armExt)) stop('define arm extent')
  if (length(armExt)!=2) stop('arm extension need a mean and sd')
  if (is.null(armElbow)) stop('define number of arms breaks')
  if (length(armElbow)!=1) stop('arm elbow need a single value')
  if (is.null(armSwing)) stop('define a swing value')
  if (length(armSwing)!=1) stop('arm swing need a single value')
  if (is.null(armTrig)) stop('define the arm behaviour')
  if (length(armTrig)!=2) stop('arm behavior need a mean and sd')
  
  new('genericShape',
      majorAxis = majorAxis,
      minorAxis = minorAxis,
      roundness = roundness,
      nArms = nArms,
      fixedArms = fixedArms,
      orientation = orientation,
      armExt = armExt,
      armElbow = armElbow,
      armSwing = armSwing,
      armTrig = armTrig)
}