#' @include genericMarker.R
#' @include genericShape.R
#' 
#' 
#' @export
cellPrototype<-methods::setClass(Class = 'cellPrototype',
                                 slots = c(name = 'character',
                                           cytoplasm = 'genericShape_NULL',
                                           nucleus = 'genericShape_NULL',
                                           organelle = 'genericShape_NULL',
                                           markers = 'list'),
                                 prototype = list( name = 'aNameForThisCell',
                                                   cytoplasm = NULL,
                                                   nucleus = NULL,
                                                   organelle = NULL,
                                                   markers = list()))

methods::setMethod(f = 'names',
                   signature = 'cellPrototype',
                   definition = function(x){
                     slot(x,'name')
                   })

methods::setMethod(f = 'names<-',
                   signature = 'cellPrototype',
                   definition = function(x,value){
                     slot(x,'name')<-value
                   })


#' @export
cell<-methods::setClass(Class = 'cell',
                        slots = c(name = 'character',
                                  cytoplasm = 'list',
                                  nucleus = 'list',
                                  organelle = 'list',
                                  markers = 'list'),
                        prototype = list( name = 'aNameForThisCell',
                                          cytoplasm = list(),
                                          nucleus = list(),
                                          organelle = list(),
                                          markers = list()))

methods::setMethod(f = 'names',
                   signature = 'cell',
                   definition = function(x){
                     slot(x,'name')
                   })

methods::setMethod(f = 'names<-',
                   signature = 'cell',
                   definition = function(x,value){
                     slot(x,'name')<-value
                   })


methods::setMethod(f = 'bh_create',
                   signature = 'cellPrototype',
                   definition = function(x,
                                         lox,
                                         loy,
                                         shrinkage,
                                         setEscape,
                                         ...){
                     
                     
                     if (!is.null(x@cytoplasm)) {
                       cytoplasm<-try(bh_create(x@cytoplasm))
                       if (inherits(cytoplasm,'try-error')) return(0)
                       cytoplasm<-sapply(cytoplasm,'+',c(lox,loy),simplify = F,USE.NAMES = T)
                       newLo<-try(sf::st_point_on_surface(cytoplasm$outline))
                       if (inherits(newLo,'try-error')) return(0)
                     } else {
                       cytoplasm<-list()
                       cytoplasm$stem<-sf::st_sfc(NULL)
                       cytoplasm$branch<-sf::st_sfc(NULL)
                       cytoplasm$outline<-sf::st_sfc(NULL)}
                     
                     if (!is.null(x@nucleus)) {
                       nucleus<-try(bh_create(x@nucleus))
                       if (inherits(nucleus,'try-error')) return(0)
                       nucleus<-sapply(nucleus,'+',newLo,simplify = F,USE.NAMES = T)
                       cytoAndNucleus<-try(sf::st_union(cytoplasm$outline,nucleus$outline))
                       if (inherits(cytoAndNucleus,'try-error')) return(0)
                       area_cytoplasm<-sf::st_area(cytoplasm$outline)
                       area_cytoAndNucleus<-sf::st_area(cytoAndNucleus)
                       safeEscape<-setEscape
                       while (length((sf::st_within(nucleus$outline,
                                                    cytoplasm$outline)[[1]]))==0){
                       
                         area_nucleus<-sf::st_area(nucleus$outline)
                         # browser()
                         Aratio<-area_cytoplasm/area_cytoAndNucleus*shrinkage
                         # if (Aratio>1) Aratio<-1/Aratio
                         cNuc<-sf::st_centroid(nucleus$outline)
                         nucleus<-sapply(nucleus,'-',cNuc,simplify = F,USE.NAMES = T)
                         nucleus<-sapply(nucleus,'*',Aratio,simplify = F,USE.NAMES = T)
                         nucleus<-sapply(nucleus,'+',cNuc,simplify = F,USE.NAMES = T)
                         cytoAndNucleus<-try(sf::st_union(cytoplasm$outline,nucleus$outline))
                         if (inherits(cytoAndNucleus,'try-error')) return(0)
                         area_cytoAndNucleus<-sf::st_area(cytoAndNucleus)
                         safeEscape <- safeEscape-1
                         if (safeEscape<0) return(0)
                       }
                     } else {
                       nucleus<-list()
                       nucleus$stem<-sf::st_sfc(NULL)
                       nucleus$branch<-sf::st_sfc(NULL)
                       nucleus$outline<-sf::st_sfc(NULL)}
                     
                     if (!is.null(x@organelle)) {
                       organelle<-try(bh_create(x@organelle))
                       if (inherits(organelle,'try-error')) return(0)
                       newLo<-sf::st_point_on_surface(cytoNotNucleus)
                       organelle<-sapply(organelle,'+',c(lox,loy),simplify = F,USE.NAMES = T)
                     } else {
                       organelle<-list()
                       organelle$stem<-sf::st_sfc(NULL)
                       organelle$branch<-sf::st_sfc(NULL)
                       organelle$outline<-sf::st_sfc(NULL)}
                     
                     newCell<-new('cell',
                                  name = x@name,
                                  cytoplasm = cytoplasm,
                                  nucleus = nucleus,
                                  organelle = organelle,
                                  markers = x@markers)
                     
                   })

#' Define a new cell
#' 
#' helper function to create a new cell prototype.
#' @param cytoplasm genericShape.
#' @param nucleus genericShape.
#' @param organelle genericSape. (not implemented yet)
#' @param markers list, makers.
#' @return An object of class cellPrototype
#' @export
bh_defineCell<-function(name = NULL,
                        cytoplasm = NULL,
                        nucleus = NULL,
                        organelle = NULL,
                        markers = NULL){
  
  if (is.null(cytoplasm)) stop('provide a name for this cell')
  if (is.null(cytoplasm)) stop('a cell needs at least a cytoplasm')
  if (is.null(markers)) stop('give some markers')
  if (!is.list(markers)) stop('wrap markers in a list')
  
  new('cellPrototype',
      name = name,
      cytoplasm = cytoplasm,
      nucleus = nucleus,
      organelle = organelle,
      markers = markers)
}

#'@export
bh_clone<-function(cell,
                   lox,
                   loy){
  for (comp in c('cytoplasm','nucleus','organelle')){
    oldConf<-slot(cell,comp)
    for (st in c('stem','branch','outline')){
      # 
      oldConfSt<-oldConf[[st]]
      cntrd<-sf::st_centroid(oldConfSt)
      newConf<-oldConfSt-cntrd+c(lox,loy)
      oldConf[[st]]<-newConf}
    slot(cell,comp)<-oldConf}
  return(cell)}
