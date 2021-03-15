#' @include genericMarker.R
#' @include genericShape.R
#' 
#' 
#' @export
cellPrototype<-methods::setClass(Class = 'cellPrototype',
                                slots = c(cytoplasm = 'genericShape_NULL',
                                          nucleus = 'genericShape_NULL',
                                          organelle = 'genericShape_NULL',
                                          markers = 'list'),
                                prototype = list( cytoplasm = NULL,
                                                  nucleus = NULL,
                                                  organelle = NULL,
                                                  markers = list()))

#' @export
cell<-methods::setClass(Class = 'cell',
                                 slots = c(cytoplasm = 'list',
                                           nucleus = 'list',
                                           organelle = 'list',
                                           markers = 'list'),
                                 prototype = list( cytoplasm = list(),
                                                   nucleus = list(),
                                                   organelle = list(),
                                                   markers = list()))


methods::setMethod(f = 'bh_create',
                   signature = 'cellPrototype',
                   definition = function(x,
                                         lox,
                                         loy,
                                         ...){
                     # browser()
                     if (!is.null(x@cytoplasm)) {
                       cytoplasm<-bh_create(x@cytoplasm)
                       cytoplasm$stem<-cytoplasm$stem+c(lox,loy)
                       cytoplasm$branch<-lapply(cytoplasm$branch,'+',c(lox,loy))
                       cytoplasm$outline<-cytoplasm$outline+c(lox,loy)
                       } else cytoplasm<-list()
                     
                     if (!is.null(x@nucleus)) {
                       nucleus<-bh_create(x@nucleus)
                       nucleus$stem<-nucleus$stem+c(lox,loy)
                       nucleus$branch<-lapply(nucleus$branch,'+',c(lox,loy))
                       nucleus$outline<-nucleus$outline+c(lox,loy)
                       } else nucleus<-list()
                
                     if (!is.null(x@organelle)) {
                       organelle<-bh_create(x@organelle)
                       organelle$stem<-organelle$stem+c(lox,loy)
                       organelle$branch<-lapply(organelle$branch,'+',c(lox,loy))
                       organelle$outline<-organelle$outline+c(lox,loy)
                       } else organelle<-list()
                       
                     
                     newCell<-new('cell',
                                  cytoplasm = cytoplasm,
                                  nucleus = nucleus,
                                  organelle = organelle,
                                  markers = x@markers)
                     
                   })