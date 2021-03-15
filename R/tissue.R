#' @import raster
#' 
#' 
#' 
#' @export
tissue<-methods::setClass(Class = 'tissue',
                          contains = 'RasterStack')

#' @export
methods::setMethod('initialize',
                   'tissue',
                   function(.Object,
                            xmn = 0,
                            xmx = 500,
                            ymn = 0,
                            ymx = 500,
                            bg = 0,
                            xRes = 1,
                            yRes = 1,
                            markers = list(),
                            ...){
                     
                     .Object <- methods::callNextMethod(.Object, ...)
                     
                     if (length(markers)==0) stop ('markers is empty')
                     newRst<-lapply(markers,function(mrkr){
                       # layerName<-mrkr@Rname
                       resolution<-c(xRes,yRes)
                       out<-raster::raster(xmn=xmn,
                                      ymn=ymn,
                                      xmx=xmx,
                                      ymx=ymx,
                                      resolution=resolution,
                                      vals=bg)
                       names(out)<-mrkr@Rname
                       return(out)
                     })
                     
                     out<-raster::stack(newRst)
                  
                     
                     class(out)<-c('tissue')
                     return(out)
                   })

#' @export
takePicture<-function(tissue = NULL,
                      cells = NULL){
  
  if (is.null(cells)) stop('need cells') else {
    if (!is.list(cells)) stop('need a list of real cells')
  }
  
  for (cell in cells){
    cellTiles<-list(cytoplasm = exactextractr::exact_extract(tissue,sf::st_sfc(cell@cytoplasm$outline),include_xy=T)[[1]][,c('x','y','coverage_fraction')],
                  nucleus = exactextractr::exact_extract(tissue,sf::st_sfc(cell@nucleus$outline),include_xy=T)[[1]][,c('x','y','coverage_fraction')])
    for (lyr in 1:length(cell@markers)){
     
      compartment<-cell@markers[[lyr]]@compartment
      newCall<-rlang::call_modify(cell@markers[[lyr]]@pattern,xy=cellTiles[[compartment]][,1:2])
      # newTiles<-cell@markers[[lyr]]@pattern(xy=cellTiles[[compartment]][,1:2]) 
      newTiles<-eval(newCall)
      tissue[[lyr]][raster::cellFromXY(tissue[[lyr]],newTiles[,1:2])]<-newTiles[,3]*newTiles[,3]
    }
  }
  return(tissue)
}