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
                     attr(class(out),'package')<-getPackageName()
                     return(out)
                   })

#' Define a new tissue
#' 
#' helper function to create a new tissue.
#' @param coords numeric, x and y limits.
#' @param resolution numeric, resolution
#' @param bg numeric, value to use in the background
#' @param markers list, makers.
#' @return An object of class tissue (a rasterStack).
#' @export
bh_defineTissue<-function(coords = NULL,
                          resolution = NULL,
                          bg = NULL,
                          markers = NULL){
  
  if (is.null(coords)) stop('give coords')
  if (length(coords)!=4) stop('coords needs four values')
  if (is.null(resolution)) {
    resolution<-c(1,1)
    warning('resolution is not specified... assuming 1')} else {
      if (length(resolution)<2) {
        resolution<-c(resolution,resolution)
        warning('single value for resolution... assuming the same for both axis')}
    }
  if (is.null(bg)) {
    bg<-0
    warning('background is not specified... assuming 0')}
  if (is.null(markers)) stop('give markers')
  if (!is.list(markers)) stop('wrap markers in a list')
    
  new('tissue',
      xmn = coords[1],
      xmx = coords[2],
      ymn = coords[3],
      ymx = coords[4],
      bg = bg,
      xRes = resolution[1],
      yRes = resolution[2],
      markers = markers)
}
                        




#' @export
takePicture<-function(tissue = NULL,
                      cells = NULL){
  
  if (is.null(cells)) stop('need cells') else {
    if (!is.list(cells)) stop('need a list of real cells')
  }
  
  for (cell in cells){
    cellTiles<-list(cytoplasm = exactextractr::exact_extract(tissue,sf::st_sfc(cell@cytoplasm$outline),include_xy=T)[[1]][,c('x','y','coverage_fraction')],
                  nucleus = exactextractr::exact_extract(tissue,sf::st_sfc(cell@nucleus$outline),include_xy=T)[[1]][,c('x','y','coverage_fraction')])
    for (i in 1:length(cell@markers)){
      compartment<-cell@markers[[i]]@compartment
      lyr<-cell@markers[[i]]@Rname
      newCall<-rlang::call_modify(cell@markers[[i]]@pattern,xy=cellTiles[[compartment]][,1:2])
      # newTiles<-cell@markers[[lyr]]@pattern(xy=cellTiles[[compartment]][,1:2]) 
      newTiles<-eval(newCall)
      newTiles[newTiles[,3]<0,3]<-0
      tissue[[lyr]][raster::cellFromXY(tissue[[lyr]],newTiles[,1:2])]<-newTiles[,3]*cellTiles[[compartment]][,3]
    }
  }
  return(tissue)
}