#' @import raster

NULL

#' Tissue
#' 
#' A tissue is essentialy a RasterStack object, for more info **[raster::raster]**
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
                       raster::crs(out)<-NA
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
#' @param coords numeric vector, x and y limits... like (xMin,xMax,yMin,yMax).
#' @param resolution numeric vector, resolution on x and y.
#' @param bg numeric scalar, value to use in the background
#' @param markers list, makers.
#' @return An object of class **[tissue]**.
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



#' Define a new tissue
#' 
#' helper function to create a new tissue.
#' @param coords numeric vector, x and y limits... like (xMin,xMax,yMin,yMax).
#' @param resolution numeric vector, resolution on x and y.
#' @param bg numeric scalar, value to use in the background
#' @param markers list, makers.
#' @return An object of class **[tissue]**.
#' @export
bh_engrave<-function(tissue = NULL,
                     cells = NULL){
  
  if (is.null(cells)) stop('need cells') else {
    if (!is.list(cells)) stop('need a list of real cells')
  }
  
  for (cell in cells){
    
    cellTiles<-sapply(c('cytoplasm',
                        'nucleus',
                        'organelle'),
                      function(comp){
                        
                        newComp<-slot(cell,comp)$outline
                        
                        if (length(newComp)!=0) {
                          if (!sf::st_is_empty(newComp) & sf::st_is(newComp,c('POLYGON','MULTIPOLYGON'))){
                            cellTiles<-exactextractr::exact_extract(tissue,sf::st_sf(newComp),include_xy=T)
                            cellTiles<-do.call(rbind,cellTiles)
                            cellTiles<-cellTiles[,c('x','y','coverage_fraction')]
                          }
                        }
                      },simplify = F,USE.NAMES = T)
    
    for (i in 1:length(cell@markers)){
      compartment<-cell@markers[[i]]@compartment
      lyr<-cell@markers[[i]]@Rname
      if (!is.null(cellTiles[[compartment]])){
        if (nrow(cellTiles[[compartment]])!=0){
          newCall<-rlang::call_modify(cell@markers[[i]]@pattern,xy=cellTiles[[compartment]][,1:2])
          newTiles<-eval(newCall)
          newTiles[newTiles[,3]<0,3]<-0
          # tissue[[lyr]][raster::cellFromXY(tissue[[lyr]],newTiles[,1:2])]<-newTiles[,3]*cellTiles[[compartment]][,3]
          
          tissue[[lyr]][raster::cellFromXY(tissue[[lyr]],newTiles[,1:2])]<-tissue[[lyr]][raster::cellFromXY(tissue[[lyr]],newTiles[,1:2])] + newTiles[,3]*cellTiles[[compartment]][,3]
          
        }
      }
    }
  }
  return(tissue)
}

#' @export
bh_familyPicture<-function(tissue = NULL,
                           sfc = NULL,
                           compartment = 'cytoplasm'){
  
  if (is.null(tissue)) stop('give tissue')
  if (is.null(sfc)) stop ('give sfc')
  if (is.null(compartment)) stop('give compartment')
  
  if (!any(compartment %in% sfc$compartment)) stop('cannot find compartment')
  
  rasterID<-raster::raster(ext = raster::extent(tissue[[1]]),resolution = raster::res(tissue[[1]]),crs='',vals=0)
  rasterCov<-rasterID
  
  sfc<-sfc[sfc$compartment==compartment,]
  
  sfc_valid<-sfc[!sf::st_is_empty(sfc),]
  
  tiles<-exactextractr::exact_extract(x = rasterID,
                                      y = sfc_valid,
                                      include_xy=T)
  
  for (i in 1:length(tiles)){
    newID<-sfc_valid$seqid[i]
    
    tileID<-raster::cellFromXY(object = rasterID,
                               xy = tiles[[i]][,c('x','y')])
    rasterID[tileID]<-newID
    rasterCov[tileID]<-tiles[[i]][,'coverage_fraction']
  }
  
  out<-list(ID=rasterID,
            cov=rasterCov)
  out<-raster::stack(out)
  
  return(out)
}

#' @export
bh_saveFamily<-function(familyPicture = NULL,
                        familyName = NULL,
                        filePath = NULL){
  
  if (is.null(familyPicture)) stop('give picture')
  if (is.null(familyName)) stop('give name')
  if (is.null(filePath)) stop('give folder')
  if(!dir.exists(filePath)) stop('cannot  find folder')
  
  newDir<-file.path(filePath, familyName)
  if (!dir.exists(newDir)) dir.create(newDir)
  nms<-names(familyPicture)
  newStk<-vector(mode = 'list',length = length(nms))
  names(newStk)<-nms
  
  for (i in nms){
    newStk[[i]]<-raster::writeRaster(familyPicture[[i]],
                                     file.path(newDir,
                                               paste0(i,'.nc')),
                                     format = 'CDF',
                                     overwrite = T)
    raster::crs(newStk[[i]])<-''
  }
  newStk<-raster::stack(newStk)
  return(newStk)
}

#' @export
bh_loadFamily<-function(filePath = NULL){
  
  if (is.null(filePath)) stop('give folder')
  if(!dir.exists(filePath)) stop('cannot  find folder')
  
  newFiles<-list.files(filePath,
                       full.names = T)
  nms<-sub(pattern = '.nc',
           replacement = "",
           list.files(filePath, full.names = F))
  newStk<-vector(mode = 'list',length = length(newFiles))
  names(newStk)<-nms
  for (i in 1:length(newFiles)){
    newStk[[i]]<-raster::raster(newFiles[[i]])
    raster::crs(newStk[[i]])<-''
  }
  newStk<-raster::stack(newStk)
  return(newStk)
}

#' @export
bh_extractFamily<-function(family = NULL,
                           sfc = NULL,
                           sfcPrimaryKey = NULL){
  
  if (is.null(family)) stop('give family')
  if (is.null(sfc)) stop('give sfc')
  if (is.null(sfcPrimaryKey)) stop('give primary key')
  if (!is.character(sfcPrimaryKey)) stop('primary key must be character')
  
  sfc<-sfc[!sf::st_is_empty(sfc),]
  sfc<-sfc[sf::st_geometry_type(sfc)=='POLYGON' |
             sf::st_geometry_type(sfc)=='MULTIPOLYGON',]
  
  out<-exactextractr::exact_extract(x = family,
                                    y = sfc,
                                    include_xy=T)
  sfcClean<-sf::st_drop_geometry(sfc)
  
  out<-lapply(1:length(out),function(i){
    data.frame(sfcClean[i,sfcPrimaryKey],
               out[[i]])
  })
  
  out<-do.call(rbind.data.frame,out)
  names(out)[1]<-sfcPrimaryKey
  
  return(out)
}
