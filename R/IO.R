#' XY list
#' Wrapper around rasterToPoint, produces a XYZ list.
#' @param tissue tissue object
#' @return array
#' @export
bh_asXYZ<-function(tissue=NULL){
  if (is.null(tissue)) stop('give a tissue')
  
  out<-raster::rasterToPoints(tissue)
  return(out)
}

#' Geometry data frame
#' Wrapper around st_sf, produces a data frame with geometries.
#' @param cells a list of cells as produced by populate().
#' @return array
#' @export
bh_asSFC<-function(cells=NULL){
  
  if (is.null(cells)) stop('give cells as list')
  if (is.null(names(cells))) stop('names to cells')
  
  whatCells<-names(cells)
  out<-lapply(1:length(whatCells),function(wc){
    
    if (length(cells[[wc]])==0) {
      out<-sf::st_sf(cell = whatCells[[wc]],
                     ID = paste0(whatCells[[wc]],'_0'),
                     compartment=c('nucleus','cytoplasm','organelle'),
                     geometry=rep(sf::st_sfc(NULL),3))
      return(out)
    }
    
    out<-lapply(1:length(cells[[wc]]),function(ic){

      newCell<-c(cells[[wc]][[ic]]@nucleus$outline,
                 cells[[wc]][[ic]]@cytoplasm$outline,
                 cells[[wc]][[ic]]@organelle$outline)

            newcell<-sf::st_sf(cell=whatCells[[wc]],
                         ID = paste0(whatCells[[wc]],'_',ic),
                         compartment=c('nucleus','cytoplasm','organelle'),
                         geometry=newCell)
    })
    return(out)
  })
  
  out<-do.call(dplyr::bind_rows,out)
  seqid<-1:length(unique(out$ID))
  seqMultiplier<-rle(as.character(out$ID))
  repV<-Vectorize(rep)
  newIDs<-repV(times = seqMultiplier$lengths, x = seqid)
  out<-dplyr::bind_cols(out,seqid = newIDs)
  return(out)
}


#' XY list
#' Wrapper around st_sf, produces a data frame with geometries.
#' @param sfc a sf object
#' @param filePath character, full path and file name. Extension .sqlite
#'   will be added if missing.
#' @return array
#' @export
bh_saveSFC<-function(sfc = NULL,
                     filePath = NULL){
  
  if (is.null(sfc)) stop('give cells as list')
 if (is.null(filePath)) stop('give a full path file name')
  if (!grepl('.sqlite',filePath)) filePath<-paste0(filePath,'.sqlite')
  
  out<-try(sf::st_write(sfc,filePath,append=F))
  if (inherits(out,'try-error')) return(1) else return(0)
}


#' Save Tiff image
#' Wrapper around raster::writeRaster
#' @param raster a Raster Layer
#' @param filePath character, full path and file name. Extension .tiff
#'   will be added if missing.
#' @export
bh_saveTiff<-function(raster = NULL,
                     filePath = NULL){
  
  if (is.null(raster)) stop('give raster')
  if (is.null(filePath)) stop('give a full path file name')
  if (!grepl('.tif',filePath)) filePath<-paste0(filePath,'.tif')
  
  out<-try(raster::writeRaster(raster,filePath,overwrite=T))
  if (inherits(out,'try-error')) return(1) else return(0)
}