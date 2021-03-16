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

#' XY list
#' Wrapper around rasterToPoint, produces a XYZ list.
#' @param tissue tissue object
#' @return array
#' @export
bh_asSFC<-function(cells=NULL,
                   compartments=NULL){
  
  if (is.null(cells)) stop('give cells as list')
  if (is.null(names(cells))) stop('names to cells')
  if (is.null(compartments)) stop('specify a list of compartments')
  
  cellNames<-names(cells)
  out_cls<-lapply(cellNames,function(nc){
    cls<-cells[[nc]]
    out_cl<-lapply(cls,function(cl){
      out_cmp<-lapply(compartments,function(cmp){
        out<-slot(cl,cmp)$outline
        })
      out_cmp<-sf::st_multipolygon(out_cmp)
      })
    out_cl<-sf::st_sfc(out_cl)
    out_cl<-sf::st_sf(data.frame(label = nc,geom = out_cl,stringsAsFactors = F))
    })
  out_cls<-do.call(rbind.data.frame,out_cls)
}

  