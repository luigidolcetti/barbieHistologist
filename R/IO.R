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
      
      is_empty<-sapply(out_cmp,sf::st_is_empty,simplify = T,USE.NAMES = F)
      is_type<-as.character(sapply(out_cmp,sf::st_geometry_type,simplify = T,USE.NAMES = F))
      if (any(is_empty)){
        
        out_cmp<-out_cmp[!is_empty]}
        # out_cmp[is_empty]<-sf::st_polygon(
        #   list(
        #     matrix(c(0,0,1,0,1,1,0,1,0,0),ncol=2,byrow = T)))}
      if (any(is_type=='MULTIPOLYGON')){
        ww<-which(is_type=='MULTIPOLYGON')
        new_out_cmp<-list()
        new_out_cmp_i<-1
        for (i in ww){
          out_cmp[[i]]<-sf::st_cast(
            sf::st_sfc(out_cmp[[i]]),
            'POLYGON')
          for (ii in 1:length(out_cmp[[i]])){
            new_out_cmp[[new_out_cmp_i]]<-out_cmp[[i]][[ii]]
            new_out_cmp_i<-new_out_cmp_i+1
          }
        }
        out_cmp<-append(out_cmp[!(is_type=='MULTIPOLYGON')],new_out_cmp)
      }
      out_cmp<-sf::st_multipolygon(out_cmp)
      # if (any(do.call(sf::st_is_empty,out_cmp))){
      #   out_cmp<-sf::st_multipolygon(
      #     list(
      #       sf::st_polygon(
      #         list(
      #           matrix(c(0,0,1,0,1,1,0,1,0,0),ncol=2,byrow = T)))))
      # } else {
      # if (all(as.character(do.call(sf::st_geometry_type,out_cmp))=='POLYGON')){
      #   out_cmp<-sf::st_multipolygon(out_cmp)}
      # }
    })
    out_cl<-sf::st_sfc(out_cl)
    out_cl<-sf::st_sf(data.frame(label = nc,geom = out_cl,stringsAsFactors = F))
  })
  out_cls<-do.call(rbind.data.frame,out_cls)
}

