#' Populate by clipping
#' 
#' Populate comes in two versions, one is by clipping, meaning that cells are 
#'   formed without any knowledge of the surrounding cells and the overlaps will be 
#'   just clipped at the end. On the contrary, **[bh_populate_byInteract()]** creates
#'   new cells considering neighboring cells.
#' @param cellPrototype list of objects of class cellPrototype
#' @param proportion numeric vector, proportion (0-1) of that specific cell
#' @param maxCloning numeric scalar, max of times the same cell can be cloned 
#'   around the area (may speed up the process).
#' @param tissue object of class **[tissue]**
#' @param cropToMeasure logical scalar, should cell outline cropped to the extent
#'   of tissue?
#' @param require_cytoplasm logical scalar, is cytoplasm otline required? cell 
#'   created without cytoplasm will be discarded. The same apply to require_nucleus
#'   and require_organelle (not implemented at the moment).
#' @param require_nucleus logical scalar.
#' @param require_organelle logical scalar.
#' @param areaTresh numeric scalar, proportion (0-1) of tissue covered by cells.
#'   Limit is calculated on the number of pixels that has either complete or 
#'   partial coverage.
#' @return a list containing cells spatially distributed.
#' @export
bh_populate_byClip<-function(cellPrototype = NULL,
                             proportion = NULL,
                             maxCloning = 10,
                             tissue = NULL,
                             cropToMesure=T,
                             require_cytoplasm = T,
                             require_nucleus = T,
                             require_organelle = F,
                             areaTresh = NULL){
  
  areaTot<-0
  cloningI<-maxCloning+1
  lim <-raster::extent(tissue)
  res <-raster::res(tissue)
  minResHalf<-min(res)
  cellNames<-sapply(cellPrototype,names,simplify = T,USE.NAMES = F)
  
  possibleXY<-raster::raster(vals=0,ext=lim,resolution=res)
  
  valFreq<-raster::freq(possibleXY)
  val1<-valFreq[valFreq[,1]==1,2]
  if (length(val1)==0) val1<-0
  val0<-valFreq[valFreq[,1]==0,2]
  if (length(val0)==0) val0<-0
  
  
  newCellList<-sapply(1:length(cellNames),function(clnms){
    vector('list',round(val0*proportion[clnms]))},simplify = F,USE.NAMES = F)
  
  names(newCellList)<-cellNames
  newCellType_counter<-rep(1,length(cellNames))
  names(newCellType_counter)<-cellNames
  
  while((val1/(val0+val1))<areaTresh){
    
    if (cloningI>maxCloning){
      
      newCellType<-sample(1:length(cellPrototype),1,prob = proportion)
      
      position<-raster::Which(possibleXY==0,cells=T)
      position<-sample(position,1)
      position<-try(raster::xyFromCell(object = tissue,cell = position),silent = T)
      
      if (inherits(position,'try-error')){
        cat('X--> bad position\n')
        next}
      
      newCellClone<-try(bh_create(cellPrototype[[newCellType]],
                                  lox = position[,'x'],
                                  loy = position[,'y']),silent = T)
      
      if (inherits(newCellClone,'try-error')){
        cat('XX-> bad cell\n')
        cloningI<-maxCloning+1
        next}
      
      if (class(newCellClone)!='cell') {
        cat('XX-> bad cell\n')
        cloningI<-maxCloning+1
        next}
      
      newCell<-newCellClone
      
      newComp<-try(sapply(c('cytoplasm',
                            'nucleus',
                            'organelle'),
                          function(cmp){
                            
                            newComp<-slot(newCell,cmp)$outline
                          },USE.NAMES = T,simplify = F),silent = T)
      
      if (inherits(newComp,'try-error')){ 
        cat('x-->bad cell\n')
        cloningI<-maxCloning+1
        next
      }
      
      unionCell<-try(sf::st_union(do.call(c,newComp),crs = 'NA',precision = minResHalf),silent = T)
      
      if (inherits(unionCell,'try-error')){
        cat('xx->bad union cell\n')
        cloningI<-maxCloning+1
        next
      }
      
      area_unionCell<-try(sf::st_area(unionCell),silent = T)
      
      if (inherits(area_unionCell,'try-error')){
        cat('XX-> bad area union\n')
        cloningI<-maxCloning+1
        next}
      
      area_cell<-try(sf::st_area(slot(newCell,'cytoplasm')$outline),silent=T)
      
      if (inherits(area_cell,'try-error')){
        cat('X--> bad area cytoplasm\n')
        cloningI<-maxCloning+1
        next}
      
      
      if (area_unionCell>area_cell){
        cat('xx->bad area cell\n')
        cloningI<-maxCloning+1
        next
      }
      
      cloningI <-1
      cat('---> created new cell clone\n')
      
    } else {
      
      position<-raster::Which(possibleXY==0,cells=T)
      position<-sample(position,1)
      position<-raster::xyFromCell(object = tissue,cell = position)
      
      newCell<-try(bh_clone(cell = newCellClone,
                            lox = position[,'x'],
                            loy = position[,'y']))
      
      if (inherits(newCell,'try-error')){ 
        cat('x-->bad cloning\n')
        cloningI<-maxCloning+1
        next
      }
      
      newComp<-try(sapply(c('cytoplasm',
                            'nucleus',
                            'organelle'),
                          function(cmp){
                            
                            newComp<-slot(newCell,cmp)$outline
                          },USE.NAMES = T,simplify = F),silent = T)
      
      if (inherits(newComp,'try-error')){ 
        cat('x-->bad cell\n')
        cloningI<-maxCloning+1
        next
      }
      
      unionCell<-try(sf::st_union(do.call(c,newComp),crs = 'NA',precision = minResHalf),silent = T)
      
      if (inherits(unionCell,'try-error')){ 
        cat('x-->bad cell\n')
        cloningI<-maxCloning+1
        next
      }
      
      cloningI <- cloningI+1
    }
    
    xyCoverage<-try(exactextractr::exact_extract(possibleXY,sf::st_sfc(unionCell),include_xy=T))
    
    if (inherits(xyCoverage,'try-error')){ 
      cat('x-->bad coverage\n')
      next
    }
    
    xyCoverage<-try(raster::cellFromXY(possibleXY,xyCoverage[[1]][,c('x','y')]))
    
    if (inherits(xyCoverage,'try-error')){ 
      cat('x-->bad coverage\n')
      next
    }
    
    possibleXY[xyCoverage]<-1
    newCellList[[newCellType]][[newCellType_counter[newCellType]]]<-newCell
    newCellType_counter[newCellType]<-newCellType_counter[newCellType]+1
    
    valFreq<-raster::freq(possibleXY)
    val1<-valFreq[valFreq[,1]==1,2]
    if (length(val1)==0) val1<-0
    val0<-valFreq[valFreq[,1]==0,2]
    if (length(val0)==0) val0<-0
    
    cat(paste0('area ratio: ',scales::label_percent(accuracy = 0.1)(((val1)/(val1+val0))),'\n'))
    
  }  
  
  whichToGetRid<-lapply(newCellList,function(x){
    sapply(x,is.null,simplify = T)})
  
  out_newCellList<-lapply(1:length(newCellList),function(x){
    newCellList[[x]][!whichToGetRid[[x]]]})
  
  names(out_newCellList)<-names(newCellList)
  
  bodyOnly<-lapply(names(out_newCellList),function(Ncls){
    
    cls<-out_newCellList[[Ncls]]
    out<-lapply(cls,function(cl){
      slot(cl,'cytoplasm')$outline})
    out<-do.call(c,out)
    
    if (is.null(out)) out<-sf::st_sfc(NULL)
    out<-sf::st_sf(data.frame(ID=1:length(out),
                              cell=Ncls,
                              geom=out))
    return(out)
  })
  
  bodyOnly<-do.call(dplyr::bind_rows,bodyOnly)
  if (!all(sf::st_is_valid(bodyOnly))) bodyOnly<-sf::st_make_valid(bodyOnly)
  molds<-try(sf::st_difference(sf::st_buffer(bodyOnly,0)),silent = T)
  
  if (inherits(molds,'try-error')) stop('populate failed retry')
  
  
  if (!all(sf::st_is_valid(molds))) {
    molds<-sf::st_make_valid(molds)
    
  }
  
  if (cropToMesure) molds<-sf::st_crop(molds,lim)
  
  if (any(sf::st_is_valid(molds))) {
    molds<-sf::st_make_valid(molds)
    
  }
  
  for (i in 1:nrow(molds)){
    mold<-molds[i,'geometry']
    index<-unlist(sf::st_drop_geometry(molds[i,'ID']))
    cell<-as.character(unlist(sf::st_drop_geometry(molds[i,'cell'])))
    for (compartment in c('cytoplasm','nucleus','organelle')){
      oldGeom<-slot(out_newCellList[[cell]][[index]],compartment)$outline
      
      if (!sf::st_is_valid(oldGeom)) oldGeom<-sf::st_make_valid(oldGeom)
      
      newGeom<-try(sf::st_intersection(oldGeom,mold),silent = T)
      
      if (inherits(newGeom,'try-error')) {
        newGeom<-sf::st_sfc(NULL)}
      
      if (length(newGeom)==0) newGeom<-sf::st_sfc(NULL) else{
        if (any(class(newGeom)=='sfc_GEOMETRYCOLLECTION' |
                class(newGeom)=='sfc_MULTIPOLYGON')) {
          
          newGeom<-newGeom[sf::st_is(newGeom,c('POLYGON','MULTIPOLYGON')),]
          newGeom<-sf::st_cast(newGeom,'POLYGON')
          stArea<-sf::st_area(newGeom)
          maxStArea<-which.max(stArea)
          newGeom<-newGeom[maxStArea,]
        }
      }
     
      
      if (length(sf::st_geometry_type(newGeom))==0) newGeom<-sf::st_sfc(NULL)
      
      if (!sf::st_is_valid(newGeom)) newGeom<-sf::st_make_valid(newGeom)
      
      slot(out_newCellList[[cell]][[index]],compartment)$outline<-newGeom
    }
  }
  
  whichToGetRid<-lapply(out_newCellList,function(x){
    sapply(x,function(x){
      status_cy<-sf::st_is_empty(x@cytoplasm$outline)
      status_nu<-sf::st_is_empty(x@nucleus$outline)
      status_or<-sf::st_is_empty(x@organelle$outline)
      type_cy<-sf::st_geometry_type(x@cytoplasm$outline)
      type_nu<-sf::st_geometry_type(x@nucleus$outline)
      type_or<-sf::st_geometry_type(x@organelle$outline)
      valid_cy<-sf::st_is_valid(x@cytoplasm$outline)
      valid_nu<-sf::st_is_valid(x@nucleus$outline)
      valid_or<-sf::st_is_valid(x@organelle$outline)
      if (status_cy & require_cytoplasm) return(T)
      if (status_nu & require_nucleus) return(T)
      if (status_or & require_organelle) return(T)
      if (!status_cy & type_cy != 'POLYGON' & require_cytoplasm) return(T)
      if (!status_nu & type_nu != 'POLYGON' & require_nucleus) return(T)
      if (!status_or & type_or != 'POLYGON' & require_organelle) return(T)
      if (!valid_cy | !valid_nu | !valid_or) return(T)
      if (is.null(x)) return(T)
      return(F)
    },simplify = T)
  })
  
  oldNames<-names(out_newCellList)
  
  out_newCellList<-lapply(1:length(out_newCellList),function(x){
    out_newCellList[[x]][!whichToGetRid[[x]]]})
  
  names(out_newCellList)<-oldNames
  
  return(out_newCellList)
}

#' Populate by interact
#' 
#' Populate comes in two versions, one is by clipping, meaning that cells are 
#'   formed without any knowledge of the surrounding cells and the overlaps will be 
#'   just clipped at the end. On the contrary, **[bh_populate_byInteract()]** creates
#'   new cells considering neighboring cells.
#' @param cellPrototype list of objects of class cellPrototype
#' @param proportion numeric vector, proportion [0-1] of that specific cell
#' @param maxCloning numeric scalar, max of times the same cell can be cloned 
#'   around the area (may speed up the process).
#' @param tissue object of class **[tissue]**
#' @param cropToMeasure logical scalar, should cell outline cropped to the extent
#'   of tissue?
#' @param require_cytoplasm logical scalar, is cytoplasm otline required? cell 
#'   created without cytoplasm will be discarded. The same apply to require_nucleus
#'   and require_organelle (not implemented at the moment).
#' @param require_nucleus logical scalar.
#' @param require_organelle logical scalar.
#' @param areaTresh numeric scalar, proportion (0-1) of tissue covered by cells.
#'   Limit is calculated on the number of pixels that has either complete or 
#'   partial coverage.
#' @return a list containing cells spatially distributed.
#' @export
bh_populate_byInteract<-function(cellPrototype = NULL,
                                 proportion = NULL,
                                 tissue = NULL,
                                 cropToMesure=T,
                                 require_cytoplasm = T,
                                 require_nucleus = T,
                                 require_organelle = F,
                                 areaTresh = NULL){
  
  areaTot<-0
  lim <-raster::extent(tissue)
  res <-raster::res(tissue)
  minResHalf<-min(res)
  cellNames<-sapply(cellPrototype,names,simplify = T,USE.NAMES = F)
  
  blob<-sf::st_sfc(NULL)
  bBox<-sf::st_sfc(sf::st_polygon(list(matrix(c(lim[1],lim[3],
                                                lim[1],lim[4],
                                                lim[2],lim[4],
                                                lim[2],lim[3],
                                                lim[1],lim[3]),
                                              ncol=2,
                                              byrow=T))))
  antiBlob<-sf::st_difference(bBox,blob)
  
  possibleXY<-raster::raster(vals=0,ext=lim,resolution=res)
  
  valFreq<-raster::freq(possibleXY)
  val1<-valFreq[valFreq[,1]==1,2]
  if (length(val1)==0) val1<-0
  val0<-valFreq[valFreq[,1]==0,2]
  if (length(val0)==0) val0<-0
  
  
  newCellList<-sapply(1:length(cellNames),function(clnms){
    vector('list',round(val0*proportion[clnms]))},simplify = F,USE.NAMES = F)
  
  names(newCellList)<-cellNames
  newCellType_counter<-rep(1,length(cellNames))
  names(newCellType_counter)<-cellNames
  
  while((val1/(val0+val1))<areaTresh){
    
    newCellType<-sample(1:length(cellPrototype),1,prob = proportion)
    
    position<-raster::Which(possibleXY==0,cells=T)
    position<-sample(position,1)
    
    position<-try(raster::xyFromCell(object = tissue,cell = position),silent = T)
    if (inherits(position,'try-error')){
      cat('X--> bad position\n')
      next}
    
    
    newCell<-try(bh_create(cellPrototype[[newCellType]],
                           lox = position[,'x'],
                           loy = position[,'y'],
                           constrained = 'outside',
                           constrainTo = list(outline=blob)),silent = T)
    if (inherits(newCell,'try-error')){
      cat('XX-> bad cell\n')
      next}
    
    if (require_cytoplasm){
      test<-slot(newCell,'cytoplasm')$outline
      if (length(test)==0) {cat('XX-> no cytoplasm\n'); next}
      if (is.null(test)) {cat('XX-> no cytoplasm\n'); next}
      if (sf::st_is_empty(test)) {cat('XX-> no cytoplasm\n'); next}
    }
    
    if (require_nucleus){
      test<-slot(newCell,'nucleus')$outline
      if (length(test)==0)  {cat('XX-> no nucleus\n'); next}
      if (is.null(test))  {cat('XX-> no nucleus\n'); next}
      if (sf::st_is_empty(test))  {cat('XX-> no nucleus\n'); next}
     
    }
    
    if (require_organelle){
      test<-slot(newCell,'organelle')$outline
      if (length(test)==0) {cat('XX-> no organelle\n'); next}
      if (is.null(test)) {cat('XX-> no organelle\n'); next}
      if (sf::st_is_empty(test)) {cat('XX-> no organelle\n'); next}
      
    }
    
    blob<-try(sf::st_union(blob,newCell@cytoplasm$outline))
    if (inherits(blob,'try-error')){
      cat('XX-> bad blob union\n')
      next}
    
    antiBlob<-try(sf::st_difference(bBox,blob))
    if (inherits(antiBlob,'try-error')){
      cat('XX-> bad blob union\n')
      next}
    
    xyCoverage<-try(exactextractr::exact_extract(possibleXY,blob,include_xy=T))
    
    if (inherits(xyCoverage,'try-error')){ 
      cat('x-->bad coverage\n')
      next
    }
    
    xyCoverage<-try(raster::cellFromXY(possibleXY,xyCoverage[[1]][,c('x','y')]))
    
    if (inherits(xyCoverage,'try-error')){ 
      cat('x-->bad coverage\n')
      next
    }
    
    possibleXY[xyCoverage]<-1
    newCellList[[newCellType]][[newCellType_counter[newCellType]]]<-newCell
    newCellType_counter[newCellType]<-newCellType_counter[newCellType]+1
    
    valFreq<-raster::freq(possibleXY)
    val1<-valFreq[valFreq[,1]==1,2]
    if (length(val1)==0) val1<-0
    val0<-valFreq[valFreq[,1]==0,2]
    if (length(val0)==0) val0<-0
    
    cat(paste0('area ratio: ',scales::label_percent(accuracy = 0.1)(((val1)/(val1+val0))),'\n'))
    
  }  
  
  whichToGetRid<-lapply(newCellList,function(x){
    sapply(x,is.null,simplify = T)})
  
  out_newCellList<-lapply(1:length(newCellList),function(x){
    newCellList[[x]][!whichToGetRid[[x]]]})
  
  names(out_newCellList)<-names(newCellList)
  
  if (cropToMesure) {
    for (i in names(out_newCellList)){
      for (ii in seq_along(out_newCellList[[i]])){
        if (!sf::st_is_empty(out_newCellList[[i]][[ii]]@nucleus$outline)){
          out_newCellList[[i]][[ii]]@nucleus$outline<-sf::st_crop(out_newCellList[[i]][[ii]]@nucleus$outline,lim)
          if (length(out_newCellList[[i]][[ii]]@nucleus$outline)==0) out_newCellList[[i]][[ii]]@nucleus$outline<-sf::st_sfc(NULL)
          }
        if (!sf::st_is_empty(out_newCellList[[i]][[ii]]@cytoplasm$outline)){
          out_newCellList[[i]][[ii]]@cytoplasm$outline<-sf::st_crop(out_newCellList[[i]][[ii]]@cytoplasm$outline,lim)
          if (length(out_newCellList[[i]][[ii]]@cytoplasm$outline)==0) out_newCellList[[i]][[ii]]@cytoplasm$outline<-sf::st_sfc(NULL)
        }
        if (!sf::st_is_empty(out_newCellList[[i]][[ii]]@organelle$outline)){
          out_newCellList[[i]][[ii]]@organelle$outline<-sf::st_crop(out_newCellList[[i]][[ii]]@organelle$outline,lim)
          if (length(out_newCellList[[i]][[ii]]@organelle$outline)==0) out_newCellList[[i]][[ii]]@organelle$outline<-sf::st_sfc(NULL)
        }
      }
    }
  }
  
  return(out_newCellList)
}

#' Save a population.
#' 
#' Just save a population created by either **[bh_populate_byClip]** or
#'   **[bh_populate_byInteract]**
#' @param x list of cells
#' @param file character scalar
#' @export
bh_savePopulation<-function(x = NULL,file = NULL){
  
  if (is.null(x)) stop('give a population')
  if (is.null(file)) stop('give a full path filename')
  
  out<-try(saveRDS(x,file))
  if (inherits(out,'try-error')) stop('cannot save this pop, check file path')
}

#' Load a population.
#' 
#' Just load a population saved to disk via **[bh_savePopulation()]**.
#' @param file character scalar
#' @export
bh_loadPopulation<-function(file = NULL){
  
  if (is.null(file)) stop('give a full path filename')
  
  out<-try(readRDS(file))
  if (inherits(out,'try-error')) stop('cannot load this pop, check file path')
  return(out)
}
