#' @export
bh_populate<-function(cellPrototype = NULL,
                      cellNames = NULL,
                      proportion = NULL,
                      tissue = NULL,
                      areaTresh = NULL){
  
  blob<-sf::st_polygon(list(matrix(c(0,0,1,0,1,1,0,1,0,0),ncol=2,byrow = T)))
  areaTot<-0
  lim <-raster::extent(tissue)
  res <-raster::res(tissue)
  possibleXY<-raster::raster(vals=0,ext=lim,resolution=res)
  areaTresh<-(lim[2]-lim[1])*(lim[4]-lim[3])*areaTresh
  # cellNames<-lapply(cellPrototype,function(x) deparse(substitute(x)))
  newCellList<-sapply(1:length(cellNames),function(clnms){
    vector('list',round(areaTresh*proportion[clnms]))},simplify = F,USE.NAMES = F)
  names(newCellList)<-cellNames
  newCellType_counter<-rep(1,length(cellNames))
  names(newCellType_counter)<-cellNames
  newBar<-utils::txtProgressBar(min=0,max=areaTresh,initial=0,style=3)
  while(areaTot<areaTresh){
    
    newCellType<-sample(1:length(cellPrototype),1,prob = proportion)
    position<-raster::Which(possibleXY==0,cells=T)
    position<-sample(position,1)
    position<-raster::xyFromCell(object = tissue,cell = position)
    newCell<-bh_create(cellPrototype[[newCellType]],
                       lox = position[,'x'],
                       loy = position[,'y'])
    newComp<-try(sapply(c('cytoplasm',
                          'nucleus',
                          'organelle'),
                        function(cmp){
                          
                          newComp<-slot(newCell,cmp)$outline
                          if (is.null(newComp)) return(newComp) else {
                            if (sf::st_geometry_type(newComp)=='POLYGON'){
                              newComp<-sf::st_difference(newComp,blob) 
                              if (!sf::st_geometry_type(newComp)=='POLYGON') {
                                
                                
                                newComp<-sf::st_sfc(newComp)
                                newComp<-sf::st_cast(newComp,'POLYGON')
                                newComp_list<-vector(mode = 'list',length = length(newComp))
                                for (i in 1:length(newComp)){
                                  newComp_list[[i]]<-newComp[[i]]}
                                area_list<-sapply(newComp_list,sf::st_area,simplify = T)
                                newComp<-newComp_list[[which.max(area_list)]]
                                
                              }
                            }}
                          return(newComp)
                        },USE.NAMES = T,simplify = F))
    
    if (!inherits(newComp,'try-error')){
      for (cmp in c('cytoplasm',
                    'nucleus',
                    'organelle')){
        slot(newCell,cmp)$outline<-newComp[[cmp]]}
      
      whichTemp<-sapply(newComp,is.null,simplify = T,USE.NAMES = F)
      
      if (length(which(!whichTemp))>1) {
      unionCell<-sf::st_union(sf::st_sfc(newComp[!whichTemp]))[[1]]} else {
        unionCell<-newComp[!whichTemp]
      }
      
      if (is.null(blob)) {blob<-unionCell} else {
        
        dummyBlob<-blob
        dummyBlob<-try(sf::st_union(dummyBlob,unionCell))
        
        if (!inherits(dummyBlob,'try-error')){
          blob<-dummyBlob
          
          areaTot<-sf::st_area(blob)
          
          xyCoverage<-exactextractr::exact_extract(possibleXY,sf::st_sfc(blob),include_xy=T)
          xyCoverage<-raster::cellFromXY(possibleXY,xyCoverage[[1]][,c('x','y')])
          possibleXY[xyCoverage]<-1
          
          newCellList[[newCellType]][[newCellType_counter[newCellType]]]<-newCell
          newCellType_counter[newCellType]<-newCellType_counter[newCellType]+1
          
          setTxtProgressBar(newBar,areaTot)
          
          
        }
      }
    }
  }
  close(newBar)
  whichToGetRid<-lapply(newCellList,function(x){
    sapply(x,is.null,simplify = T)})
  out_newCellList<-lapply(1:length(newCellList),function(x){
    newCellList[[x]][!whichToGetRid[[x]]]})
  names(out_newCellList)<-names(newCellList)
  return(out_newCellList)
}
