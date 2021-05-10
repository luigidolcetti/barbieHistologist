.envelope<-function(x){
  newXBranchList<-sf::st_coordinates(x)
  newXBranchList<-split(newXBranchList[,1:2],newXBranchList[,4])
  newXBranchList<-lapply(newXBranchList,matrix,ncol=2,byrow=F)
  # centroid<-sf::st_coordinates(centroid)
  
  outLine<-lapply(newXBranchList,function(x){
    lx<-x[c(T,F,F),,drop=F]
    rx<-x[c(F,F,T),,drop=F]
    rx<-rx[nrow(rx):1,]
    out<-rbind(lx,rx)
    })
  outLine<-do.call(rbind,outLine)
  outLine<-rbind(outLine,outLine[1,])
  outLine<-sf::st_polygon(list(outLine))
  outLine<-sf::st_buffer(outLine,0)
  if (!sf::st_is_valid(outLine)){
    outLine<-sf::st_make_valid(outLine)}
  return(outLine)
}
