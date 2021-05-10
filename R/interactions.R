.interact<-function(x,y){
  yOutline<-y$outline
  xInputStem<-x$stem
  xInputBranch<-x$branch
  xOrg<-x$centroid
  
  xInputStemSep<-sf::st_cast(xInputStem,'LINESTRING')
  xTrimStemSep<-sf::st_difference(xInputStemSep,yOutline)
  newSfc<-vector('list',length(xInputStemSep))
  for (i in 1:length(xTrimStemSep)){
    if (sf::st_geometry_type(xTrimStemSep[[i]])=='MULTILINESTRING'){
      TEMP<-sf::st_cast(sf::st_sfc(xTrimStemSep[[i]]),'LINESTRING')
      TEMP_which<-which(sf::st_intersects(xOrg,TEMP,sparse = F))
      newSfc[[i]]<-TEMP[[TEMP_which]]
    } else {
      newSfc[[i]]<-xTrimStemSep[[i]]
    }
  }
  
  xTrimStemSep<-sf::st_sfc(newSfc)
  orgLenght<-sf::st_length(xInputStemSep)
  newLength<-sf::st_length(xTrimStemSep)
  
  newXstem<-xInputStemSep
  newXbranch<-x$branch
  for (i in 1:length(xInputStemSep)){
    shiftLineString<-xInputStemSep[[i]]

    shiftLineString<-((shiftLineString-xOrg[[1]])*(newLength[i]/orgLenght[i]))+xOrg[[1]]
    
    newXstem[[i]]<-shiftLineString
    newXBranch[[i]]<-((newXBranch[[i]]-xOrg[[1]])*(newLength[i]/orgLenght[i]))+xOrg[[1]]
  }
  
  compactXstem<-sf::st_cast(newXstem,'MULTILINESTRING')

  newOutline<-.envelope(newXBranch)

  out<-list(centroid = xOrg,
            stem = compactXstem,
            branch = newXbranch,
            outline = newOutline)
  return(out)
  }


.enclose<-function(x,y){
  if (length(sf::st_contains_properly(x$outline,y$outline)[[1]])!=0) return(x)
  if (length(sf::st_covered_by(x$centroid,y$outline))==0) return(0)
  
  xInputStem<-x$stem
  xInputBranch<-x$branch
  xOrg<-x$centroid
  yOutline<-sf::st_cast(y$outline,'LINESTRING')
  
  xInputStemSep<-sf::st_cast(xInputStem,'LINESTRING')
  
  
  xInputStempIntersectMatrix<-sf::st_intersects(xInputStemSep,yOutline)
  
  newSfc<-vector('list',length(xInputStemSep))
  
  capLength<-lapply(xInputStemSep,function(xISS){
    oldLength<-sf::st_length(xISS)
    xInputStempIntersect<-sf::st_intersection(xISS,yOutline)
    if (length(xInputStempIntersect)==0) newLength<-oldLength else {
      
      xInputStempIntersect<-sf::st_cast(sf::st_sfc(xInputStempIntersect),'POINT')
      dst<-sf::st_distance(x$centroid,xInputStempIntersect)
      newLength<-min(dst)
    }
    return(list(old=oldLength,new=newLength))
    })
  
  newXbranch<-lapply(1:length(xInputBranch),function(i){
    ((xInputBranch[[i]]-xOrg[[1]])*(capLength[[i]]$new/capLength[[i]]$old))+xOrg[[1]]
  })
  
  newXbranch<-sf::st_sfc(newXbranch)
  
  newXstem<-lapply(1:length(xInputStemSep),function(i){
    ((xInputStemSep[[i]]-xOrg[[1]])*(capLength[[i]]$new/capLength[[i]]$old))+xOrg[[1]]
  })
  
  compactXstem<-sf::st_sfc(newXstem)
  compactXstem<-sf::st_cast(compactXstem,'MULTILINESTRING')
  
  
  newOutline<-.envelope(newXbranch)
  
  out<-list(centroid = xOrg,
            stem = compactXstem,
            branch = newXbranch,
            outline = newOutline)
  return(out)
  
}
