#' @export
methods::setGeneric(name = 'bh_translate',
                    signature = 'x',
                    function(x,...){})
#' @export
methods::setMethod(f = 'bh_translate',
                   signature = 'list',
                   definition = function(x,targetCenter){

                     cntr<-x$centroid
                     out<-lapply(x,function(x){
                       # cntr<-sf::st_centroid(x)
                       out<-x-cntr+targetCenter
                     })
                     
                     out<-new('simpleShape',
                              centroid = out$centroid,
                              stem = out$stem,
                              branch = out$branch,
                              outline = out$outline)
                     
                     return(out)
                   })

#' @export
methods::setMethod(f = 'bh_translate',
                   signature = 'simpleShape',
                   definition = function(x,targetCenter){
                     
                     cntr<-x$centroid
                     out<-lapply(x,function(x){
                       # cntr<-sf::st_centroid(x)
                       out<-x-cntr+targetCenter
                     })
                     
                     out<-new('simpleShape',
                              centroid = out$centroid,
                              stem = out$stem,
                              branch = out$branch,
                              outline = out$outline)
                     
                     return(out)
                   })

.rot = function(a) matrix(c(cos(a), sin(a)), ncol = 2, nrow = 1)
