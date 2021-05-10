library(barbieHistologist)
cytop2<-bh_defineShape(majorAxis = c(6,1),
                       minorAxis = c(6,1),
                       roundness = c(0.5,0),
                       nArms = 40,
                       fixedArms = T,
                       armExt = c(3,0.1),
                       armElbow = 6,
                       armSwing = 1,
                       armTrig = c(-0.33,0.1))
nuc1<-bh_defineShape(majorAxis = c(5,0.1),
                     minorAxis = c(1,0.1),
                     roundness = c(0.3,0),
                     nArms = 40,
                     fixedArms = T,
                     armExt = c(1,0),
                     armElbow = 2,
                     armSwing = 0,
                     armTrig = c(-1,0.1))

CD45<-bh_defineMarker(markerName = 'cd45',
                      channelName = 'ch1',
                      patternFunction = '.pattern.random',
                      patternModifier = list(mean = 5, sd=5),
                      compartment = 'cytoplasm')

TEMP<-bh_create(cytop2)
TEMP1<-bh_create(nuc1)

TEMP1<-barbieHistologist:::bh_translate(TEMP1,c(7,0))

TEMP2<-.enclose(TEMP1,TEMP)


plot(TEMP$outline)
plot(TEMP1$outline,add=T)
plot(TEMP1$stem,add=T)
plot(TEMP1$centroid,add=T)

TEMP2<-.enclose(TEMP1,TEMP)
TEMP1<-barbieHistologist:::bh_translate(TEMP,c(5,5))

plot(TEMP$outline)
plot(TEMP1$outline,add=T)

TEMP2<-.interact(TEMP1,TEMP)

plot(NA,xlim=c(-10,10),ylim=c(-10,10))
plot(TEMP$outline,add=T)
plot(TEMP2$outline,add=T)
plot(TEMP1$outline,add=T)
plot(TEMP$branch)
plot(TEMP$stem,add=T)
plot(TEMP$outline,add=T)
plot(TEMP2$branch,add=F)
plot(TEMP2$stem,add=T)
plot(TEMP2$outline,add=T)
plot(TEMP1$stem,add=T)
plot(TEMP1$branch,add=T)
plot(TEMP1$outline,add=T)

plot(TEMP1$outline)
plot(TEMP2$outline,add=T,col='red')

TEMP3<-sf::st_union(TEMP$outline,TEMP1$outline)
TEMP3<-sf::st_buffer(TEMP3,0,0)
plot(TEMP3)

TEMP5<-bh_defineCell(name = 'aaa',
                     cytoplasm = cytop2,
                     nucleus = nuc1,
                     markers = list(CD45))

TEMP6<-bh_create(TEMP5,lox=1,loy=1)
plot(TEMP6@cytoplasm$outline)
plot(TEMP6@nucleus$outline,add=T)
