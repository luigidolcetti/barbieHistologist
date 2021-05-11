library(barbieHistologist)
cytop1<-bh_defineShape(majorAxis = c(2,1),
                       minorAxis = c(2,1),
                       roundness = c(1,0),
                       nArms = 40,
                       fixedArms = T,
                       armExt = c(1,0.1),
                       armElbow = 2,
                       armSwing = 1,
                       armTrig = c(-0.33,0.1))
cytop2<-bh_defineShape(majorAxis = c(6,1),
                       minorAxis = c(6,1),
                       roundness = c(1,0),
                       nArms = 6,
                       fixedArms = T,
                       armExt = c(1,0.1),
                       armElbow = 2,
                       armSwing = 0,
                       armTrig = c(-0.33,0.1))
nuc1<-bh_defineShape(majorAxis = c(1,0.1),
                     minorAxis = c(1,0.1),
                     roundness = c(1,0),
                     nArms = 40,
                     fixedArms = T,
                     armExt = c(1,0),
                     armElbow = 2,
                     armSwing = 0,
                     armTrig = c(-1,0.1))
nuc2<-bh_defineShape(majorAxis = c(3,0.1),
                     minorAxis = c(3.5,0.1),
                     roundness = c(1,0),
                     orientation = NULL,
                     nArms = 40,
                     fixedArms = T,
                     armExt = c(1,0),
                     armElbow = 2,
                     armSwing = 0,
                     armTrig = c(-1,0.1))
CD45<-bh_defineMarker(markerName = 'cd45',
                      channelName = 'ch1',
                      patternFunction = '.pattern.random',
                      patternModifier = list(mean = 5, sd=0),
                      compartment = 'cytoplasm')
CD46<-bh_defineMarker(markerName = 'CD46',
                      channelName = 'ch1',
                      patternFunction = '.pattern.random',
                      patternModifier = list(mean = 5, sd=0),
                      compartment = 'nucleus')

cell1<-bh_defineCell(name = 'aaa',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     markers = list(CD45,CD46))

cell2<-bh_defineCell(name = 'bbb',
                     cytoplasm = cytop2,
                     nucleus = nuc2,
                     markers = list(CD45,CD46))

TEMP_tissue<-bh_defineTissue(coords = c(0,100,0,100),
                             resolution = c(1,1),
                             bg = 0,
                             markers = list(CD45))

TEMP<-bh_create(cell2,lox=1,loy=1)
TEMP1<-bh_translate(TEMP,c(3,3))
plot(TEMP@cytoplasm$outline)


TEMP8<-bh_populate_byInteract(cellPrototype = list(cell1,cell2),
                              proportion = c(0.5,0.5),
                              tissue = TEMP_tissue,
                              cropToMesure = T,
                              areaTresh = 0.90)

TEMP8<-bh_populate_byClip(cellPrototype = list(cell1,cell2),
                          proportion = c(0.01,0.99),
                          maxCloning = 3,
                          
                          tissue = TEMP_tissue,
                          cropToMesure = F,
                          areaTresh = 0.99)


engr_tissue<-bh_engrave(TEMP_tissue,TEMP8$aaa)
engr_tissue<-bh_engrave(engr_tissue,TEMP8$bbb)
raster::plot(engr_tissue)
GEOM_list<-bh_asSFC(cells = TEMP8)
plot(GEOM_list[3],col=NA,border=c('blue','red','black'),add=T)
plot(GEOM_list[1])
