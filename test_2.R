library(barbieHistologist)

CD45<-bh_defineMarker(markerName = 'cd45',
                      channelName = 'ch1',
                      patternFunction = '.pattern.random',
                      patternModifier = list(mean = 5, sd=5),
                      compartment = 'cytoplasm')
CD4<-bh_defineMarker(markerName = 'cd4',
                     channelName = 'ch2',
                     patternFunction = '.pattern.skin2',
                     patternModifier = list(mean = 2, sd=3,k=5,X0=0.2),
                     compartment = 'cytoplasm')
CD8<-bh_defineMarker(markerName = 'cd8',
                     channelName = 'ch3',
                     patternFunction = '.pattern.skin2',
                     patternModifier = list(mean = 4, sd=4,k=3,X0=0.3 ),
                     compartment = 'cytoplasm')
CD11c<-bh_defineMarker(markerName = 'cd11c',
                       channelName = 'ch4',
                       patternFunction = '.pattern.skin2',
                       patternModifier = list(mean = 4, sd=2,k=3,X0=0.7 ),
                       compartment = 'cytoplasm')
CD1000Low<-bh_defineMarker(markerName = 'cd1000',
                           channelName = 'ch5',
                           patternFunction = '.pattern.random',
                           patternModifier = list(mean = 2, sd=3),
                           compartment = 'cytoplasm')
CD1000high<-bh_defineMarker(markerName = 'cd1000',
                            channelName = 'ch5',
                            patternFunction = '.pattern.random',
                            patternModifier = list(mean = 8, sd=4),
                            compartment = 'cytoplasm')
DNA<-bh_defineMarker(markerName = 'dna',
                     channelName = 'ch6',
                     patternFunction = '.pattern.random',
                     patternModifier = list(mean = 10, sd=5),
                     compartment = 'nucleus')









cytop1<-bh_defineShape(majorAxis = c(6,0.5),
                       minorAxis = c(6,0.5),
                       roundness = c(1,0),
                       nArms = 8,
                       armExt = c(1.2,0.1),
                       armElbow = 2,
                       armSwing = 0,
                       armTrig = c(-2,0.1))

cytop2<-bh_defineShape(majorAxis = c(4,0.5),
                       minorAxis = c(4,0.5),
                       roundness = c(1,0),
                       nArms = 8,
                       armExt = c(3,0.1),
                       armElbow = 4,
                       armSwing = 2,
                       armTrig = c(-2,0.1))

nuc1<-bh_defineShape(majorAxis = c(3,0.001),
                     minorAxis = c(3,0.001),
                     roundness = c(1,0),
                     nArms = 8,
                     armExt = c(1,0),
                     armElbow = 2,
                     armSwing = 0,
                     armTrig = c(-1,0.1))


# SS1<-bh_create(cytop2)
# plot(SS1$outline)
# lapply(SS1$branch,plot,add=T)
# SS2<-bh_create(nuc1)
# plot(SS2$outline,add=T)

cell1<-bh_defineCell(name = 'CD4_lo',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD4,DNA,CD1000Low))
cell2<-bh_defineCell(name = 'CD4_hi',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD4,DNA,CD1000high))

cell3<-bh_defineCell(name = 'CD8_lo',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD8,DNA,CD1000Low))
cell4<-bh_defineCell(name = 'CD8_hi',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD8,DNA,CD1000high))

cell5<-bh_defineCell(name = 'dendr',
                     cytoplasm = cytop2,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,DNA,CD11c))


tissue1<-bh_defineTissue(coords = c(0,250,0,250),
                         resolution = 1,
                         bg = 0,
                         markers = list(CD45,CD4,CD8,CD11c,CD1000high,DNA))

TEMP_population<-bh_populate(cellPrototype = list(cell1,cell2,cell3,cell4,cell5),
                             proportion = c(0.3,0.15,0.2,0.2,0.1),
                             tissue = tissue1,
                             maxCloning = 15,
                             areaTresh = 0.99)


TEMP_pic<-bh_engrave(tissue = tissue1,
                     cells = TEMP_population$CD4_lo)
TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                     cells = TEMP_population$CD4_hi)
TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                     cells = TEMP_population$CD8_lo)
TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                     cells = TEMP_population$CD8_hi)
TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                     cells = TEMP_population$dendr)

par(mar=c(0,0,0,0))
raster::plot(TEMP_pic,col=gray.colors(n = 255,0,1))
plot(TEMP1['cell'],col=NA,border='red',add=T)

TEMP_mod<-TEMP_pic

raster::plot(TEMP_mod$x.ch1.ch1,col=gray.colors(n = 255,0,1),colNA='blue')

TEMP_mod$x.ch1.ch1<-bh_modifier(TEMP_mod$x.ch1.ch1,wMatrix = matrix(1,5,5),fun=.modifier.multDiv,fun.param = list(quantity=5))
TEMP_mod$x.ch1.ch1<-bh_modifier(TEMP_mod$x.ch1.ch1,wMatrix = matrix(1,5,5),fun = mean)

raster::plot(TEMP_mod$x.ch1.ch1,col=gray.colors(n = 255,0,1),colNA='blue')

raster::plot(TEMP_mod$x.ch2.ch2,col=gray.colors(n = 255,0,1))

TEMP_mod$x.ch2.ch2<-bh_modifier(TEMP_mod$x.ch2.ch2,wMatrix = matrix(1,5,5),fun = mean)
TEMP_mod$x.ch2.ch2<-bh_modifier(TEMP_mod$x.ch2.ch2,wMatrix = matrix(1,5,5),fun=.modifier.multDiv,fun.param = list(quantity=3))

raster::plot(TEMP_mod$x.ch2.ch2,col=gray.colors(n = 255,0,1))

raster::plot(TEMP_mod$x.ch3.ch3,col=gray.colors(n = 255,0,1))

TEMP_mod$x.ch3.ch3<-bh_modifier(TEMP_mod$x.ch3.ch3,wMatrix = matrix(c(0,1,0,1,0,1,0,1,0),3,3),fun = mean)

raster::plot(TEMP_mod$x.ch3.ch3,col=gray.colors(n = 255,0,1))

raster::plot(TEMP_mod$x.ch6.ch6,col=gray.colors(n = 255,0,1))
TEMP_mod$x.ch6.ch6<-bh_modifier(TEMP_mod$x.ch6.ch6,wMatrix = .matrix.gauss(TEMP_mod$x.ch2.ch2,1),fun = mean)
raster::plot(TEMP_mod$x.ch6.ch6,col=gray.colors(n = 255,0,1),colNA='blue')

raster::plot(TEMP_mod$x.ch4.ch4,col=gray.colors(n = 255,0,1))
TEMP_mod$x.ch4.ch4<-bh_modifier(TEMP_mod$x.ch4.ch4,wMatrix = .matrix.gauss(TEMP_mod$x.ch2.ch2,1),fun = median)
raster::plot(TEMP_mod$x.ch4.ch4,col=gray.colors(n = 255,0,1))

raster::plot(TEMP_mod$x.ch5.ch5,col=gray.colors(n = 255,0,1))
TEMP_mod$x.ch5.ch5<-bh_modifier(TEMP_mod$x.ch5.ch5,wMatrix = matrix(1,11,11),fun = .modifier.multDiv,fun.param = list(quantity=3))
TEMP_mod$x.ch5.ch5<-bh_modifier(TEMP_mod$x.ch5.ch5,wMatrix = matrix(1,5,5),fun = median)
raster::plot(TEMP_mod$x.ch5.ch5,col=gray.colors(n = 255,0,1))

raster::plot(TEMP_mod,col=gray.colors(n = 255,0,1))

for (i in names(TEMP_mod)){
  raster::writeRaster(TEMP_mod[[i]],
                      filename = file.path("C:/Users/k1343421/Documents/BH",paste0(i,'_training.tiff')),overwrite=T)
}



TEMP2<-bh_asXYZ(tissue = TEMP_mod)

colnames(TEMP2)<-c('X','Y','ch1-CD45(ch1)','ch2-CD4(ch2)','ch3-CD8(ch3)','ch4-CD11c(ch4)','ch5-CDx(ch5)','ch6-DNA(ch6)')
write.table(TEMP2,"C:/Users/k1343421/Documents/BH/BH_TEST.txt",
            row.names = F,
            quote = F,
            sep='\t')

TEMP1<-bh_asSFC(cells = TEMP_population)

plot(TEMP1['cell'],col=NA)
sf::st_write(TEMP1,"C:/Users/k1343421/Documents/BH/ground_TEST.sqlite",)
