library(barbieHistologist)

CD45<-bh_defineMarker(markerName = 'cd45',
                      channelName = 'ch1',
                      patternFunction = '.pattern.random',
                      patternModifier = list(mean = 3, sd=3),
                      compartment = 'cytoplasm')
CD4<-bh_defineMarker(markerName = 'cd4',
                     channelName = 'ch2',
                     patternFunction = '.pattern.skin',
                     patternModifier = list(mean = 10, sd=8,trsh=0.3),
                     compartment = 'cytoplasm')
CD8<-bh_defineMarker(markerName = 'cd8',
                     channelName = 'ch3',
                     patternFunction = '.pattern.random',
                     patternModifier = list(mean = 4, sd=2 ),
                     compartment = 'cytoplasm')
DNA<-bh_defineMarker(markerName = 'dna',
                     channelName = 'ch4',
                     patternFunction = '.pattern.random',
                     patternModifier = list(mean = 3, sd=3),
                     compartment = 'nucleus')

cytop1<-bh_defineShape(majorAxis = c(6,0.5),
                   minorAxis = c(6,0.5),
                   roundness = c(0.7,0),
                   nArms = 8,
                   armExt = c(2,0.1),
                   armElbow = 5,
                   armSwing = 1,
                   armTrig = c(-2,0.1))

cytop2<-bh_defineShape(majorAxis = c(4,0.5),
                       minorAxis = c(4,0.5),
                       roundness = c(1,0),
                       nArms = 8,
                       armExt = c(1,0.1),
                       armElbow = 2,
                       armSwing = 1,
                       armTrig = c(-2,0.1))

nuc1<-bh_defineShape(majorAxis = c(1.5,0.001),
                       minorAxis = c(1.5,0.001),
                       roundness = c(1,0),
                       nArms = 8,
                       armExt = c(1,0),
                       armElbow = 2,
                       armSwing = 0,
                       armTrig = c(-2,0.1))

             
SS1<-bh_create(cytop2)
plot(SS1$outline)
lapply(SS1$branch,plot,add=T)
SS2<-bh_create(nuc1)
plot(SS2$outline,add=T)

cell1<-bh_defineCell(cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD4,DNA))

cell2<-bh_defineCell(cytoplasm = cytop2,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD8,DNA))

tissue1<-bh_defineTissue(coords = c(0,100,0,100),
                         resolution = 1,
                         bg = 0,
                         markers = list(CD45,CD4,CD8,DNA))

TEMP_population<-bh_populate(cellPrototype = list(cell1,cell2),
                  cellNames = c('CD4','CD8'),
                  proportion = c(0.5,0.5),
                  tissue = tissue1,
                  areaTresh = 0.8)


TEMP_pic<-takePicture(tissue = tissue1,
                  cells = TEMP_population$CD4)
TEMP_pic<-takePicture(tissue = TEMP_pic,
                  cells = TEMP_population$CD8)

TEMP_pic$x.ch1.ch1<-raster::focal(TEMP_pic$x.ch1.ch1,raster::focalWeight(TEMP_pic$x.ch2.ch2, 0.5, "Gauss"))
raster::plot(TEMP_pic,col=gray.colors(n = 255,0,1))
TEMP1<-bh_asSFC(cells = TEMP_population,
                compartments = list('cytoplasm','nucleus'))

plot(TEMP1,add=T,col=NA,border='red')
plot(TEMP1)
TEMP2<-bh_asXYZ(tissue = TEMP)

colnames(TEMP2)<-c('X','Y','ch1-CD45(ch1)','ch2-CD4(ch2)','ch3-CD8(ch3)','ch4-DNA(ch4)')
write.table(TEMP,"C:/Users/k1343421/Documents/IMC/BarbieHistologist/BH.txt",
            row.names = F,
            quote = F,
            sep='\t')

