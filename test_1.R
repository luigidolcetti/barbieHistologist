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
                     patternFunction = '.pattern.skin',
                     patternModifier = list(mean = 4, sd=2 ,trsh=0.2),
                     compartment = 'cytoplasm')
DNA<-bh_defineMarker(markerName = 'dna',
                     channelName = 'ch4',
                     patternFunction = '.pattern.random',
                     patternModifier = list(mean = 3, sd=3),
                     compartment = 'nucleus')

cytop1<-bh_defineShape(majorAxis = c(4,1),
                   minorAxis = c(4,1),
                   roundness = c(0.2,0),
                   nArms = 8,
                   armExt = c(4,0.2),
                   armElbow = 5,
                   armSwing = 1,
                   armTrig = c(-2,0.1))

cytop2<-bh_defineShape(majorAxis = c(4,1),
                       minorAxis = c(4,1),
                       roundness = c(1,0),
                       nArms = 12,
                       armExt = c(1.1,0.1),
                       armElbow = 2,
                       armSwing = 1,
                       armTrig = c(-2,0.1))

nuc1<-bh_defineShape(majorAxis = c(3,0.001),
                       minorAxis = c(3,0.001),
                       roundness = c(1,0),
                       nArms = 8,
                       armExt = c(1,0),
                       armElbow = 2,
                       armSwing = 0,
                       armTrig = c(-2,0.1))

             
SS1<-bh_create(cytop2)
plot(SS1$outline)
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

tissue1<-bh_defineTissue(coords = c(0,200,0,200),
                         resolution = 1,
                         bg = 0,
                         markers = list(CD45,CD4,CD8,DNA))

real_CD4<-vector(mode = 'list',length = 5)
for (i in 1:5){
real_CD4[[i]]<-bh_create(cell1,lox=sample(1:200,1),loy=sample(1:200,1))
}

real_CD8<-vector(mode = 'list',length = 10)
for (i in 1:10){
  real_CD8[[i]]<-bh_create(cell2,lox=sample(1:200,1),loy=sample(1:200,1))
}

TEMP<-takePicture(tissue = tissue1,
                  cells = real_CD4)
TEMP<-takePicture(tissue = TEMP,
                  cells = real_CD8)


raster::plot(TEMP$x.ch1.ch1,col=gray.colors(n = 255,0,1))
TEMP1<-bh_asSFC(cells = list(CD4=real_CD4,
                             CD8=real_CD8),
                compartments = list('cytoplasm','nucleus'))
plot(TEMP1,add=T,col=NA,border='red')

TEMP2<-bh_asXYZ(tissue = TEMP)

colnames(TEMP2)<-c('X','Y','ch1-CD45(ch1)','ch2-CD4(ch2)','ch3-CD8(ch3)','ch4-DNA(ch4)')
write.table(TEMP,"C:/Users/k1343421/Documents/IMC/BarbieHistologist/BH.txt",
            row.names = F,
            quote = F,
            sep='\t')

