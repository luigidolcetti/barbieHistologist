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
                     patternModifier = list(mean = 4, sd=4,k=1,X0=0.5 ),
                     compartment = 'cytoplasm')
CD11c<-bh_defineMarker(markerName = 'cd11c',
                       channelName = 'ch4',
                       patternFunction = '.pattern.random',
                       patternModifier = list(mean = 4, sd=100),
                       compartment = 'cytoplasm')
CD1Low<-bh_defineMarker(markerName = 'cd1000',
                        channelName = 'ch5',
                        patternFunction = '.pattern.random',
                        patternModifier = list(mean = 2, sd=3),
                        compartment = 'cytoplasm')
CD1high<-bh_defineMarker(markerName = 'cd1000',
                         channelName = 'ch5',
                         patternFunction = '.pattern.random',
                         patternModifier = list(mean = 8, sd=4),
                         compartment = 'cytoplasm')
DNA<-bh_defineMarker(markerName = 'dna',
                     channelName = 'ch6',
                     patternFunction = '.pattern.random',
                     patternModifier = list(mean = 10, sd=5),
                     compartment = 'nucleus')









cytop1<-bh_defineShape(majorAxis = c(3,0.5),
                       minorAxis = c(3,0.5),
                       roundness = c(1,0),
                       nArms = 8,
                       armExt = c(1,0.1),
                       armElbow = 2,
                       armSwing = 0,
                       armTrig = c(-2,0.1))

cytop2<-bh_defineShape(majorAxis = c(2,0.5),
                       minorAxis = c(2,0.5),
                       roundness = c(1,0),
                       nArms = 8,
                       armExt = c(3,0.1),
                       armElbow = 4,
                       armSwing = 2,
                       armTrig = c(-2,0.1))

nuc1<-bh_defineShape(majorAxis = c(1.5,0.1),
                     minorAxis = c(1.5,0.1),
                     roundness = c(1,0),
                     nArms = 8,
                     armExt = c(1,0),
                     armElbow = 2,
                     armSwing = 0,
                     armTrig = c(-1,0.1))


cell1<-bh_defineCell(name = 'CD4_lo',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD4,DNA,CD1Low))
cell2<-bh_defineCell(name = 'CD4_hi',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD4,DNA,CD1high))

cell3<-bh_defineCell(name = 'CD8_lo',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD8,DNA,CD1Low))
cell4<-bh_defineCell(name = 'CD8_hi',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD8,DNA,CD1high))

cell5<-bh_defineCell(name = 'dendr',
                     cytoplasm = cytop2,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,DNA,CD11c))


  tissue1<-bh_defineTissue(coords = c(0,50,0,50),
                           resolution = 1,
                           bg = 0,
                           markers = list(CD45,CD4,CD8,CD11c,CD1high,DNA))
  
  TEMP_population<-bh_populate(cellPrototype = list(cell1,cell2,cell3,cell4,cell5),
                               proportion = c(0.2,0.2,0.2,0.2,0.2),
                               tissue = tissue1,
                               maxCloning = 3,
                               areaTresh = 0.5)
  GEOM_list<-bh_asSFC(cells = TEMP_population)
  
  TEMP<-bh_familyPicture(tissue = tissue1,
                         sfc = GEOM_list)
 