library(barbieHistologist)

M1<-new('genericMarker',
  marker = 'a',
        channel = 'ch1',
        pattern = call('.pattern.random',mean=1),
        compartment = 'cytoplasm')

M3<-new('genericMarker',
        marker = 'c',
        channel = 'ch3',
        pattern = call('.pattern.skin'),
        compartment = 'cytoplasm')

M2<-new('genericMarker',
  marker = 'b',
        channel = 'ch2',
        pattern = call('.pattern.concentric'),
        compartment = 'nucleus')






TEMP_tissue<-new('tissue',
                 xmn = 0,
                 xmx = 500,
                 ymn = 0,
                 ymx = 500,
                 bg = 0,
                 xRes = 1,
                 yRes = 1,
                 markers = list(M1,M2,M3))

raster::plot(TEMP_tissue)

TEMP_Cy<-new('genericShape',
             majorAxis = c(3,0.1),
             minorAxis = c(4,0.1),
             roundness = c(0.5,0.1),
             nArms = 4,
             armExt = c(3,0.1),
             armElbow = 4,
             armSwing = 1,
             armTrig = c(-5,1))

TEMP_Nc<-new('genericShape',
             majorAxis = c(3,0.1),
             minorAxis = c(3,0.1),
             roundness = c(1,0.1),
             nArms = 4,
             armExt = c(0.2,0.1),
             armElbow = 3,
             armSwing = 1,
             armTrig = c(-1,1))

TEMP<-bh_create(TEMP_Cy)
plot(TEMP$outline)
TEMP<-bh_create(TEMP_Nc)
plot(TEMP$outline,add=T)


TEMP_Cell<-new('cellPrototype',
               cytoplasm = TEMP_Cy,
               nucleus = TEMP_Nc,
               organelle = NULL,
               markers = list(M1,M2,M3))

TEMP_real<-vector(mode = 'list',length = 10)
for (i in 1:10){
TEMP_real[[i]]<-bh_create(x = TEMP_Cell,
      
                     lox = sample(50:450,1),
                     loy = sample(50:450,1))
}

TEMP_rst<-takePicture(TEMP_tissue,TEMP_real)
raster::plot(TEMP_rst,col=gray.colors(255,start=0,end=1))

