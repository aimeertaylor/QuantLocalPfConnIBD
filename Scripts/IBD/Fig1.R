############################################################################################
# 1) Map based on fraction of IBD proportions > 0.5 (2 plots, combined in manuscript)
# This map was created using QGIS (QGIS Development Team 2017) and the raster (Hijmans 2016),
# maptools (Bivand and Lewin-Koh 2017) and igraph (Csardi and Nepusz 2006) packages
############################################################################################
rm(list = ls())

library(maptools)
library(raster)
library(igraph)

# The compass rose was adapted from an R-sig-Geo post by Jim Lemon in 2010
# source: https://stat.ethz.ch/pipermail/r-sig-geo/2010-February/007564.html
compassRose<-function(x,y,rot=0,cex=1) {
  oldcex<-par(cex=cex)
  mheight<-strheight("M")
  xylim<-par("usr")
  plotdim<-par("pin")
  xmult<-(xylim[2]-xylim[1])/(xylim[4]-xylim[3])*plotdim[2]/plotdim[1]
  point.angles<-seq(0,7*pi/4,by=pi/4)+pi*rot/180
  crspans<-rep(c(mheight*3,mheight/2),4)
  xpoints<-cos(point.angles)*crspans*xmult+x
  ypoints<-sin(point.angles)*crspans+y
  polygon(xpoints,ypoints)
  txtxpoints<-cos(point.angles[c(1,3,5,7)])*1.33*crspans[1]*xmult+x
  txtypoints<-sin(point.angles[c(1,3,5,7)])*1.33*crspans[1]+y
  text(txtxpoints,txtypoints,c("E","N","W","S"), cex = 1.2)
  par(oldcex)
}

# Load shape files (downloaded from http://gadm.org/)
thailand <- getData("GADM",country="Thailand",level=0)
myanmar <- getData("GADM",country="Myanmar",level=0)
tm_border <- readShapeLines("/Users/aimeet/Documents/BroadLaptop/TM_border/QGIS/TM_border.shp")

# Load data
load('../../RData/Barcode_threshold.RData')

# Vertices
site_comps <- unique(Barcode$Site_comparison)

# Store in which to store edge weights
mylinks <- data.frame(array(dim = c(length(site_comps), 3),
                            dimnames = list(site_comps,c('site1', 'site2','weight_barcode'))))
mylinks[,c('site1','site2')] <- do.call(rbind, strsplit(site_comps, split = '_'))

# Extract edge weights and directionality
for(site_comp in site_comps){
  ind_barcode <- Barcode$Site_comparison == site_comp
  p_barIBD <- mean(Barcode$ProbIBD_93_tail[ind_barcode], na.rm = TRUE)
  mylinks[site_comp, 'weight_barcode'] <- p_barIBD
}

# Generate net, remove loops
net <- graph.data.frame(mylinks[-c(1,5,8,10),], directed = FALSE)

# Coordinates
lonlat <- read.csv('../../RawData/SMRU_clinics_lonlat.csv')
write.table(lonlat, file = '../../RawData/SMRU_clinics_lonlat.txt',
            quote = FALSE, row.names = FALSE)
rownames(lonlat) <- lonlat[,'clinic']
layout_lonlat <- as.matrix(lonlat[rownames(as.matrix(V(net))), c('lon','lat')])

# IBD barcode
par(mfrow = c(1,1))
E(net)$width <- E(net)$weight_barcode * 5000
#E(net)$curved <- 0.69*c(-1, 0.8, -1.2, 1, -1, 0.6, 1.2, 1, -1.2,  1)
E(net)$curved <- 0.59*c(0.8, -1.3, 1, 0.6, 1.3, -1.3)
E(net)$loop.angle <- 1.68*pi
E(net)$color <- adjustcolor('black', alpha.f = 0.75)

# Change labels, text colour etc.
print(V(net)) # Print order
V(net)$label <- c('Mawker Thai', 'Mae Kon Ken', 'Maela', 'Wang Pha')
V(net)$label.color <- "black"
V(net)$label.dist <- c(1.4, 2, 1.2, 1.4)
V(net)$label.degree <- c(pi, 2*pi, 1.5*pi, 2*pi)
V(net)$label.cex <- 1.2
V(net)$color <- "black"
V(net)$size <- 2 # Change to clinic size? Or sample size? Proves prblematic!
V(net)$frame.color <- "white"


# Plot context
png(file = './Network.png', height = 7, width = 6, units = 'in', res = 400)
par(mfrow = c(1,1), family = "serif", mar = c(0,0,0,0))
# Plot border with net
plot(tm_border, col = 'darkblue', lwd = 2)
plot(net, layout=layout_lonlat, add = TRUE, rescale = FALSE)
compassRose(x = 98.25, y = 17)
# Scalebar below doesn't work. Used to calibrate segment length
#scalebar(d = 10, xy = c(98.202, 16.45), type = "line", lonlat = TRUE, lwd = 2, cex = 1.2)
text(x = 98.25, y = 16.492, labels = '10 km', cex = 1.2, pos = 3)
segments(y0 = 16.5, y1 = 16.5, x0 = 98.202, x1 = 98.296, lwd = 2, col = 'darkblue')
dev.off()

pdf(file = './Map.pdf', height = 7, width = 6)
par(mfrow = c(1,1), family = "serif", mar = c(0,0,0,0))
plot(thailand, col = 'gray', xlim = c(95, 102), ylim = c(5,28), border = 'white')
plot(myanmar, add = TRUE, col = 'gray', border = 'white')
text(x = c(102,96), y = c(16,22), labels = c('Thailand','Myanmar'), col = 'white', cex = 2)
plot(tm_border, col = 'darkblue', lwd = 2, add = TRUE)
dev.off()

