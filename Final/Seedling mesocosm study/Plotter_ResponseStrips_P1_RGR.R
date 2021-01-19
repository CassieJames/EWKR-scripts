
yStrings <- c('Height','MassAbove','MassBelow','MassTotal','MassAB','nLeaves','LeafArea','RootLength','Coppice')
			 
yLabels <- c('Height (cm)', 'Above-ground mass (g)', 'Below-ground mass (g)','Total mass (g)', 'Above:Below mass (g)','Leaf count',expression('Total leaf area (cm'^2*')'),'Root length (cm)','Coppicing branch count')
			 
xLabels <- c('DD','FD','FF')

tiff('Fig_TreatmentMeansRGR_P1.tiff', width = 16, height = 16, units = 'cm', res = 300)

myPar <- par(mfrow = c(3,3), mar = c(3,5,2,1), oma = c(3,1,1,1))
pCount <- 1

for (Y in yStrings) {
  Ymean <- SLstatsP1RGR[[paste(Y,'Mean',sep='')]]
  YSE <- SLstatsP1RGR[[paste(Y,'SE',sep='')]]
  X <- as.numeric(SLstatsP1RGR$Treat_newP1) # Numeric to stop automatic boxplot when parsed factor
  
  spp <- unique(SLstatsP1RGR$Species)
  myPch <- c(21,22,25)
  myBg <- c('black','green','red')#c(grey(0.8))
  pointCex <- 1.6
  
  yMax <- max(Ymean+2*YSE)
  yMin <- min(Ymean-2*YSE)
  plot(X, y = Ymean, axes = FALSE, xlab = '', ylab = '', type = 'n', ylim = c(yMin,yMax), xlim = c(0.5,3.5))
  # Lines:
  #for (s in spp) {
  #  yL <- Ymean[SLstatsP1RGR$Species == s]
  #  xL <- X[SLstatsP1RGR$Species == s]
  #  lines(xL, yL, col = grey(0.7))
  #}
  # Points:
  arrows(X, Ymean, X, Ymean+1.96*YSE, length = 0, angle = 90)
  arrows(X, Ymean, X, Ymean-1.96*YSE, length = 0, angle = 90)
  points(X, Ymean, pch = myPch, cex = pointCex, bg = myBg)
  axis(1, labels = xLabels, at = 1:3, las = 3)
  axis(2)
  l = letters[pCount]
  title(main = paste('(',l,')', sep = ''), font.main = 1, adj = 0)
  title(ylab = yLabels[pCount])
  box(bty = 'l')
  
  if (pCount == 8) {mtext('Watering treatment', side = 1, line = 4, at = 2, cex = 0.8)}
  
  pCount <- pCount+1
}
par(myPar)
dev.off()