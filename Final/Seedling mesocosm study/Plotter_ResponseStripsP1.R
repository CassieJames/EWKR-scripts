yStrings <- c('nLeaves','Coppice','LeafArea','Height','RootLength','MassAbove','MassBelow','MassAB')
yLabels <- c('Leaf count', 'Coppicing branch count', expression('Leaf area (cm'^2*')'), 'Height (cm)',
             'Root length (cm)', 'Above-ground mass (g)', 'Below-ground mass (g)', 'Above:Below mass (g)')
xLabels <- c('FD','DD','FF')

tiff('Fig_TreatmentMeansP1.tiff', width = 16, height = 12, units = 'cm', res = 300)

myPar <- par(mfrow = c(2,4), mar = c(3,5,2,1), oma = c(3,1,1,1))
pCount <- 1

for (Y in yStrings) {
  Ymean <- SLstats[[paste(Y,'Mean',sep='')]]
  YSE <- SLstats[[paste(Y,'SE',sep='')]]
  X <- as.numeric(SLstats$Treat_newP1) # Numeric to stop automatic boxplot when parsed factor
  
  spp <- unique(SLstats$Species)
  myPch <- c(21,22,25)
  myBg <- c('black','green','red')#c(grey(0.8))
  pointCex <- 1.6
  
  yMax <- max(Ymean+2*YSE)
  plot(X, y = Ymean, axes = FALSE, xlab = '', ylab = '', type = 'n', ylim = c(0,yMax), xlim = c(0.5,3.5))
  # Lines:
  #for (s in spp) {
  #  yL <- Ymean[SLstats$Species == s]
  #  xL <- X[SLstats$Species == s]
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
  
  if (pCount == 6) {mtext('Watering treatment', side = 1, line = 4, at = 8, cex = 0.8)}
  
  pCount <- pCount+1
}
par(myPar)
dev.off()