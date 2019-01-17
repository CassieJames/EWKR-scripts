#### Model comparisons


png("Candidate models comparison.png",width=35, height=45, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1),mfrow=c(3,2))
		
		boxplot(cv ~  model, data = tmp.abund.wet, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-likelihood", xlab = "")
		axis(1, at = 1:6, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.abund.wet$model)), par("usr")[3] - 0.09,srt = 30, adj = 1, label = levels(tmp.abund.wet$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.abund.wet), tmp.abund.wet$b, function(x) lines(1:6, tmp.abund.wet[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

		boxplot(cv ~  model, data = tmp.rich, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-likelihood", xlab = "")
		axis(1, at = 1:6, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.rich$model)), par("usr")[3] - 0.015,srt = 30, adj = 1, label = levels(tmp.rich$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.rich), tmp.rich$b, function(x) lines(1:6, tmp.rich[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")
		
		boxplot(cv ~  model, data = tmp.ter, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-likelihood", xlab = "")
		axis(1, at = 1:6, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.ter$model)), par("usr")[3] - 0.05,srt = 30, adj = 1, label = levels(tmp.ter$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.ter), tmp.ter$b, function(x) lines(1:6, tmp.ter[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")
		
		boxplot(cv ~  model, data = tmp.rich.ter, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-likelihood", xlab = "")
		axis(1, at = 1:6, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.rich.ter$model)), par("usr")[3] -0.008,srt = 30, adj = 1, label = levels(tmp.rich.ter$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.rich.ter), tmp.rich.ter$b, function(x) lines(1:6, tmp.rich.ter[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")
		
		boxplot(cv ~  model, data = tmp.ExoticAbund, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:6, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.ExoticAbund$model)), par("usr")[3] - 0.09,srt = 30, adj = 1, label = levels(tmp.ExoticAbund$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.ExoticAbund), tmp.ExoticAbund$b, function(x) lines(1:6, tmp.ExoticAbund[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")
		
		
dev.off()