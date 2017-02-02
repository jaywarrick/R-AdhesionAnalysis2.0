rm(list=ls())
source('D:/GitHub/R-General/.Rprofile')
source('D:/GitHub/R-Adhesion2/ParticleTracking/R/PackageFunctions.R')
source('D:/GitHub/R-CPD/cpd/R/cpd.R')
library(parallel)

# library(foreign)
# fileTable <- read.arff('G:/Jex/Brightfield Adhesion II/temp/JEXData0000000121.arff')
# duh <- data.table(read.arff(fileTable$Value[1]))
# duh <- reorganize(duh)
# convertColsToNumeric(duh, exceptions=c('Surface'))
# duh <- duh[order(Time)]
# lapply(duh, mode)
# duh <- duh[Point==716]
# View(duh)
# plot(duh$Time, duh$Mean, type='l')

dt <- getData(dbPath="G:/Jex/Brightfield Adhesion II", ds="20170127", x=0, y=0, type='File', name='Compiled Table')
dt2 <- getData(dbPath="G:/Jex/Brightfield Adhesion II", ds="20170127", x=0, y=0, type='File', name='Compiled Table Staining')

x <- fread(dt$fileList[1])
y <- fread(dt2$fileList[1])

convertColsToNumeric(x, exceptions=c('Surface'))
lapply(x, function(q){paste(class(q), mode(q))})
convertColsToNumeric(y, exceptions=c('Surface'))
lapply(y, function(q){paste(class(q), mode(q))})

x <- x[order(Surface,Rep,Point,Time)]
y <- y[order(Surface,Rep,Point)]

duh <- x[Surface=='EpCAM' & Rep==1 & Point==716]
duh <- duh[order(Time)]
plot(duh$Freq, duh$FilteredSignal, type='l')
lines(duh$Freq, 5000*duh$FilteredSignal, col='red')

cl <- assignToClusters(log10(y$Blue - min(y$Blue) + 1))
plotClusters(cl$data$x, cl$data$Cluster.Clean)
y[,Cluster:=cl$data$Cluster.Clean]

x[,Time:=(Time-1)*0.035]
x[,Freq:=getFrequencies(Time)$f]

# Look up the cell type 1 is LO blue, 2 is HI blue
x[, CellType:=y[y$Surface==.BY[1] & y$Rep==.BY[2] & y$Point==.BY[3],'Cluster'][1], by=.(Surface,Rep,Point)]
x[, Color:=c('red','blue')[match(CellType, c(1,2))]]

# Analyze the Signals
x[,maxSignal:=quantile(Mean,0.95), by=.(Surface, Rep, Point)]
x <- x[maxSignal > 4000]
x[,Signal:=Mean/maxSignal, by=.(Surface, Rep, Point)]
x[,FilteredSignal:=filterVector(Signal, W=0.05), by=.(Surface, Rep, Point)]

set.seed(12345)
temp <- x[Surface=='EpCAM' & Rep==1,]
chosenIds <- sample(unique(temp$Point), 100)
temp <- temp[Point %in% chosenIds]

temp <- temp[order(Surface,Rep,Point,Time)]
plot(c(), c(), xlab='Frequency [Hz]', ylab='Normalized Adhesion Signal', xlim=rev(range(temp$Freq)), ylim=c(0,1.5), type='l')
temp[,data.table.lines(x=Freq, y=FilteredSignal, type='l', col=Color[1]), by=.(Surface,Rep,Point)]

for(p in unique(temp$Point))
{
	temp2 <- temp[Point==p]
	plot(temp2$Freq, temp2$FilteredSignal, type='l', xlim=rev(range(temp2$Freq)), ylim=c(0,1.5), color=temp2$Color[1], main=paste(p, ', [', temp2$X[1], ',', temp2$Y[1], ']'), xlab='Frequency [Hz]', ylab='Normalized Adhesion Signal')
	print(paste0(p, ': ', as.numeric(quantile(temp2$Signal, 0.95))))
}

getSignalMat <- function(signal)
{
	return(repmat(as.matrix(signal), m=1, n=length(signal)))
}

getPredMat <- function(freq)
{
	ret <- matrix(rep(1.0, length(freq)^2), ncol=length(freq))
	ret[upper.tri(ret)] <- 0
	return(ret)
}

# par is c(A) where A is the amplitude of the HI signal
getDiff <- function(par, signal, predMat, retIndex=FALSE)
{
	signalMat <- getSignalMat(signal)
	diff.2 <- (signalMat-par[1]*predMat)^2
	if(retIndex)
	{
		return(which.min(colSums(diff.2)))
	}
	else
	{
		return(min(colSums(diff.2)))
	}
}

# par is c(A) where A is the amplitude of the HI signal
getDiffWrapper <- function(amplitudes, freq, signal, ct, c, by)
{
	print(by)
	predMat <- getPredMat(freq)
	diffs <- sapply(as.list(amplitudes), getDiff, signal, predMat)
	amp <- amplitudes[which.min(diffs)]
	freq <- freq[getDiff(amp, signal, predMat, retIndex=TRUE)]
	return(list(amp=amp, freq=freq, CellType=ct, Color=c))
}

# Test getting the best fit freq
l(amp, freq) %=% getDiffWrapper(seq(0.5, 1.5, by=0.05), freq=duh$Freq, signal=duh$FilteredSignal)

summary <- temp[,getDiffWrapper(seq(0.5, 1.5, by=0.05), freq=duh$Freq, signal=duh$FilteredSignal, ct=CellType[1], c=Color[1], by=.BY), by=.(Surface, Rep, Point)]

getDiff(1.1, duh$FilteredSignal, getPredMat(duh$FilteredSignal))

optim(par=c(1), method='L-BFGS-B', lower=c(0.5, 1.5), getDiff, duh$FilteredSignal, getPredMat(duh$FilteredSignal))

getFit <- function(freq, signal)
{
	initialGuessI <- which(signal > 0.5)[1]
	initialGuess <- c(1,freq[initialGuessI])
}
