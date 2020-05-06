
CircularTukeyDepth <- function(x){

AM <- circular::circular((median.circular(x)+pi), modulo="2pi")
x <- as.vector(na.omit(replace(as.vector(x),as.vector(x)==as.vector(AM), NA)))
x2 <- as.matrix(sort(circular( (x-AM), modulo="2pi")))

AnticlockRank <- as.matrix(seq(1,length(x2), by=1))
ClockRank     <- as.matrix(rev(seq(1,length(x2), by=1)))
Combined <- cbind(AnticlockRank,ClockRank)
Tukeyway <- numeric(length(x2))
for(i in 1:length(x2)){
  Tukeyway[i]   <- Combined[i,][which.min((Combined[i,]))]
}
OuterInward   <-  as.matrix(Tukeyway)
TukeyRanking <- as.matrix(cbind(circular((x2+AM), modulo="2pi"),OuterInward))
colnames(TukeyRanking) <- c("observations",  "depth")

data <- TukeyRanking
data <- as.matrix(data)
CTM <- which(data[,2]>=which.max(data[,2]))
CTM <- circular(mean( circular(data[c(CTM), 1], modulo = "2pi")), modulo="2pi")

n <- length(x)
depthofmedian <- round(((1+n)/2)-0.1)
depthofquartiles <- (1+depthofmedian)/2

if (depthofquartiles%%1==0) {
  quartiles <- which(data[,2] == round(1+depthofmedian)/2)
  qA <-  circular(as.vector(data[quartiles[1],1]),modulo="2pi")
  qC <-  circular(as.vector(data[quartiles[2],1]),modulo="2pi")
}
else  {
  depthq1 <- depthofquartiles+0.5
  depthq2 <- depthofquartiles-0.5
  q1 <- which(data[,2] == depthq1)
  q2 <- which(data[,2] == depthq2)
  qA  <- circular(mean(circular(data[c(q1[1],c(q2[1])),1], modulo="2pi")), modulo = "2pi")
  qC  <- circular(mean(circular(data[c(q1[2],c(q2[2])),1], modulo="2pi")), modulo = "2pi")
}

IQRdepth <- which(data[,2] >= depthofquartiles)
IQR <- c(data[IQRdepth,1],qA,qC)
IQRange <- range(c(qA,qC))

TukeyDepthValues = list()
  TukeyDepthValues$depth = TukeyRanking
  TukeyDepthValues$median = as.circular(CTM)
  TukeyDepthValues$iqr = IQRange

  return(TukeyDepthValues)
}

