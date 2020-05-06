#The CircularBoxplot function produces a box-and-whisker-plot  for circular data, according to the procedure described
#in Buttarazzi D., Pandolfo G., Porzio G.C. (2018). A boxplot for circular data, Biometrics, under review.

CircularBoxplot <- function(A, template="degrees", place="none", units="degrees", marg = "large", shrink = 1.5, H=FALSE, stack=FALSE, constant="optimal") {

  # checking if package are installed, if not they will be installed
  #library(circular) #| install.packages("circular", dep=T)
  #library(plotrix)  #| install.packages("plotrix", dep=T)

  #Check if Median is uniquely defined
  if(is.na(median.circular(A))==T){stop("The median is not unique for this data set. \n \ The circular boxplot is not drawn.")}


  if(constant=="optimal"){
    conc     <- A1inv(rho.circular(A))
    q1       <- qvonmises(0.25, mu=circular(0), kappa = conc)
    me       <- qvonmises(0.5 , mu=circular(0), kappa = conc)
    q3       <- qvonmises(0.75, mu=circular(0), kappa = conc)
    box      <- range(c(q1,me,q3))
    q9965    <- qvonmises(1-(0.007/2), mu=circular(0),kappa = conc)
    q0035    <- qvonmises((0.007/2), mu=circular(0),kappa = conc)
    constant <- range(c(q9965,q3))/box
    }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))


  if(marg == "small") {par(oma=c(0,0,0,0))} else if(marg == "large"){par(mai=c(0.0,0.0,0,0))}

  if(place == "outside") {place="outside"} else if(place == "inside"){place="inside"}


  # inspect and re-specify the input circular vector A.

  if(!is.circular(A)){stop("argument A must be entered as a vector of class circular")}
  set1  <- conversion.circular(A, units = "radians", modulo="2pi", zero=0, rotation="counter")


  # median and IQR

    x <- set1
    AM <- circular((median(x)+pi), modulo="2pi")
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
      #cat("Quartiles
      #  ")
      #print(c(qA,qC))

    }
    else  {
      depthq1 <- depthofquartiles+0.5
      depthq2 <- depthofquartiles-0.5
      q1 <- which(data[,2] == depthq1)
      q2 <- which(data[,2] == depthq2)
      qA  <- mean(circular(data[c(q1[1],c(q2[1])),1], modulo="2pi"))
      qC  <- mean(circular(data[c(q1[2],c(q2[2])),1], modulo="2pi"))
      #cat("Quartiles
      #  ")
      #print(c(qA,qC))

    }


    IQRdepth <- which(data[,2] >= depthofquartiles)
    IQR <- c(data[IQRdepth,1],qA,qC)
    IQRange <- range(c(qA,qC))
    #cat("RANGE
    #    ")
    #print(IQRange)
    set_1 <- set1
    fi <- as.circular(CTM)
    #cat("MEDIAN
    #    ")
    #print(fi)

  #drawing the template
  plot(circular(NA, modulo = "2pi"), cex=0.5, axes=FALSE, shrink=shrink, template=NULL , control.circle=circle.control(col="gray60", lty=2, lwd=0.5))


  #create a general detailed additional template, in degrees or radians, to be put inside or outside the boxplot when template=NULL

  if (is.null(template)){

    if (place=="none"){}

    if((place == "outside") && (units=="radians")){
      draw.arc(0,0,1.15,0,2*pi, col="burlywood4")

      draw.radial.line(1.09,1.15,center=c(0,0), rad(0),   col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0), rad(45),  col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0), rad(90),  col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0), rad(135), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0), rad(180), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0), rad(225), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0), rad(270), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0), rad(315), col="burlywood4")

      draw.radial.line(0,1.08,center=c(0,0),rad(0), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(45), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(90), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(135), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(180), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(225), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(270), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(315), col="azure2",lty=2)
      labelsrad = c(expression(0,frac(pi,4),frac(pi,2),frac(3*pi,4),pi,frac(5*pi,4),frac(3*pi,2),frac(7*pi,4)))
      cosCoord <- cos(rad(circular(c(0,45,90,135,180,225,270,315))))
      sinCoord <- sin(rad(circular(c(0,45,90,135,180,225,270,315))))
      labCoord <- cbind(cosCoord, sinCoord)
      text(1.3*circular(labCoord[,1], units = "radians"),1.3*circular(labCoord[,2], units = "radians"), labels=labelsrad, col="burlywood4", cex=0.6)
    }

    else if((place=="inside") && (units=="radians")){

      draw.arc(0,0,0.7,0,2*pi, col="burlywood4")

      draw.radial.line(0,0.93,center=c(0,0),rad(0), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(45), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(90), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(135), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(180), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(225), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(270), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(315), col="azure2",lty=2)

      draw.radial.line(0.64,0.7,center=c(0,0),rad(0), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(45), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(90), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(135), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(180), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(225), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(270), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(315), col="burlywood4")
      labelsrad = c(expression(0,frac(pi,4),frac(pi,2),frac(3*pi,4),pi,frac(5*pi,4),frac(3*pi,2),frac(7*pi,4)))
      cosCoord <- cos(rad(circular(c(0,45,90,135,180,225,270,315))))
      sinCoord <- sin(rad(circular(c(0,45,90,135,180,225,270,315))))
      labCoord <- cbind(cosCoord, sinCoord)
      text(0.45*circular(labCoord[,1], units = "radians"),0.45*circular(labCoord[,2], units = "radians"), labels=labelsrad, col="burlywood4", cex=0.6)

    }


    if((place == "outside") && (units=="degrees")){
      draw.arc(0,0,1.15,0,2*pi, col="burlywood4")

      draw.radial.line(1.09,1.15,center=c(0,0),rad(0), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0),rad(45), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0),rad(90), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0),rad(135), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0),rad(180), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0),rad(225), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0),rad(270), col="burlywood4")
      draw.radial.line(1.09,1.15,center=c(0,0),rad(315), col="burlywood4")

      draw.radial.line(0,1.08,center=c(0,0),rad(0), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(45), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(90), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(135), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(180), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(225), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(270), col="azure2",lty=2)
      draw.radial.line(0,1.08,center=c(0,0),rad(315), col="azure2",lty=2)
      labelsdeg = c("0","45","90","135","180","225","270","315")

      cosCoord <- cos(rad(circular(c(0,45,90,135,180,225,270,315))))
      sinCoord <- sin(rad(circular(c(0,45,90,135,180,225,270,315))))
      labCoord <- cbind(cosCoord, sinCoord)
      text(1.23*circular(labCoord[,1], units = "radians"), 1.23*circular(labCoord[,2], units = "radians"), labels=labelsdeg, col="burlywood4", cex=0.6)
    }

    else if((place=="inside") && (units=="degrees")){

      draw.arc(0,0,0.7,0,2*pi, col="burlywood4")

      draw.radial.line(0,0.93,center=c(0,0),rad(0), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(45), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(90), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(135), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(180), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(225), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(270), col="azure2",lty=2)
      draw.radial.line(0,0.93,center=c(0,0),rad(315), col="azure2",lty=2)

      draw.radial.line(0.64,0.7,center=c(0,0),rad(0), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(45), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(90), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(135), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(180), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(225), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(270), col="burlywood4")
      draw.radial.line(0.64,0.7,center=c(0,0),rad(315), col="burlywood4")

      labelsdeg = c("0","45","90","135","180","225","270","315")
      cosCoord <- cos(rad(circular(c(0,45,90,135,180,225,270,315))))
      sinCoord <- sin(rad(circular(c(0,45,90,135,180,225,270,315))))
      labCoord <- cbind(cosCoord, sinCoord)
      text(0.5*circular(labCoord[,1], units = "radians"), 0.5*circular(labCoord[,2], units = "radians"), labels=labelsdeg, col="burlywood4", cex=0.6)

    }

  }

  else if(template=="degrees"){
    labelsdeg = c("0","90","180","270")
    cosCoord <- cos(rad(circular(c(0,90,180,270))))
    sinCoord <- sin(rad(circular(c(0,90,180,270))))
    labCoord <- cbind(cosCoord, sinCoord)
    text(0.82*circular(labCoord[,1], units = "radians"),0.82*circular(labCoord[,2], units = "radians"), labels=labelsdeg, cex=0.6)
  }

  else if(template=="radians"){
    labelsrad = c(expression(0,frac(pi,2),pi,frac(3*pi,2)))
    cosCoord <- cos(rad(circular(c(0,90,180,270))))
    sinCoord <- sin(rad(circular(c(0,90,180,270))))
    labCoord <- cbind(cosCoord, sinCoord)
    text(0.65*circular(labCoord[,1], units = "radians"),0.65*circular(labCoord[,2], units = "radians"), labels=labelsrad, cex=0.6)
  }

  else if(template=="geographics"){
    labelsgeo = c("E","N","W","S")
    cosCoord <- cos(rad(circular(c(0,90,180,270))))
    sinCoord <- sin(rad(circular(c(0,90,180,270))))
    labCoord <- cbind(cosCoord, sinCoord)
    text(0.82*circular(labCoord[,1], units = "radians"),0.82*circular(labCoord[,2], units = "radians"), labels=labelsgeo, cex=0.6)
  }

  #drawing the plot

  if(H==T){points(circular(set_1),  cex=0.75)}
  else if(H==F){points(circular(NA, modulo = "2pi"))}
  points(circular(IQR, modulo = "2pi"), cex=1.1, col="white")
  points(0,0,pch=21, bg=4,  cex=1.1)

  # controlling wrap-around effect in case of median at pi (180?)

  if (rad(round(deg(circular((fi+pi), modulo="2pi"))))==0){
    fi <- pi
    AM<- 2*pi}

  else{AM <- rad(round(deg(circular((fi+pi), modulo="2pi"))))}

  # controlling wrap-around effect

  if (range(circular(IQR, modulo = "2pi"))< ((2*pi)/(2*(constant + (1/2)))) ) {

    if (fi<pi) {
      setAnti <- subset(IQR, IQR>=fi & IQR<=AM)
      setClock<- subset(IQR, IQR<=fi | IQR>=AM)
      QAnti   <- rad(round(deg(circular(max(setAnti), modulo="2pi"))))
      Qc      <- QAnti-rad(round(deg(range(circular(IQR, modulo = "2pi")))))
      QClock  <- rad(round(deg(circular(Qc, modulo="2pi"))))

      #box in pale gray
      grid <- seq(Qc, QAnti, by=0.001)
      ngrid <- length(grid)
      for(i in 1:ngrid){
        draw.radial.line(0.9,1.1, center=c(0,0), grid[i], col="gray80" , lwd=2)
      }


      draw.arc(0,0,0.9,QAnti,Qc,col=1,lwd=2)
      draw.arc(0,0,1.1,QAnti,Qc,col=1,lwd=2)

      draw.radial.line(0.9,1.1,center=c(0,0),QAnti,col=1,lwd=2)
      draw.radial.line(0.9,1.1,center=c(0,0),QClock,col=1,lwd=2)

      d <- (rad(round(deg(range(circular(IQR, modulo = "2pi"))))))


      fA<- rad(round(deg(QAnti + d*constant)))
      fC<- rad(round(deg(QClock - d*constant)))

      semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi | as.vector(set_1)>=AM)
      semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi & as.vector(set_1)<=AM)

      semicircleClock <- c(semicircleClock, QClock)
      semicircleAnti <- c(semicircleAnti, QAnti)
      if (fC<0) {
        swc <- subset(semicircleClock, semicircleClock>= rad(round(deg(circular(fC, modulo="2pi"))))| semicircleClock<= QClock)
        swc <- c(swc, QClock)
        whiskerC <- range(circular(swc))
        wC <- QClock-whiskerC
        draw.arc(0,0,1,wC,QClock,col=1,lwd=2, lty=1)
        draw.radial.line(0.95,1.05,center=c(0,0),wC,col=1,lwd=2)
        faroutClock  <- subset(semicircleClock, semicircleClock>=AM & semicircleClock<rad(round(deg(circular(fC, modulo="2pi")))))
      }
      else if (fC>=0 & QClock>=pi){
        swc <- subset(semicircleClock, semicircleClock>=fC)
        swc <- c(swc, QClock)
        wC <- min(swc)
        draw.arc(0,0,1,wC,QClock,col=1,lwd=2, lty=1)
        draw.radial.line(0.95,1.05,center=c(0,0),wC,col=1,lwd=2)
        faroutClock <- subset(semicircleClock, semicircleClock>=AM & semicircleClock<fC)
      }
      else if (fC>=0 & QClock<pi){
        swc <- subset(semicircleClock, semicircleClock>=fC)
        swc <- c(swc, QClock)
        wC <- min(swc)
        draw.arc(0,0,1,wC,QClock,col=1,lwd=2, lty=1)
        draw.radial.line(0.95,1.05,center=c(0,0),wC,col=1,lwd=2)
        faroutClock <- subset(semicircleClock, semicircleClock>=AM | semicircleClock<fC)
      }

      swa <- subset(semicircleAnti, semicircleAnti<=fA)
      swa <- c(swa, QAnti)
      wA <- max(swa)
      draw.arc(0,0,1,wA,QAnti,col=1,lwd=2, lty=1)
      draw.radial.line(0.95,1.05,center=c(0,0),wA,col=1,lwd=2)
      faroutAnti  <- subset(semicircleAnti, semicircleAnti>fA)
    }

    if (fi==pi) {
      setAnti <- subset(IQR, IQR>=fi & IQR<=2*pi)
      setClock<- subset(IQR, IQR<=fi | IQR>=0)
      QAnti   <- rad(round(deg(circular(max(setAnti), modulo="2pi"))))
      Qc      <- QAnti-rad(round(deg(range(circular(IQR, modulo = "2pi")))))
      QClock  <- rad(round(deg(circular(Qc, modulo="2pi"))))

      grid <- seq(Qc, QAnti, by=0.001)
      ngrid <- length(grid)
      for(i in 1:ngrid){
        draw.radial.line(0.9,1.1, center=c(0,0), grid[i], col="gray80" , lwd=2)
      }


      draw.arc(0,0,1.1,QAnti,Qc,col=1,lwd=2)
      draw.arc(0,0,0.9,QAnti,Qc,col=1,lwd=2)

      draw.radial.line(0.9,1.1,center=c(0,0),QAnti,col=1,lwd=2)
      draw.radial.line(0.9,1.1,center=c(0,0),QClock,col=1,lwd=2)

      # defining the whiskers
      d <- (rad(round(deg(range(circular(IQR, modulo = "2pi"))))))
      #cat("IQ-RANGE")
      #print(d)
      fA<- rad(round(deg(QAnti + d*constant)))
      fC<- rad(round(deg(QClock - d*constant)))

      semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi | as.vector(set_1)>= 0)
      semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi & as.vector(set_1)<= 2*pi)

      semicircleClock <- c(semicircleClock, QClock)
      semicircleAnti <- c(semicircleAnti, QAnti)
      if (fC<0) {
        swc <- subset(semicircleClock, semicircleClock>= rad(round(deg(circular(fC, modulo="2pi"))))| semicircleClock<= QClock)
        swc <- c(swc, QClock)
        whiskerC <- range(circular(swc))
        wC <- QClock-whiskerC

        #drawing the whiskers
        draw.arc(0,0,1,wC,QClock,col=1,lwd=2, lty=1)
        draw.radial.line(0.95,1.05,center=c(0,0),wC,col=1,lwd=2)
        faroutClock  <- subset(semicircleClock, semicircleClock>=0 & semicircleClock<rad(round(deg(circular(fC, modulo="2pi")))))
      }
      else if (fC>=0){
        swc <- subset(semicircleClock, semicircleClock>=fC)
        swc <- c(swc, QClock)
        wC <- min(swc)
        draw.arc(0,0,1,wC,QClock,col=1,lwd=2, lty=1)
        draw.radial.line(0.95,1.05,center=c(0,0),wC,col=1,lwd=2)
        faroutClock <- subset(semicircleClock, semicircleClock>=0 | semicircleClock<fC)
      }
      swa <- subset(semicircleAnti, semicircleAnti<=fA)
      swa <- c(swa, QAnti)
      wA <- max(swa)
      draw.arc(0,0,1,wA,QAnti,col=1,lwd=2, lty=1)
      draw.radial.line(0.95,1.05,center=c(0,0),wA,col=1,lwd=2)
      faroutAnti  <- subset(semicircleAnti, semicircleAnti>fA)
      #print(faroutAnti)
    }

    else if (fi>pi) {
      setAnti <- subset(IQR, IQR>=fi | IQR<=AM)
      setClock<- subset(IQR, IQR<=fi & IQR>=AM)
      QClock   <- min(setClock)
      Qa      <- QClock+range(circular(IQR, modulo = "2pi"))
      QAnti  <- rad(round(deg(circular(Qa, modulo="2pi"))))

      grid <- seq(QClock, Qa, by=0.001)
      ngrid <- length(grid)
      for(i in 1:ngrid){
        draw.radial.line(0.9,1.1, center=c(0,0), grid[i], col="gray80" , lwd=2)
      }


      draw.arc(0,0,1.1,QClock,Qa,col=1,lwd=2)
      draw.arc(0,0,0.9,QClock,Qa,col=1,lwd=2)

      draw.radial.line(0.9,1.1,center=c(0,0),QAnti,col=1,lwd=2)
      draw.radial.line(0.9,1.1,center=c(0,0),QClock,col=1,lwd=2)

      # defining the whiskers
      d <- range(circular(IQR, modulo = "2pi"))
      fC<- rad(round(deg(QClock - d*constant)))
      fA<- rad(round(deg(QAnti + d*constant)))
      semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi & as.vector(set_1)>=AM)
      semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi | as.vector(set_1)<=AM)

      semicircleClock <- c(semicircleClock, QClock)
      semicircleAnti <- c(semicircleAnti, QAnti)
      swc <- subset(semicircleClock, semicircleClock>=fC)
      swc <- c(swc, QClock)
      wC <- min(swc)
      draw.arc(0,0,1,wC,QClock,col=1,lwd=2, lty=1)
      draw.radial.line(0.95,1.05,center=c(0,0),wC,col=1,lwd=2)
      faroutClock <- subset(semicircleClock, semicircleClock<fC )

      if (fA>2*pi ) {
        swa <- subset(semicircleAnti, semicircleAnti<= rad(round(deg(circular(fA, modulo="2pi")))) | semicircleAnti>= rad(round(deg(circular(QAnti, modulo="2pi")))))
        swa <- c(swa, QAnti)
        whiskerA <- range(circular(swa))
        wA <- QAnti+whiskerA

        # drawing the whiskers
        draw.arc(0,0,1,QAnti,wA,col=1,lwd=2, lty=1)
        draw.radial.line(0.95,1.05,center=c(0,0),wA,col=1,lwd=2)
        faroutAnti <- subset(semicircleAnti, semicircleAnti> circular(fA, modulo = "2pi") & semicircleAnti <= AM )
      }
      else if (fA<=2*pi & QAnti>=pi) {
        swa <- subset(semicircleAnti, semicircleAnti<=fA)
        swa <- c(swa, QAnti)
        wA <- max(swa)
        draw.arc(0,0,1,wA,QAnti,col=1,lwd=2, lty=1)
        draw.radial.line(0.95,1.05,center=c(0,0),wA,col=1,lwd=2)
        faroutAnti <- subset(semicircleAnti, semicircleAnti>fA | semicircleAnti<= AM  )
      }
      else if (fA<=2*pi & QAnti<pi) {
        swa <- subset(semicircleAnti, semicircleAnti<=fA)
        swa <- c(swa, QAnti)
        wA <- max(swa)
        draw.arc(0,0,1,wA,QAnti,col=1,lwd=2, lty=1)
        draw.radial.line(0.95,1.05,center=c(0,0),wA,col=1,lwd=2)
        faroutAnti <- subset(semicircleAnti, semicircleAnti>fA & semicircleAnti<= AM  )
      }
    }

    # plotting and printing far out values

    faroutvalues1 <- c(faroutClock, faroutAnti)
    compare <- set_1
    faroutvalues2 <- compare[compare %in% faroutvalues1]
    faroutvalues <- as.circular(circular(faroutvalues2), modulo="2pi")
    farout<- as.matrix(faroutvalues2)
    colnames(farout) <- c("Far out values")

    if(H==T){
      points(faroutvalues, cex=0.8, col="white")
      points(faroutvalues, cex=0.8, pch=8) }

    else if(H==F){
      if(stack==T){points(faroutvalues, cex=0.6, stack=stack, bins=500, sep=0.1, pch=8)  }
      else if(stack==F) {points(faroutvalues, cex=0.6, stack=stack, pch=8)}
    }
  }





  #####from here on is in case the range(box)>= (360/2(c+1/2))


  else {
    if (fi<=pi) {
      setAnti <- subset(IQR, IQR>=fi & IQR<=AM)
      setClock<- subset(IQR, IQR<=fi | IQR>=AM)
      QAnti   <- rad(round(deg(circular(max(setAnti), modulo="2pi"))))
      Qc      <- QAnti-range(circular(IQR, modulo = "2pi"))
      QClock  <- rad(round(deg(circular(Qc, modulo="2pi"))))

      grid <- seq(Qc, QAnti, by=0.001)
      ngrid <- length(grid)
      for(i in 1:ngrid){
        draw.radial.line(0.9,1.1, center=c(0,0), grid[i], col="gray80" , lwd=2)
      }

      draw.arc(0,0,1.1,QAnti,Qc,col=1,lwd=2)
      draw.arc(0,0,0.9,QAnti,Qc,col=1,lwd=2)

      draw.radial.line(0.9,1.1,center=c(0,0),QAnti,col=1,lwd=2)
      draw.radial.line(0.9,1.1,center=c(0,0),QClock,col=1,lwd=2)

      # defining the whiskers
      semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi | as.vector(set_1)>=AM)
      semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi & as.vector(set_1)<=AM)

      semicircleClock <- c(semicircleClock, QClock)
      semicircleAnti <- c(semicircleAnti, QAnti)
      if (QClock <= pi){
        swc <- subset(semicircleClock, semicircleClock >= AM | semicircleClock <= QClock)
      }
      else if (QClock > pi){
        swc <- subset(semicircleClock, semicircleClock >= AM & semicircleClock <= QClock)
      }

      whiskerC <- range(circular(swc))
      wC <- QClock-whiskerC

      # drawing the whiskers
      draw.arc(0,0,1,wC,QClock,col=1,lwd=2, lty=1)
      draw.radial.line(0.95,1.05,center=c(0,0),wC,col=1,lwd=2)
      swa <- subset(semicircleAnti, semicircleAnti<=AM)
      wA <- max(swa)
      draw.arc(0,0,1,wA,QAnti,col=1,lwd=2, lty=1)
      draw.radial.line(0.95,1.05,center=c(0,0),wA,col=1,lwd=2)

    }

    else if (fi>=pi) {
      setAnti <- subset(IQR, IQR>=fi | IQR<=AM)
      setClock<- subset(IQR, IQR<=fi & IQR>=AM)
      QClock   <- min(setClock)
      Qa      <- QClock+range(circular(IQR, modulo = "2pi"))
      QAnti  <- rad(round(deg(circular(Qa, modulo="2pi"))))

      grid <- seq(QClock, Qa, by=0.001)
      ngrid <- length(grid)
      for(i in 1:ngrid){
        draw.radial.line(0.9,1.1, center=c(0,0), grid[i], col="gray80" , lwd=2)
      }


      draw.arc(0,0,1.1,QClock,Qa,col=1,lwd=2)
      draw.arc(0,0,0.9,QClock,Qa,col=1,lwd=2)

      draw.radial.line(0.9,1.1,center=c(0,0),QAnti,col=1,lwd=2)
      draw.radial.line(0.9,1.1,center=c(0,0),QClock,col=1,lwd=2)

      semicircleClock <- subset(as.vector(set_1),as.vector(set_1)<=fi & as.vector(set_1)>=AM)
      semicircleAnti <- subset(as.vector(set_1),as.vector(set_1)>=fi | as.vector(set_1)<=AM)
      semicircleClock <- c(semicircleClock, QClock)
      semicircleAnti <- c(semicircleAnti, QAnti)
      swc <- subset(semicircleClock, semicircleClock>=AM)
      wC <- min(swc)

      # drawing the whiskers
      draw.arc(0,0,1,wC,QClock,col=1,lwd=2, lty=1)
      draw.radial.line(0.95,1.05,center=c(0,0),wC,col=1,lwd=2)

      if(QAnti<pi){
        swa <- subset(semicircleAnti, semicircleAnti<=AM & semicircleAnti >= QAnti)
      }
      else if(QAnti>pi){
        swa <- subset(semicircleAnti, semicircleAnti<=AM | semicircleAnti >= QAnti)
      }
      whiskerA <- range(circular(swa))
      wA <- QAnti+whiskerA
      draw.arc(0,0,1,wA,QAnti,col=1,lwd=2, lty=1)
      draw.radial.line(0.95,1.05,center=c(0,0),wA,col=1,lwd=2)
    }

  }
  gradi <- (as.matrix(deg(data[,1])))
  output <- as.matrix(cbind(data,gradi))
  colnames(output) <- c("Obs.Radians", "Ranking", "Obs.Degrees")
  #print(output)

  # drawing an arrow indicating the median
  draw.radial.line(0.905,1.095,center=c(0,0),CTM,col=4,lwd=2)
  arrows.circular(CTM,0.78, col=4, angle=15)

  # output object
  out = list()
  if(exists("faroutvalues")==TRUE){
    if(length(faroutvalues)!=0){out$farout = faroutvalues}
    else{out$farout = c("no far out values detected")}
  }
  else{out$farout = c("no far out values detected")}

  out$constant = constant
  return(invisible(out))

   }



