dateTimeToMJD <- function(dateTime, timeSystem="UTC") {
    if(!(timeSystem %in% c("UTC", "UT1", "TT", "TDB"))) {
        stop(strwrap("Please choose a time system from UTC, UT1, TT and TDB"))
    }
    nanotimeDateTime <- nanotime(dateTime, format="%Y-%m-%d %H:%M:%E9S", tz = "UTC")
    if(is.na(nanotimeDateTime)) {
        stop("Please provide a complete date-time string in 
             year-month-day hour:minute:second format.")
    }
    dateTimeComponents <- unclass(as.POSIXlt(nanotimeDateTime, tz="UTC"))
    year <- dateTimeComponents$year + 1900
    month <- dateTimeComponents$mon + 1
    day <- dateTimeComponents$mday
    hour <- dateTimeComponents$hour
    minute <- dateTimeComponents$min
    second <- dateTimeComponents$sec
    #MJD <- iauCal2jd(year, month, day, hour, minute, second)$DATE
    UTCSecondsJ2000 <- as.numeric(difftime(as.POSIXct(dateTime, tz="UTC"), J2000_POSIXct, units="secs"))
    MJD_UTC <- UTCSecondsJ2000/86400 + MJD_J2000
    MJD <- MJD_UTC
    if(timeSystem == "TT") {
        hasData()
        MJD <- MJDUTCtoMJDTT(MJD)
    } else if(timeSystem == "UT1") {
        hasData()
        MJD <- MJDUTCtoMJDUT1(MJD)
    } else if(timeSystem == "TDB") {
        hasData()
        MJD <- Mjday_TDB(MJDUTCtoMJDTT(MJD))
    }
    return(MJD)
}

delaunayVariables <- function(CenturiesJ2000_TDB, degreesOutput=FALSE) {
    rev <- if(degreesOutput) 360 else (2*pi)
    delaunayL <- 134.96340251 + (CenturiesJ2000_TDB * (1717915923.2178 + 
                    CenturiesJ2000_TDB * (31.8792 + CenturiesJ2000_TDB *
                    (0.051635 + CenturiesJ2000_TDB * - 0.00024470))))/3600
    delaunayLprime <- 357.52910918 + (CenturiesJ2000_TDB * (129596581.0481 + 
                     CenturiesJ2000_TDB * (-0.5532 + CenturiesJ2000_TDB *
                     (0.000136 + CenturiesJ2000_TDB * - 0.00001149))))/3600
    delaunayF <- 93.27209062 + (CenturiesJ2000_TDB * (1739527262.8478 + 
                   CenturiesJ2000_TDB * (-12.7512 + CenturiesJ2000_TDB *
                   (-0.001037 + CenturiesJ2000_TDB * 0.00000417))))/3600
    delaunayD <- 297.85019547 + (CenturiesJ2000_TDB * (1602961601.2090 + 
                   CenturiesJ2000_TDB * (-6.3706 + CenturiesJ2000_TDB *
                   (0.006593 + CenturiesJ2000_TDB * -0.00003169))))/3600
    delaunayOmega <- 125.04455501 + (CenturiesJ2000_TDB * (-6962890.5431 + 
                     CenturiesJ2000_TDB * (7.4722 + CenturiesJ2000_TDB *
                     (0.007702 + CenturiesJ2000_TDB * -0.00005939))))/3600
    delaunayValues <- c(
        delaunayL,
        delaunayLprime,
        delaunayF,
        delaunayD,
        delaunayOmega
    )
    if(!degreesOutput) delaunayValues <- deg2rad(delaunayValues)
    #delaunayValues <- rem(delaunayValues, rev)
    return(delaunayValues)
}

delaunayToDoodsonVariables <- function(delaunayVariables, GMST, degreesInput=FALSE) {
    halfRev <- if(degreesInput) 180 else pi
    doodsonS <- delaunayVariables[3] + delaunayVariables[5]
    doodsonH <- doodsonS - delaunayVariables[4]
    doodsonP <- doodsonS - delaunayVariables[1]
    doodsonNprime <- -delaunayVariables[5]
    doodsonPs <- doodsonS - delaunayVariables[4] - delaunayVariables[2]
    doodsonTau <- GMST + halfRev - doodsonS
    doodsonVars <- c(doodsonTau, doodsonS, doodsonH, doodsonP, doodsonNprime, doodsonPs)
    return(doodsonVars)
}

xprod <- function(...) {
    ## Code by M. Lundberg and G.V.Welland
    args <- list(...)
    # Check for valid arguments
    if (length(args) == 0) {
        stop("No data supplied")
    }
    len <- unique(sapply(args, FUN=length))
    if (length(len) > 1) {
        stop("All vectors must be the same length")
    }
    if (len != length(args) + 1) {
        stop("Must supply N-1 vectors of length N")
    }
    # Compute generalized cross product by taking the determinant of sub-matrices
    m <- do.call(rbind, args)
    sapply(seq(len),
           FUN=function(i) {
               det(m[,-i,drop=FALSE]) * (-1)^(i+1)
           })
}

JPLephemerides <- function(MJD, timeSystem="UTC", centralBody="SSB", derivatives="acceleration") {
    hasData()
    if(!(timeSystem %in% c("UTC", "UT1", "TT", "TDB"))) {
        stop(strwrap("Please choose a time system from UTC, UT1, TT and TDB"))
    }
    if(!(derivatives %in% c("none", "velocity", "acceleration"))) {
        stop(strwrap("Please choose a derivatives level from none, velocities or acceleration"))
    }
    if(length(centralBody) != 1 | !(centralBody %in% c("SSB", "Mercury", "Venus", "Earth",
                                                       "Moon", "Jupiter", "Saturn", "Uranus",
                                                       "Neptune", "Pluto"))) {
        stop(strwrap("Please choose a single central body from SSB, Mercury, Venus, Earth, Moon, Jupiter, 
                      Saturn, Uranus, Neptune or Pluto"))
    }
    if(timeSystem == "UTC") {
        MJD <- Mjday_TDB(MJDUTCtoMJDTT(MJD))
    } else if(timeSystem == "UT1") {
        MJD <- Mjday_TDB(MJDUTCtoMJDTT(MJDUT1toMJDUTC(MJD)))
    } else if(timeSystem == "TT") {
        MJD <- Mjday_TDB(MJD)
    }
    derivativesOrder <- switch(derivatives, velocity=1, acceleration=2, 0)
    return(JPLephemeridesDE440(MJD, centralBody=centralBody, derivativesOrder=derivativesOrder))
}

registerPlanetaryEphemerides <- function(ephemeridesFile) {
    if(!file.exists(ephemeridesFile)) {
        stop("The file ", ephemeridesFile, " does not exist.")
    }
    if(!endsWith(ephemeridesFile, ".bsp")) {
        stop("The file ", ephemeridesFile, " is not a binary SPICE ephemerides file.")
    }
    if(!grepl("de[0-9]{3}", ephemeridesFile)) {
        stop("The file ", ephemeridesFile, " does not seem to be a JPL DE ephemerides file.")
    }
    ephemeridesVersion <- gsub("^de", "", tolower(substr(ephemeridesFile, 1, 5)))
    checkDEversion(ephemeridesVersion)
    DEfolder <- file.path(asteRiskDataFolder, "DEkernels")
    if(!dir.exists(DEfolder)) {
        dir.create(DEfolder, recursive = TRUE)
    }
    file.copy(ephemeridesFile, file.path(DEfolder, paste("de", ephemeridesVersion, ".bsp", sep="")),
              overwrite=TRUE)
    message("The file ", basename(ephemeridesFile), " has been registered as a JPL DE ephemerides file.")
    return(TRUE)
}

listAvailablePlanetaryEphemerides <- function() {
    DEfolder <- file.path(asteRiskDataFolder, "DEkernels")
    if(!dir.exists(DEfolder)) {
        dir.create(DEfolder, recursive = TRUE)
    }
    availableDEfilenames <- list.files(DEfolder)
    availableDEfilenames <- availableDEfilenames[endsWith(availableDEfilenames, ".bsp")]
    availableDEfilenames <- availableDEfilenames[grepl("de[0-9]{3}", availableDEfilenames)]
    return(availableDEfilenames)
}

.JPLplanetaryEphemerides <- function(ephemerisTimes, ephemeridesVersion) {
    ephemeridesVersion <- gsub("^de", "", tolower(ephemeridesVersion))
    checkDEversion(ephemeridesVersion)
    if(!exists(paste("de", ephemeridesVersion, sep=""), envir = loadedPlanetaryEphemerides)) {
        DEfolder <- file.path(asteRiskDataFolder, "DEkernels")
        if(!dir.exists(DEfolder)) {
            dir.create(DEfolder, recursive = TRUE)
        }
        ephemeridesFilename <- paste("de", ephemeridesVersion, ".bsp", sep="")
        availableDEfilenames <- list.files(DEfolder)
        DEfilePath <- file.path(DEfolder, ephemeridesFilename)
        if(!(ephemeridesFilename %in% availableDEfilenames)) {
            ephemeridesURL <- paste(JPLbinPlanetEphemeridesURL, ephemeridesFilename, sep="")
            filesize <- as.numeric(httr::headers(httr::HEAD(ephemeridesURL))[["Content-Length"]])
            filesizeMB <- filesize/1048576
            currentTimeout <- getOption("timeout")
            expectedTime <- filesizeMB/0.4
            # on tests download from JPL FTP usually goes at around 0.8 MB/s
            downloadPrompt <- readline(prompt=strwrap("Download is ", filesizeMB, 
                                " large. This might take a long time, possibly over ",
                                floor(expectedTime/120), " minutes. Alternatively, you can
                                download the required DE ephemerides file separately and register it 
                                with function registerPlanetaryEphemerides(). Would you like to
                                proceed with download anyway? (y/n): "))
            downloadPrompt <- substr(tolower(downloadPrompt), 1, 1)
            if(downloadPrompt != "y") {
                stop("Download of required planetary ephemerides has been aborted.")
            } else {
                options(timeout=expectedTime)
                download.file(ephemeridesURL, DEfilePath)
                options(timeout=currentTimeout)
                loadPlanetaryEphemerides(DEfilePath)
                DEkernel <- get(paste("de", ephemeridesVersion, sep=""), envir=loadedPlanetaryEphemerides)
                DEconstants <- get(paste("de", ephemeridesVersion, "constants", sep=""), envir=loadedPlanetaryEphemerides)
            }
        } else {
            loadPlanetaryEphemerides(DEfilePath)
            DEkernel <- get(paste("de", ephemeridesVersion, sep=""), envir=loadedPlanetaryEphemerides)
            DEconstants <- get(paste("de", ephemeridesVersion, "constants", sep=""), envir=loadedPlanetaryEphemerides)
        }
    } else {
        DEkernel <- get(paste("de", ephemeridesVersion, sep=""), envir=loadedPlanetaryEphemerides)
        DEconstants <- get(paste("de", ephemeridesVersion, "constants", sep=""), envir=loadedPlanetaryEphemerides)
    }
    segments <- DEkernel$segments
    mercurySegment <- segments[[1]]
    venusSegment <- segments[[2]]
    embSegment <- segments[[3]]
    marsSegment <- segments[[4]]
    jupiterSegment <- segments[[5]]
    saturnSegment <- segments[[6]]
    uranusSegment <- segments[[7]]
    neptuneSegment <- segments[[8]]
    plutoSegment <- segments[[9]]
    sunSegment <- segments[[10]]
    moonSegment <- segments[[11]] # frame of reference centered on EMB
    earthSegment <- segments[[12]] # frame of reference centered on EMB
    #FUTURE: ADD CORRECTIONS TO MARS, VENUS AND MERCURY BARYCENTERS?
    mercuryData <- evaluateType2SPKSegment(mercurySegment, ephemerisTimes)
    venusData <- evaluateType2SPKSegment(venusSegment, ephemerisTimes)
    embData <- evaluateType2SPKSegment(embSegment, ephemerisTimes)
    marsData <- evaluateType2SPKSegment(marsSegment, ephemerisTimes)
    jupiterData <- evaluateType2SPKSegment(jupiterSegment, ephemerisTimes)
    saturnData <- evaluateType2SPKSegment(saturnSegment, ephemerisTimes)
    uranusData <- evaluateType2SPKSegment(uranusSegment, ephemerisTimes)
    neptuneData <- evaluateType2SPKSegment(neptuneSegment, ephemerisTimes)
    plutoData <- evaluateType2SPKSegment(plutoSegment, ephemerisTimes)
    sunData <- evaluateType2SPKSegment(sunSegment, ephemerisTimes)
    moonData <- evaluateType2SPKSegment(moonSegment, ephemerisTimes)
    earthData <- evaluateType2SPKSegment(earthSegment, ephemerisTimes)
    results <- list(
        mercury=mercuryData,
        venus=venusData,
        emb=embData,
        mars=marsData,
        jupiter=jupiterData,
        saturn=saturnData,
        uranus=uranusData,
        neptune=neptuneData,
        pluto=plutoData,
        sun=sunData,
        moon=moonData,
        earth=earthData
    )
    names(results) <- c("Mercury", "Venus", "Earth-Moon Barycenter", 
                        "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", 
                        "Pluto", "Sun", "Moon", "Earth")
    return(results)
}

loadPlanetaryEphemerides <- function(DEKernelPath) {
    DEKernel <- readBinSPK(DEKernelPath)
    DEconstants <- parseBinDEKernelComments(DEKernel$comments)
    ephemeridesVersion <- gsub("^de", "", tolower(substr(basename(DEKernelPath), 1, 5)))
    assign(paste("de", ephemeridesVersion, sep=""), DEKernel, envir=loadedPlanetaryEphemerides)
    assign(paste("de", ephemeridesVersion, "constants", sep=""), DEconstants, envir=loadedPlanetaryEphemerides)
    return(TRUE)
}

# doodsonVariables <- function(MJD_TT) {
#     Tu0 <- (floor(MJD_TT)-51544.5)/36525.0
#     gmst0 <- (6.0/24 + 41.0/(24*60) + 50.54841/(24*60*60))
#     gmst0 <- gmst0 + (8640184.812866/(24*60*60))*Tu0
#     gmst0 <- gmst0 + (0.093104/(24*60*60))*Tu0*Tu0
#     gmst0 <- gmst0 + (-6.2e-6/(24*60*60))*Tu0*Tu0*Tu0
#     r <- 1.002737909350795 + 5.9006e-11*Tu0 - 5.9e-15*Tu0*Tu0
#     gmst <- rem(2*pi*(gmst0 + r * rem(MJD_TT,1)), 2*pi)
#     t = (MJD_TT-51544.5)/365250.0
#     doodsonS = (218.31664562999 + (4812678.81195750 + (-0.14663889 + ( 0.00185140 + -0.00015355*t)*t)*t)*t) * pi/180
#     doodsonH = (280.46645016002 + ( 360007.69748806 + ( 0.03032222 + ( 0.00002000 + -0.00006532*t)*t)*t)*t) * pi/180
#     doodsonP = ( 83.35324311998 + (  40690.13635250 + (-1.03217222 + (-0.01249168 +  0.00052655*t)*t)*t)*t) * pi/180
#     doodsonNprime = (234.95544499000 + (  19341.36261972 + (-0.20756111 + (-0.00213942 +  0.00016501*t)*t)*t)*t) * pi/180
#     doodsonPs = (282.93734098001 + (     17.19457667 + ( 0.04568889 + (-0.00001776 + -0.00003323*t)*t)*t)*t) * pi/180
#     doodsonTau = gmst + pi - doodsonS
#     doodsonVars <- c(doodsonTau, doodsonS, doodsonH, doodsonP, doodsonNprime, doodsonPs)
#     return(doodsonVars)
# }


# alternativeGmst <- function(MJD_TT) {
#     Tu0 <- (floor(MJD_TT)-51544.5)/36525.0
#     gmst0 <- (6.0/24 + 41.0/(24*60) + 50.54841/(24*60*60))
#     gmst0 <- gmst0 + (8640184.812866/(24*60*60))*Tu0
#     gmst0 <- gmst0 + (0.093104/(24*60*60))*Tu0*Tu0
#     gmst0 <- gmst0 + (-6.2e-6/(24*60*60))*Tu0*Tu0*Tu0
#     r <- 1.002737909350795 + 5.9006e-11*Tu0 - 5.9e-15*Tu0*Tu0
#     gmst <- rem(2*pi*(gmst0 + r * rem(MJD_TT,1)), 2*pi)
#     return(gmst)
# }

legendreNormFactor <- function(n, m) {
    if(m==0) {
        delta <- 1
    } else {
        delta <- 0
    }
    sqrt((factorial(n-m)*(2*n+1)*(2-delta))/factorial(n+m))
}

acot <- function(x) {
    return(atan(1/x))
}

acoth <- function(x) {
    return(atanh(1/x))
}

stumpff <- function(x, k) {
    # Using the closed-form solutions for all values of x
    # might be a problem for high k values
    if(k < 0) {
        stop("Negative value of k for Stumpff's functions")
    }
    if(x > 0) {
        z <- sqrt(x)
        c0 <- cos(z)
        if(k == 0) {
            return(c0)
        }
        c1 <- sin(z)/z
        if(k == 1) {
            return(c(c0, c1))
        }
    } else {
        z <- sqrt(-x)
        c0 <- cosh(z)
        if(k == 0) {
            return(c0)
        }
        c1 <- sinh(z)/z
        if(k == 1) {
            return(c(c0, c1))
        }
    }
    ck <- vector(mode="numeric", length=k+1)
    ck[1] <- c0
    ck[2] <- c1
    for(i in 2:(k)) {
        ck[i+1] <- (1/factorial(i-2) - ck[i-1])/x
    }
    return(ck)
}

lagrange_ <- function(x0, y0) {
    x0Min <- min(x0)
    x0Max <- max(x0) 
    intervalRadius <- (x0Max - x0Min)/2
    x0norm <- 2*(x0-x0Min) / (x0Max - x0Min) - 1
    interpolationPoly <- poly.calc(x0norm, y0)
    f <- function(x, derivativesOrder) {
        results <- numeric(derivativesOrder + 1)
        xNorm <- 2*(x-x0Min) / (x0Max - x0Min) - 1
        results[1] <- as.function(interpolationPoly)(xNorm)
        if(derivativesOrder > 0) {
            newPolynom <- interpolationPoly
            for(i in 1:derivativesOrder) {
                newPolynom <- deriv(newPolynom)
                results[i+1] <- as.function(newPolynom)(xNorm)/(intervalRadius^i)
            }
        }
        return(results)
    }
    return(Vectorize(f, "x"))
}

hermite_ <- function(x0, y0, yp0) {
    x0Mean <- mean(x0)
    # x0Min <- min(x0)
    # x0Max <- max(x0)
    # intervalRadius <- (x0Max - x0Min)/2
    # x0norm <- 2*(x0-x0Min) / (x0Max - x0Min) - 1
    x0 <- x0-x0Mean
    A <- polynomial(0)
    B <- polynomial(0)
    for(i in 1:length(x0)) {
        li0 <- poly.calc(x0[-i])/prod(x0[i]-x0[-i])
        li1 <- deriv(li0)
        xMinusxi <- polynomial(c(-x0[i], 1))
        A <- A + (1-2*xMinusxi*predict(li1, x0[i]))*li0^2*y0[i]
        B <- B + xMinusxi*li0^2*yp0[i]
    }
    interpolationPoly <- A + B
    f <- function(x, derivativesOrder) {
        results <- numeric(derivativesOrder + 1)
        # xNorm <- 2*(x-x0Min) / (x0Max - x0Min) - 1
        xNorm <- x - x0Mean
        results[1] <- as.function(interpolationPoly)(xNorm)
        if(derivativesOrder > 0) {
            newPolynom <- interpolationPoly
            for(i in 1:derivativesOrder) {
                newPolynom <- deriv(newPolynom)
                # results[i+1] <- as.function(newPolynom)(xNorm)/(intervalRadius^i)
                results[i+1] <- as.function(newPolynom)(xNorm)
            }
        }
        return(results)
    }
    return(Vectorize(f, "x"))
}


