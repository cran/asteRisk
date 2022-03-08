## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(message=FALSE, fig.path='figures/')
hasData <- requireNamespace("asteRiskData", quietly = TRUE)
if (!hasData) {                                                                     
  evaluateExtraChunks <- FALSE                                             
  msg <- strwrap("Note: some examples in this vignette require that the
                 `asteRiskData` package be installed. The system
                  currently running this vignette does not have that package
                  installed, so code examples will not be evaluated.")
  msg <- paste(msg, collapse="\n")
  message(msg)                                                                    
} else {
  evaluateExtraChunks <- TRUE
  library(asteRiskData)
}

## ----tidy = TRUE, eval = FALSE------------------------------------------------
#  install.packages("asteRisk")

## ----tidy = TRUE, eval = TRUE-------------------------------------------------
library(asteRisk)

## ----tidy = TRUE--------------------------------------------------------------
# Read the provided file with 29 benchmark TLE, which contains objects with
# a variety of orbital parameters

test_TLEs <- readTLE(paste0(path.package("asteRisk"), "/testTLE.txt"))
# TLE number 17 contains a state vectors for Italsat 2
test_TLEs[[17]]

# It is also possible to directly parse a character vector with 2 or 3 elements,
# where each element is a string representing a line of the TLE, to obtain the
# same result:
italsat2_lines <- c("ITALSAT 2",
"1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600",
"2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119")

italsat2_TLE <- parseTLElines(italsat2_lines)
italsat2_TLE

test_TLEs[[17]]$inclination == italsat2_TLE$inclination

## ----tidy = TRUE--------------------------------------------------------------
# Read the provided test RINEX navigation files for both GPS and GLONASS:

testGPSnav <- readGPSNavigationRINEX(paste0(path.package("asteRisk"), 
"/testGPSRINEX.txt"))
testGLONASSnav <- readGLONASSNavigationRINEX(paste0(path.package("asteRisk"), 
"/testGLONASSRINEX.txt"))

# Count the number of positional messages in each file:

length(testGPSnav$messages)
length(testGLONASSnav$messages)

## ----tidy = TRUE--------------------------------------------------------------
# Element 11 of the set of test TLE contains an orbital state vector for
# satellite MOLNIYA 1-83, launched from the USSR in 1992 and decayed in 2007

molniya <- test_TLEs[[11]]
1/molniya$meanMotion

# From the inverse of the mean motion, we can see that the orbital period is 
# approximately half a day, in accordance with a Molniya orbit
# Let´s use the SDP4 model to calculate the position and velocity of the
# satellite for a full orbital period every 10 minutes. It is important to
# provide the mean motion in radians/min, the inclination, anomaly, 
# argument of perigee and longitude of the ascending node in radians, and
# the target time as an increment in minutes for the epoch time

targetTimes <- seq(0, 720, by=10)

results_position_matrix <- matrix(nrow=length(targetTimes), ncol=3)
results_velocity_matrix <- matrix(nrow=length(targetTimes), ncol=3)

for(i in 1:length(targetTimes)) {
    new_result <- sgdp4(n0=molniya$meanMotion*((2*pi)/(1440)),
                        e0=molniya$eccentricity,
                        i0=molniya$inclination*pi/180,
                        M0=molniya$meanAnomaly*pi/180,
                        omega0=molniya$perigeeArgument*pi/180,
                        OMEGA0=molniya$ascension*pi/180,
                        Bstar=molniya$Bstar,
                        initialDateTime=molniya$dateTime, targetTime = targetTimes[i])
    results_position_matrix[i,] <- new_result[[1]]
    results_velocity_matrix[i,] <- new_result[[2]]
}
last_molniya_propagation <- new_result
results_position_matrix = cbind(results_position_matrix, targetTimes)
colnames(results_position_matrix) <- c("x", "y", "z", "time")

# Let´s verify that the SDP4 algorithm was automatically chosen

last_molniya_propagation$algorithm

## ----tidy = TRUE, eval = FALSE------------------------------------------------
#  # We can visualize the resulting trajectory using a plotly animation to confirm
#  # that indeed a full revolution was completed and that the orbit is highly
#  # eccentric.
#  
#  library(plotly)
#  library(lazyeval)
#  library(dplyr)
#  
#  # In order to create the animation, we must first define a function to create
#  # the accumulated dataframe required for the animation
#  
#  accumulate_by <- function(dat, var) {
#      var <- f_eval(var, dat)
#      lvls <- plotly:::getLevels(var)
#      dats <- lapply(seq_along(lvls), function(x) {
#          cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
#      })
#      bind_rows(dats)
#  }
#  
#  accumulated_df <- accumulate_by(as.data.frame(results_position_matrix), ~time)
#  
#  orbit_animation <- plot_ly(accumulated_df, x = ~x, y=~y, z=~z, type = "scatter3d",
#                             mode="marker", opacity=0.8, line=list(width = 6,
#                                                                   color = ~time,
#                                                                   reverscale = FALSE),
#                             frame= ~frame)
#  
#  orbit_animation <- animation_opts(orbit_animation, frame=50)
#  
#  orbit_animation <- layout(orbit_animation, scene = list(
#      xaxis=list(range=c(min(results_position_matrix[,1]), max(results_position_matrix[,1]))),
#      yaxis=list(range=c(min(results_position_matrix[,2]), max(results_position_matrix[,2]))),
#      zaxis=list(range=c(min(results_position_matrix[,3]), max(results_position_matrix[,3])))
#  ))
#  
#  orbit_animation

## ----tidy = TRUE, eval = evaluateExtraChunks----------------------------------
# Let us convert the last propagation previously calculated for the MOLNIYA 1-83
# satellite into the ITRF frame. In order to do so, it is required to provide
# a date-time string indicating the time for the newly calculated position and
# velocity. Since this was 720 minutes after the epoch for the original state
# vector, we can just add 12 hours to it 

molniya$dateTime

new_dateTime <- "2006-06-25 12:33:43"

ITRF_coordinates <- TEMEtoITRF(last_molniya_propagation$position,
                               last_molniya_propagation$velocity,
                               new_dateTime)

# Let us now convert the previously calculated set of TEME coordinates to
# geodetic latitude and longitude

geodetic_matrix <- matrix(nrow=nrow(results_position_matrix), ncol=3)

for(i in 1:nrow(geodetic_matrix)) {
    new_dateTime <- as.character(as.POSIXct(molniya$dateTime, tz="UTC") + 60*targetTimes[i])
    new_geodetic <- TEMEtoLATLON(results_position_matrix[i, 1:3]*1000,
                                 new_dateTime)
    geodetic_matrix[i,] <- new_geodetic
}

colnames(geodetic_matrix) <- c("latitude", "longitude", "altitude")

# We can then visualize the ground track of the satellite

library(ggmap)

ggmap(get_map(c(left=-180, right=180, bottom=-80, top=80))) +
  geom_segment(data=as.data.frame(geodetic_matrix), 
               aes(x=longitude, y=latitude, 
                   xend=c(tail(longitude, n=-1), NA), 
                   yend=c(tail(latitude, n=-1), NA)), 
               na.rm=TRUE) +
  geom_point(data=as.data.frame(geodetic_matrix), aes(x=longitude, y=latitude), 
             color="blue", size=0.3, alpha=0.8)

## ----tidy = TRUE, eval = evaluateExtraChunks----------------------------------
# The HPOP requires as input the satellite mass, the effective areas subjected
# to solar radiation pressure and atmospheric drag, and the drag and
# reflectivity coefficients. 
# The mass and cross-section of Molniya satellites are approximately 1600 kg and
# 15 m2, respectively. We will use the cross-section to approximate the
# effective areafor both atmospheric drag and radiation pressure.
# Regarding the drag and reflectivity coefficients, while their values vary
# for each satellite and orbit, 2.2 and 1.2 respectively can be used as
# approximations.

molniyaMass <- 1600
molniyaCrossSection <- 15
molniyaCd <- 2.2
molniyaCr <- 1.2

# As initial conditions, we will use the initial conditions provided in the
# same TLE for MOLNIYA 1-83 used previously for the SGP4/SDP4 propagator.
# We first need to calculate the initial position and velocity in the GCRF
# ECI frame of reference from the provided orbital elements. 
# As an approximation, we will use the results obtained for t = 0 with the
# SGP4/SDP4 propagator. We convert those into the GCRF frame of reference.
# It should be noted that such an approximation introduces an error due to
# a mismatch between the position derivative calculated at each propagation
# point through SGP4/SDP4 and the actual velocity of the satellite.

GCRF_coordinates <- TEMEtoGCRF(results_position_matrix[1,1:3]*1000,
                               results_velocity_matrix[1,1:3]*1000, 
                               molniya$dateTime)

initialPosition <- GCRF_coordinates$position
initialVelocity <- GCRF_coordinates$velocity

# Let´s use the HPOP to calculate the position each 2 minutes during a period
# of 3 hours

targetTimes <- seq(0, 10800, by=120)

hpop_results <- hpop(initialPosition, initialVelocity, molniya$dateTime,
                     targetTimes, molniyaMass, molniyaCrossSection,
                     molniyaCrossSection, molniyaCr, molniyaCd)

# Now we can calculate and plot the corresponding geodetic coordinates

geodetic_matrix_hpop <- matrix(nrow=nrow(hpop_results), ncol=3)

for(i in 1:nrow(geodetic_matrix_hpop)) {
    new_dateTime <- as.character(as.POSIXct(molniya$dateTime, tz="UTC") + targetTimes[i])
    new_geodetic <- GCRFtoLATLON(as.numeric(hpop_results[i, 2:4]), new_dateTime)
    geodetic_matrix_hpop[i,] <- new_geodetic
}

colnames(geodetic_matrix_hpop) <- c("latitude", "longitude", "altitude")

library(ggmap)

ggmap(get_map(c(left=-180, right=180, bottom=-80, top=80))) +
  geom_segment(data=as.data.frame(geodetic_matrix_hpop), 
               aes(x=longitude, y=latitude, 
                   xend=c(tail(longitude, n=-1), NA), 
                   yend=c(tail(latitude, n=-1), NA)), 
               na.rm=TRUE) +
  geom_point(data=as.data.frame(geodetic_matrix_hpop), aes(x=longitude, y=latitude), 
             color="blue", size=0.3, alpha=0.8)


