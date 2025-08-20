\name{registerPlanetaryEphemerides}
\alias{registerPlanetaryEphemerides}
\title{Register a Planetary Ephemerides File}
\description{
This function registers a binary SPICE ephemerides file (with .bsp extension) as a JPL DE ephemerides file. 
The file is copied into a sub directory of the user data directory for asteRisk, and a message is displayed to 
confirm the registration.
}
\usage{
registerPlanetaryEphemerides(ephemeridesFile)
}
\arguments{
  \item{ephemeridesFile}{The path to the ephemerides file to be registered. The file have a .bsp extension, and be a JPL DE ephemerides file.}
}
\value{
Returns TRUE if the file is successfully registered.
}
\examples{
\dontrun{
# Assuming "de430.bsp" is a valid ephemerides file in the current directory
registerPlanetaryEphemerides("de430.bsp")
}
}
\references{
NAIF/SPICE Documentation: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
}