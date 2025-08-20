\name{listAvailablePlanetaryEphemerides}
\alias{listAvailablePlanetaryEphemerides}
\title{List Available Planetary Ephemerides Files}
\description{
This function lists all available JPL DE ephemerides files (with .bsp extension).
More JPL DE ephemerides files can be downloaded from the JPL Horizons website (https://ssd.jpl.nasa.gov/?horizons),
and registered with the \code{\link{registerPlanetaryEphemerides}}.
}
\usage{
listAvailablePlanetaryEphemerides()
}
\arguments{
  None
}
\value{
Returns a vector of filenames of the available JPL DE ephemerides files.
}
\examples{
# List available ephemerides files
listAvailablePlanetaryEphemerides()
}
\references{
NAIF/SPICE Documentation: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
}