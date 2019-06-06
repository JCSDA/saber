# Offline or online usage

The bump.x executable can be used as standalone code, with NetCDF inputs, for the following models:
  - [ARPEGE](http://www.cnrm-game-meteo.fr/spip.php?article121&lang=en)
  - [AROME](http://www.cnrm-game-meteo.fr/spip.php?article120&lang=en)
  - [FV3](https://www.gfdl.noaa.gov/fv3)
  - [GEM](https://en.wikipedia.org/wiki/Global_Environmental_Multiscale_Model)
  - [GEOS](https://gmao.gsfc.nasa.gov/GEOS)
  - [GFS](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs)
  - [IFS](http://www.ecmwf.int/en/research/modelling-and-prediction)
  - [MPAS](https://mpas-dev.github.io)
  - [NEMO](http://www.nemo-ocean.eu)
  - RES
  - [WRF](https://www.mmm.ucar.edu/weather-research-and-forecasting-model)

It can also be used "online" within an other code, using a dedicated interface: [type_bump.F90](../../src/bump/type_bump.F90)
