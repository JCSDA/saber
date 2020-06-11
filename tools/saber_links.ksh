#!/usr/bin/env ksh
# ----------------------------------------------------------------------
# Korn shell script: saber_links
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
# ----------------------------------------------------------------------

# Parameters
modeldata=$1
testdata=$2
test=$3

# Make directory
mkdir -p ${testdata}
# Link members
i=0
ne=1
while [ ${i} -lt ${ne} ] ; do
   # Update
   let i=i+1
   echo "Link member "${i}

   # Copy and typeset
   i3=$i
   typeset -RZ3 i3
   i4=$i
   typeset -RZ4 i4
   i6=$i
   typeset -RZ4 i6

   # AROME
   if test ${test} = "bump_aro" ; then
      ne=10
      xp=7G0N
      date=20131221H00P
      for timeslot in "02" "03" "04" ; do
         ln -sf ${modeldata}/${test}/${xp}/${date}/member_${i3}/forecast/ICMSHAROM+00${timeslot}.nc ${testdata}/ens1_${timeslot}_${i6}.nc
      done
   fi

   # ARPEGE
   if test ${test} = "bump_arp" ; then
      ne=10
      xp=877D
      date=20170114H00A
      for timeslot in "00" "06" ; do
         ln -sf ${modeldata}/${test}/${xp}/${date}/ensemble4D/${i3}/ICMSHARPE+00${timeslot}.nc ${testdata}/ens1_${timeslot}_${i6}.nc
      done
   fi

   # FV3
   if test ${test} = "bump_fv3" ; then
      ne=80
      date=20190612
      hh=12
      resol=c0384
      for itile in $(seq 1 1 6); do
         ln -sf ${modeldata}/${test}/${date}_${hh}/mem${i3}/${date}.${hh}0000.bmat.fv_core.res.tile${itile}.nc ${testdata}/ens1_01_${i6}_tile${itile}.nc
      done
   fi

   # GEM
   if test ${test} = "bump_gem" ; then
      ne=10
      date=2014101706
      ln -sf ${modeldata}/${test}/${date}_006_${i4}.nc ${testdata}/ens1_00_${i6}.nc

      if test ${i} = 1 ; then
         j4=1
         typeset -RZ4 j4
         for string in 'kfc' 'kuo' ; do
            for string2 in 'BLAC62' 'BOUJO' ; do
               k4=1
               typeset -RZ4 k4
               while [ ${k4} -le 64 ] ; do
                  ln -sf  ${modeldata}/${test}/member_${string}_${string2}_${k4}.nc ${testdata}/ens1_00_${j4}_${k4}.nc
                  let k4=k4+1
               done
               let j4=j4+1
            done
         done
      fi
   fi

   # GEOS
   if test ${test} = "bump_geos" ; then
      ne=224
      ts=0
      typeset -RZ2 ts
      while [ ${ts} -lt 13 ] ; do
         ln -sf ${modeldata}/${test}/mem${i3}/e200_C180.prog.eta.20180415_${ts}z.nc4 ${testdata}/ens1_${ts}_${i6}.nc
         let ts=ts+1
      done
   fi

   # GFS
   if test ${test} = "bump_gfs" ; then
      ne=10
      date=2014040100
      ln -sf ${modeldata}/${test}/sfg_${date}_fhr06s_mem${i3}.nc4 ${testdata}/ens1_00_${i6}.nc
   fi

   # IFS
   if test ${test} = "bump_ifs" ; then
      ne=10
      ln -sf ${modeldata}/${test}/member_${i}.nc ${testdata}/ens1_01_${i6}.nc
   fi

   # MPAS
   if test ${test} = "bump_mpas" ; then
      ne=10
      ln -sf ${modeldata}/${test}/x1.40962.output.2012-06-25_21.00.00.e${i}.nc ${testdata}/ens1_01_${i6}.nc
   fi

   # NEMO
   if test ${test} = "bump_nemovar" ; then
      ne=19
      ln -sf ${modeldata}/${test}/ENSEMBLES/ECMWF/goqu/opa${i}/goqu_20110605_000000_restart.nc ${testdata}/ens1_01_${i6}.nc
   fi
   if test ${test} = "bump_cera-20c" ; then
      ne=9
      j4=$i
      typeset -RZ4 j4
      for date in "20090215" "20090216" "20090217" "20090218" "20090219" "20090221" "20090222" "20090223" "20090224" ; do
         ln -sf ${modeldata}/${test}/CERA-20C/member_${date}+00_${i}.nc ${testdata}/ens1_01_${j4}.nc
         let j4=j4+9
      done
   fi

   # RES
   if test ${test} = "bump_res" ; then
      ne=101
      ln -sf ${modeldata}/${test}/Ens_${i}.nc ${testdata}/ens1_01_${i6}.nc
   fi

   # WRF
   if test ${test} = "bump_wrf" ; then
      ne=8
      date=2017-07-28_06:00:00
      ln -sf ${modeldata}/${test}/wrfout_d01_${date}.${i3} ${testdata}/ens1_01_${i6}.nc
   fi

   # Exit
   if test ${i} = ${ne} ; then
      break
   fi
done

# Link wind fields

# Create grid
echo "Link wind fields"

# GEOS
if test ${test} = "bump_geos" ; then
   ts=0
   typeset -RZ2 ts
   while [ ${ts} -lt 12 ] ; do
      ln -sf ${modeldata}/${test}/avg/e200_C180.prog.eta.20180415_${ts}z.nc4 ${testdata}/wind_${ts}.nc
      let ts=ts+1
   done
   ln -sf ${modeldata}/${test}/avg/e200_C180.prog.eta.20180416_00z.nc4 ${testdata}/wind_${ts}.nc
fi

# Create grid
echo "Link grid"

# AROME
if test ${test} = "bump_aro" ; then
   grid=${testdata}/grid.nc
   saved_grid=${modeldata}/${test}/${xp}/grid.nc
   if test -e ${saved_grid} ; then
      # Copy grid.nc from data
      ln -sf ${saved_grid} ${grid}
   else
      # Generate grid.nc with EPyGrAM
      origin=${modeldata}/${test}/${xp}/${date}/member_001/forecast/ICMSHAROM+0003
      grid=${testdata}/grid.nc
      rm -f ${grid}
      cat<<EOFNAM >epygram_request.py
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import epygram
epygram.init_env()
r = epygram.formats.resource("${origin}", "r")
T = r.readfield("S001TEMPERATURE")
if T.spectral:
    T.sp2gp()
gd = T.geometry.dimensions
tab = T.getdata()
tab[...] = 0.0
tab[gd['Y_CIoffset']:
gd['Y_CIoffset']+2*gd['Y_Iwidth']+gd['Y_Czone'],
gd['X_CIoffset']:
gd['X_CIoffset']+2*gd['X_Iwidth']+gd['X_Czone']] = 0.5
tab[gd['Y_CIoffset']+gd['Y_Iwidth']:
gd['Y_CIoffset']+gd['Y_Iwidth']+gd['Y_Czone'],
gd['X_CIoffset']+gd['X_Iwidth']:
gd['X_CIoffset']+gd['X_Iwidth']+gd['X_Czone']] = 1.0
T.setdata(tab)
mapfac = T.geometry.map_factor_field()
rout = epygram.formats.resource("${grid}", "w", fmt="netCDF")
T.fid["netCDF"]="cmask"
mapfac.fid["netCDF"]="mapfac"
rout.behave(flatten_horizontal_grids=False)
rout.writefield(T)
rout.writefield(mapfac)
rout.close()
EOFNAM
      python epygram_request.py
      rm -f epygram_request.py
   fi
fi

# ARPEGE
if test ${test} = "bump_arp" ; then
   grid=${testdata}/grid.nc
   saved_grid=${modeldata}/${test}/${xp}/grid.nc
   if test -e ${saved_grid} ; then
      # Copy grid.nc from data
      ln -sf ${saved_grid} ${grid}
   else
      # Generate grid.nc with EPyGrAM
      origin=${modeldata}/${test}/${xp}/${date}/ensemble4D/001/ICMSHARPE+0000
      rm -f ${grid}
      cat<<EOFNAM >epygram_request.py
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import epygram
epygram.init_env()
r = epygram.formats.resource("${origin}", "r")
T = r.readfield("S001TEMPERATURE")
if T.spectral:
    T.sp2gp()
mapfac = T.geometry.map_factor_field()
rout = epygram.formats.resource("${grid}", "w", fmt="netCDF")
rout.behave(flatten_horizontal_grids=False)
mapfac.fid["netCDF"]="mapfac"
rout.writefield(mapfac)
rout.close()
EOFNAM
      python epygram_request.py
      rm -f epygram_request.py
   fi
fi

# FV3
if test ${test} = "bump_fv3" ; then
   # Copy grid.nc
   origin_h=${modeldata}/${test}/${date}_${hh}/fv3files/fv3grid_${resol}.nc4
   origin_v=${modeldata}/${test}/${date}_${hh}/fv3files/akbk64.nc4
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v flons,flats ${origin_h} ${grid}
   ncks -A -v ak,bk ${origin_v} ${grid}
fi

# GEM
if test ${test} = "bump_gem" ; then
   # Generate grid with ncks
   origin=${testdata}/ens1_00_0001.nc
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v lat,lon,lev,ap,b ${origin} ${grid}
fi

# GEOS
if test ${test} = "bump_geos" ; then
   # Generate grid with ncks and ncwa
   origin=${testdata}/wind_00.nc
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v lons,lats,lev,delp ${origin} ${grid}
fi

# GFS
if test ${test} = "bump_gfs" ; then
   # Generate grid with ncks
   origin=${testdata}/ens1_00_0001.nc
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v latitude,longitude,level,ak,bk ${origin} ${grid}
fi

# IFS
if test ${test} = "bump_ifs" ; then
   # Generate grid.nc with ncks
   origin=${testdata}/ens1_01_0001.nc
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v latitude,longitude,level ${origin} ${grid}

   # Add pressure profile to grid.nc with ncl
   # Copy the full array found on http://www.ecmwf.int/en/forecasts/documentation-and-support/${nflevg}-model-levels into an ascii file "L${nflevg}") where ${nflevg} denotes the number of levels
   nflevg=`ncdump -h ${grid} | grep "level =" | gawk '{print $3}'`
   if test -e "${modeldata}/${test}/L${nflevg}" ; then
      # Remove level 0 and extract pf
      sed '1d' ${modeldata}/${test}/L${nflevg} | gawk '{print $5}' > pf_L${nflevg}

      # Insert pf into grid.nc
      cat<<EOFNAM >pf_into_grid.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

pf = asciiread("pf_L${nflevg}",${nflevg},"float")*1.0e2
data = addfile("${grid}","w")
level = data->level
pf!0 = "level"
pf&level = level
pf@units = "Pa"
pf@long_name = "pressure at full levels"
pf@missing_value = -999
pf@_FillValue = -999
data->pf = pf

end
EOFNAM
      ncl pf_into_grid.ncl

      # Cleaning
      rm -f pf_into_grid.ncl
      rm -f pf_L${nflevg}
   else
      echo "Please copy the full array found on http://www.ecmwf.int/en/forecasts/documentation-and-support/"${nflevg}"-model-levels into an ascii file \"L"${nflevg}"\""
   fi
fi

# MPAS
if test ${test} = "bump_mpas" ; then
   origin=${modeldata}/${test}/x1.40962.init.2018-04-15_00.00.00.nc
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v latCell,lonCell ${origin} ${grid}
   ncwa -O -v pressure_base -a Time,nCells ${origin} pressure.nc
   ncks -A -v pressure_base pressure.nc ${grid}
   rm -f pressure.nc
fi

# NEMO
if test ${test} = "bump_nemovar" ; then
   origin=${modeldata}/${test}/mesh_mask
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v nav_lat,nav_lon,tmask,e1t,e2t,e3t ${origin} ${grid}
fi
if test ${test} = "bump_cera-20c" ; then
   origin=${modeldata}/${test}/mesh_mask
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v nav_lat,nav_lon,tmask,e1t,e2t,e3t ${origin} ${grid}
fi

# RES
if test ${test} = "bump_res" ; then
   origin=${modeldata}/${test}/MyGrid.nc
   grid=${testdata}/grid.nc
   rm -f ${grid}
   cp -f ${origin} ${grid}
fi

# WRF
if test ${test} = "bump_wrf" ; then
   origin=${testdata}/ens1_01_0001.nc
   grid=${testdata}/grid.nc
   rm -f ${grid}
   ncks -O -v XLONG,XLAT ${origin} ${grid}
   ncwa -O -v PB -a Time,south_north,west_east ${origin} pressure.nc
   ncks -A -v PB pressure.nc ${grid}
   rm -f pressure.nc
fi
