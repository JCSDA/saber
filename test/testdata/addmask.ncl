load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

data = addfile("grid.nc","w")
lon = data->lon
lat = data->lat
z = data->z
dims = dimsizes(lon)
nx = dims(1)
ny = dims(0)
nz = dimsizes(z)

lmask = new((/nz,ny,nx/),double)
do ix=0,nx-1
  do iy=0,ny-1
    if (sqrt(lon(iy,ix)^2+(5.0*(lat(iy,ix)*(1+0.1*sin(6.28*lon(iy,ix)/60.0))-45.0))^2).gt.100.0) then
       lmask(3,iy,ix) = 1.0
    else
       lmask(3,iy,ix) = 0.0
    end if
    if (sqrt(lon(iy,ix)^2+(5.0*(lat(iy,ix)*(1+0.1*sin(6.28*lon(iy,ix)/60.0))-45.0))^2).gt.80.0) then
       lmask(0,iy,ix) = 1.0
       lmask(2,iy,ix) = 1.0
    else
       lmask(0,iy,ix) = 0.0
       lmask(2,iy,ix) = 0.0
    end if
    if (sqrt(lon(iy,ix)^2+(5.0*(lat(iy,ix)*(1+0.1*sin(6.28*lon(iy,ix)/60.0))-45.0))^2).gt.60.0) then
       lmask(1,iy,ix) = 1.0
    else
       lmask(1,iy,ix) = 0.0
    end if
  end do
end do

lmask!0 = "nz"
lmask!1 = "ny"
lmask!2 = "nx"

data->lmask = lmask

end
