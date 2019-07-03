load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Parameters:
; -----------

system("mkdir -p fig")
home = getenv("HOME")
datadir = home + "/data/bump/nemo/nemovar"
prefix = "nemo_lct"

; Resources :
; -----------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@cnFillOn = True
res@cnFillMode = "AreaFill"
res@trGridType = "TriangularMesh"
res@cnMonoFillPattern = True
res@cnMonoFillColor = False
res@lbLabelBarOn = True
res@lbOrientation = "vertical"
res@cnInfoLabelOn = False
res@cnLineLabelsOn = False
res@cnLinesOn = False
res@cnNoDataLabelOn = False
res@cnMissingValFillColor = 0
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 1.0
res@cnMinLevelValF = -1.0
res@cnLevelSpacingF = 0.1

res@mpOutlineOn = True
res@mpDataBaseVersion = "MediumRes"
res@mpGridLatSpacingF = 20.0
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "Geophysical"
res@mpLandFillColor = -1
res@mpGridAndLimbDrawOrder = "PreDraw"
res@mpPerimOn = False
res@mpGreatCircleLinesOn = True
;res@mpProjection = "WinkelTripel"
res@mpProjection = "CylindricalEquidistant"
res@mpGridAndLimbOn = True
res@mpGridSpacingF = 360.0
res@mpCenterLonF = 0.0
res@mpCenterLatF = 0.0

; Load and plot data:
;--------------------

data = addfile(datadir + "/" + prefix + "_dirac.nc","r")

nl = data@nl
levs = new(nl,integer)
do il=0,nl-1
   levs(il) = toint(str_get_field(data@levs,il+1,":"))
end do
nlevs = min((/nl,10/))
if (nlevs.gt.1) then
   ilevs_plot = toint(fspan(0,nl-1,nlevs))
else
   ilevs_plot = (/0/)
end if
nv = data@nv
varname = new(nv,string)
do iv=0,nv-1
   varname(iv) = str_get_field(data@varname,iv+1,":")
end do
nts = data@nts
timeslot = new(nts,string)
do its=0,nts-1
   timeslot(its) = str_get_field(data@timeslot,its+1,":")
end do

lon = data->lon
lat = data->lat
res@sfXArray = lon
res@sfYArray = lat

do ilevs=0,nlevs-1
   do iv=0,nv-1
      do its=0,nts-1
         tmp = varname(iv) + "_" + sprinti("%0.2i",its+1)
         var = data->$tmp$(ilevs_plot(ilevs),:)

         output = "../fig/" + prefix + "_dirac_" + varname(iv) + "_ts" + timeslot(its) + "_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
         wks = gsn_open_wks("epsi",output)
         gsn_define_colormap(wks,"BlWhRe")

         res@tiMainString = varname(iv) + " @ lev " + levs(ilevs_plot(ilevs)) + " / ts " + timeslot(its)
         plot = gsn_csm_contour_map(wks,var,res)

         draw(plot)
         frame(wks)
         delete(wks)
         system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
      end do
   end do
end do

end
