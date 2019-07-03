load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Parameters:
; -----------

system("mkdir -p fig")
home = getenv("HOME")
datadir = home + "/data/bump/arp/877D"
prefix = "arp_877D"
type = "flt"

; Resources :
; -----------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@cnFillOn = True
res@cnFillMode = "RasterFill"
res@cnMonoFillPattern = True
res@cnMonoFillColor = False
res@lbLabelBarOn = True
res@lbOrientation = "vertical"
res@cnInfoLabelOn = False
res@cnLineLabelsOn = False
res@cnLinesOn = False
res@cnNoDataLabelOn = False
res@cnMissingValFillColor = 0

res@mpOutlineOn = True
res@mpDataBaseVersion = "MediumRes"
res@mpGridLatSpacingF = 20.0
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "Geophysical"
res@mpLandFillColor = -1
res@mpGridAndLimbDrawOrder = "PreDraw"
res@mpPerimOn = False
res@mpGreatCircleLinesOn = True
res@mpProjection = "WinkelTripel"
res@mpGridAndLimbOn = True
res@mpGridSpacingF = 360.0
res@mpCenterLonF = 0.0
res@mpCenterLatF = 0.0

xyres = True
xyres@gsnFrame = False
xyres@gsnDraw = False
xyres@tmXBLabelFontHeightF = 0.02
xyres@tmYLLabelFontHeightF = 0.02
xyres@xyLineThicknessF = 5.0
xyres@xyLineColors = (/"red","blue","blue"/)
xyres@xyDashPatterns = (/0,14,14/)
xyres@xyLineDashSegLenF = 0.7
xyres@vpWidthF = 0.4
xyres@vpHeightF = 0.8
xyres@tiXAxisString = "Iterations"

pmres = True
pmres@gsMarkerIndex = 1
pmres@gsMarkerSizeF = 0.01

plres = True
plres@gsLineThicknessF = 3.0

resP = True
resP@gsnFrame = False
resP@gsnPanelXWhiteSpacePercent = 8
resP@gsnPanelYWhiteSpacePercent = 5

; Load and plot data:
;--------------------

data = addfile(datadir + "/" + prefix + "_displ_diag.nc","r")

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
nts = data@nts
timeslot = new(nts,string)
do its=0,nts-1
   timeslot(its) = str_get_field(data@timeslot,its+1,":")
end do

lon_c2 = data->lon_c2
lat_c2 = data->lat_c2
dims = dimsizes(lon_c2)
nc2 = dims(1)
if (type.eq."raw") then
   lon_c2_dst = data->lon_c2_raw
   lat_c2_dst = data->lat_c2_raw
   dist_c2_dst = data->dist_c2_raw
end if
if (type.eq."flt") then
   lon_c2_dst = data->lon_c2_flt
   lat_c2_dst = data->lat_c2_flt
   dist_c2_dst = data->dist_c2_flt
end if
dot = new(nc2,graphic)
line = new(nc2,graphic)
dist_min = 0.0
dist_max = max(dist_c2_dst)
valid = data->valid
dist = data->dist
rhflt = data->rhflt
dims_pro = dimsizes(valid)
niter = dims_pro(2)

do ilevs=0,nlevs-1
   do its=1,nts-1
      output = "../fig/" + prefix + "_displ_profiles_" + type + "_ts" + timeslot(its) + "_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
      wks = gsn_open_wks("epsi",output)
      gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

      xyres@trXMinF = 1
      xyres@trXMaxF = niter
      plot_prof = new(3,graphic)

      xyres@trYMinF = 0.0
      xyres@trYMaxF = max(rhflt)*1.05
      xyres@tiXAxisString = "Filtering length-scale (km)"
      plot_prof(0) = gsn_csm_xy(wks,fspan(1,niter,niter),rhflt(its-1,ilevs_plot(ilevs),:),xyres)

      xyres@trYMinF = 0.0
      xyres@trYMaxF = 105.0
      xyres@gsnYRefLine = 100.0
      xyres@tiXAxisString = "Valid points (%)"
      plot_prof(1) = gsn_csm_xy(wks,fspan(1,niter,niter),valid(its-1,ilevs_plot(ilevs),:)*100.0,xyres)
      delete(xyres@gsnYRefLine)

      xyres@trYMinF = 0.0
      xyres@trYMaxF = max(dist)*1.05
      xyres@tiXAxisString = "Displacement distance (km)"
      plot_prof(2) = gsn_csm_xy(wks,fspan(1,niter,niter),dist(its-1,ilevs_plot(ilevs),:),xyres)

      gsn_panel(wks,plot_prof,(/1,3/),resP)
      frame(wks)
      delete(wks)
      system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
   end do
end do

do ilevs=0,nlevs-1
   do its=1,nts-1
      output = "../fig/" + prefix + "_displ_diag_" + type + "_ts" + timeslot(its) + "_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
      wks = gsn_open_wks("epsi",output)
      gsn_define_colormap(wks,"MPL_brg")

      res@tiMainString = "Displacement @ lev " + levs(ilevs_plot(ilevs)) + " / ts " + timeslot(its)
      plot = gsn_csm_map(wks,res)

      dist_color = 2+toint((dist_c2_dst(its-1,ilevs_plot(ilevs),:)-dist_min)/(dist_max-dist_min)*(61.0-2.0))
      do ic2=0,nc2-1
         pmres@gsMarkerColor = min((/dist_color(ic2),61/))
         plres@gsLineColor = min((/dist_color(ic2),61/))
         dot(ic2) = gsn_add_polymarker(wks,plot,lon_c2(ilevs_plot(ilevs),ic2),lat_c2(ilevs_plot(ilevs),ic2),pmres)
         line(ic2) = gsn_add_polyline(wks,plot,(/lon_c2(ilevs_plot(ilevs),ic2),lon_c2_dst(its-1,ilevs_plot(ilevs),ic2)/),(/lat_c2(ilevs_plot(ilevs),ic2),lat_c2_dst(its-1,ilevs_plot(ilevs),ic2)/),plres)
      end do

      draw(plot)
      frame(wks)
      delete(wks)
      system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
   end do
end do

end
