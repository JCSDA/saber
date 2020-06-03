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
diag = "cor_01_02_01_01"
nl0_plot = 5

; Resources:
; ----------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@cnFillMode = "RasterFill"
res@cnFillOn = True
res@cnMonoFillPattern = True
res@cnMonoFillColor = False
res@cnInfoLabelOn = False
res@cnLineLabelsOn = False
res@cnLinesOn = False
res@cnNoDataLabelOn = False
res@cnMissingValFillColor = 0
res@lbLabelBarOn = True
res@lbOrientation = "vertical"

xyres = True
xyres@gsnFrame = False
xyres@gsnDraw = False
xyres@tmXBLabelFontHeightF = 0.02
xyres@tmYLLabelFontHeightF = 0.02
xyres@tiMainFontHeightF = 0.03
xyres@xyLineThicknessF = 5.0
xyres@xyLineColors = (/"red"/)
xyres@xyDashPatterns = (/0/)
xyres@vpWidthF = 0.4
xyres@vpHeightF = 0.8

resP = True
resP@gsnFrame = False
resP@gsnPanelXWhiteSpacePercent = 8
resP@gsnPanelYWhiteSpacePercent = 5

; Load and plot data:
; -------------------

data = addfile(datadir + "/" + prefix + "_diag.nc","r")
disth = data->disth
vunit = data->vunit
tmp = diag + "_l0rl0_to_l0"
l0rl0_to_l0 = data->$tmp$

nc3 = dimsizes(disth)
dims = dimsizes(l0rl0_to_l0)
nl0 = dims(0)
nl0r = dims(1)
if (nl0.gt.1) then
   il0_plot = toint(fspan(0,nl0-1,nl0_plot))
else
   il0_plot = (/0/)
end if

tmp = diag + "_raw"
raw = new((/nl0,nl0,nc3/),float)
raw = -999.0
raw@_FillValue = -999.0
do jl0r=0,nl0r-1
   do il0=0,nl0-1
      jl0 = l0rl0_to_l0(il0,jl0r)-1
      raw(il0,jl0,:) = tofloat(data->$tmp$(il0,jl0r,:))
      raw(jl0,il0,:) = raw(il0,jl0,:)
   end do
end do
maxval = max(abs(raw))

tmp = diag + "_fit"
isfit = isfilevar(data,tmp)
if (isfit) then
   fit = new((/nl0,nl0,nc3/),float)
   fit = -999.0
   fit@_FillValue = -999.0
   do jl0r=0,nl0r-1
      do il0=0,nl0-1
         jl0 = l0rl0_to_l0(il0,jl0r)-1
         fit(il0,jl0,:) = tofloat(data->$tmp$(il0,jl0r,:))
         fit(jl0,il0,:) = fit(il0,jl0,:)
      end do
   end do
   maxval = max((/maxval,max(abs(fit))/))

   tmp = diag + "_fit_rh"
   fit_rh = tofloat(data->$tmp$)
   tmp = diag + "_fit_rv"
   fit_rv = tofloat(data->$tmp$)
end if

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = maxval
res@cnMinLevelValF = -maxval
res@cnLevelSpacingF = 0.05*maxval

plot = new(2,graphic)

res@sfXArray = disth
res@sfYArray = vunit

do il0=0,nl0_plot-1
   output = "../fig/" + prefix + "_" + diag + "_lev" + sprinti("%0.3i",il0_plot(il0)+1)

   wks = gsn_open_wks("epsi",output)
   gsn_define_colormap(wks,"BlWhRe")

   res@tiMainString = "Raw diagnostic @ lev " + (il0_plot(il0)+1)
   plot(0) = gsn_csm_contour(wks,raw(:,il0_plot(il0),:),res)
   if (isfit) then
      res@tiMainString = "Fit @ lev " + (il0_plot(il0)+1)
      plot(1) = gsn_csm_contour(wks,fit(:,il0_plot(il0),:),res)
   end if

   gsn_panel(wks,plot,(/1,2/),resP)
   frame(wks)
   delete(wks)
   system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
end do

delete(res@sfXArray)
res@sfXArray = vunit

output = "../fig/" + prefix + "_" + diag + "_zs"

wks = gsn_open_wks("epsi",output)
gsn_define_colormap(wks,"BlWhRe")

res@tiMainString = "Raw diagnostic"
plot(0) = gsn_csm_contour(wks,raw(:,:,0),res)
if (isfit) then
   res@tiMainString = "Fit"
   plot(1) = gsn_csm_contour(wks,fit(:,:,0),res)
end if

gsn_panel(wks,plot,(/1,2/),resP)
frame(wks)
delete(wks)
system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")

if (isfit) then
   output = "../fig/" + prefix + "_" + diag + "_rh-rv"

   wks = gsn_open_wks("epsi",output)
   gsn_define_colormap(wks,"BlWhRe")

   xyres@trXMinF = 0.0
   xyres@trYMinF = min(vunit)
   xyres@trYMaxF = max(vunit)

   xyres@trXMaxF = max(fit_rh)*1.05
   xyres@tiMainString = "Horizontal support radius"
   plot(0) = gsn_csm_xy(wks,fit_rh,vunit,xyres)

   xyres@trXMaxF = max(fit_rv)*1.05
   xyres@tiMainString = "Vertical support radius"
   plot(1) = gsn_csm_xy(wks,fit_rv,vunit,xyres)

   gsn_panel(wks,plot,(/1,2/),resP)
   frame(wks)
   delete(wks)
   system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
end if

end
