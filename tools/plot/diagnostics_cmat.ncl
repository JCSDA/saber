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
var = (/"common"/)

nvar = dimsizes(var)

; Resources :
; -----------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True
res@gsnAddCyclic = False

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
xyres@tiYAxisString = "Level"

resP = True
resP@gsnFrame = False
resP@gsnPanelXWhiteSpacePercent = 8
resP@gsnPanelYWhiteSpacePercent = 5

; Load and plot data:
;--------------------

do ivar=0,nvar-1
   data = addfile(datadir + "/" + prefix + "_cmat_" + var + ".nc","r")
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
   displ_diag = (data@displ_diag.eq.".true.")
   if (ivar.eq.0) then
      lon = data->lon
      lat = data->lat
      res@sfXArray = lon
      res@sfYArray = lat
   end if
   coef_ens = tofloat(data->coef_ens)
   coef_sta = tofloat(data->coef_sta)
   rh0 = tofloat(data->rh0)*6371.229
   rv0 = tofloat(data->rv0)
   if (displ_diag) then
      displ_lon = tofloat(data->displ_lon)
      displ_lat = tofloat(data->displ_lat)
   end if

   ; Plot maps
   do ilevs=0,nlevs-1
      if (min(coef_ens(ilevs_plot(ilevs),:)).ne.max(coef_ens(ilevs_plot(ilevs),:))) then
         output = "../fig/" + prefix + "_coef_ens_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
         wks = gsn_open_wks("epsi",output)
         gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

         res@tiMainString = "Ensemble coefficient @ lev " + levs(ilevs_plot(ilevs))
         plot = gsn_csm_contour_map(wks,coef_ens(ilevs_plot(ilevs),:),res)

         draw(plot)
         frame(wks)
         delete(wks)
         system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
      end if

      if (min(coef_sta(ilevs_plot(ilevs),:)).ne.max(coef_sta(ilevs_plot(ilevs),:))) then
         output = "../fig/" + prefix + "_coef_sta_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
         wks = gsn_open_wks("epsi",output)
         gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

         res@tiMainString = "Static coefficient @ lev " + levs(ilevs_plot(ilevs))
         plot = gsn_csm_contour_map(wks,coef_sta(ilevs_plot(ilevs),:),res)

         draw(plot)
         frame(wks)
         delete(wks)
         system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
      end if

      if (min(rh0(ilevs_plot(ilevs),:)).ne.max(rh0(ilevs_plot(ilevs),:))) then
         output = "../fig/" + prefix + "_rh0_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
         wks = gsn_open_wks("epsi",output)
         gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

         res@tiMainString = "Horizontal support radius (km) @ lev " + levs(ilevs_plot(ilevs))
         plot = gsn_csm_contour_map(wks,rh0(ilevs_plot(ilevs),:),res)

         draw(plot)
         frame(wks)
         delete(wks)
         system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
      end if

      if (min(rv0(ilevs_plot(ilevs),:)).ne.max(rv0(ilevs_plot(ilevs),:))) then
         output = "../fig/" + prefix + "_rv0_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
         wks = gsn_open_wks("epsi",output)
         gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

         res@tiMainString = "Vertical support radius @ lev " + levs(ilevs_plot(ilevs))
         plot = gsn_csm_contour_map(wks,rv0(ilevs_plot(ilevs),:),res)

         draw(plot)
         frame(wks)
         delete(wks)
         system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
      end if

      if (displ_diag) then
         do its=1,nts-1
            output = "../fig/" + prefix + "_displ_lon_ts" + timeslot(its) + "_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
            wks = gsn_open_wks("epsi",output)
            gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

            res@tiMainString = "Displaced longitude @ lev " + levs(ilevs_plot(ilevs)) + " / ts " + timeslot(its)
            plot = gsn_csm_contour_map(wks,displ_lon(its-1,ilevs_plot(ilevs),:),res)

            draw(plot)
            frame(wks)
            delete(wks)
            system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")

            output = "../fig/" + prefix + "_displ_lat_ts" + timeslot(its) + "_lev" + sprinti("%0.3i",levs(ilevs_plot(ilevs)))
            wks = gsn_open_wks("epsi",output)
            gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

            res@tiMainString = "Displaced latitude @ lev " + levs(ilevs_plot(ilevs)) + " / ts " + timeslot(its)
            plot = gsn_csm_contour_map(wks,displ_lat(its-1,ilevs_plot(ilevs),:),res)

            draw(plot)
            frame(wks)
            delete(wks)
            system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
         end do
      end if
   end do

   if (nl.gt.1) then
      ; Plot profiles
      coef_ens_prof = new((/3,nl/),float)
      coef_sta_prof = new((/3,nl/),float)
      rh0_prof = new((/3,nl/),float)
      rv0_prof = new((/3,nl/),float)
      do il=0,nl-1
         lavg = avg(coef_ens(il,:))
         lstd = stddev(coef_ens(il,:))
         coef_ens_prof(:,il) = (/lavg,lavg-lstd,lavg+lstd/)
         lavg = avg(coef_sta(il,:))
         lstd = stddev(coef_sta(il,:))
         coef_sta_prof(:,il) = (/lavg,lavg-lstd,lavg+lstd/)
         lavg = avg(rh0(il,:))
         lstd = stddev(rh0(il,:))
         rh0_prof(:,il) = (/lavg,lavg-lstd,lavg+lstd/)
         lavg = avg(rv0(il,:))
         lstd = stddev(rv0(il,:))
         rv0_prof(:,il) = (/lavg,lavg-lstd,lavg+lstd/)
      end do

      output = "../fig/" + prefix + "_cmat_profiles"
      wks = gsn_open_wks("epsi",output)
      gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

      xyres@trYMinF = min(levs)
      xyres@trYMaxF = max(levs)
      plot_prof = new(4,graphic)

      xyres@trXMinF = -0.05
      xyres@trXMaxF = 1.05
      xyres@tiXAxisString = "Ensemble coefficient"
      plot_prof(0) = gsn_csm_xy(wks,coef_ens_prof,levs,xyres)

      xyres@tiXAxisString = "Static coefficient"
      plot_prof(1) = gsn_csm_xy(wks,coef_sta_prof,levs,xyres)

      xyres@trXMinF = 0.0
      xyres@trXMaxF = max(rh0_prof)*1.1
      xyres@tiXAxisString = "Horizontal length-scale (km)"
      plot_prof(2) = gsn_csm_xy(wks,rh0_prof,levs,xyres)

      xyres@trXMaxF = max(rv0_prof)*1.1
      xyres@tiXAxisString = "Vertical length-scale"
      plot_prof(3) = gsn_csm_xy(wks,rv0_prof,levs,xyres)

      gsn_panel(wks,plot_prof,(/2,2/),resP)
      frame(wks)
      delete(wks)
      system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
   end if
end do

end
