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
nproc = 3
mpicom = 2
strategy = "common"

; Resources:
; ----------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True
res@gsnAddCyclic = False

res@cnFillOn = True
res@cnFillMode = "RasterFill"
res@cnMonoFillPattern = True
res@cnMonoFillColor = False
res@cnInfoLabelOn = False
res@cnLineLabelsOn = False
res@cnLinesOn = False
res@cnNoDataLabelOn = False
res@cnMissingValFillColor = 0
res@lbLabelBarOn = False
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

pmres = True
pmres@gsMarkerIndex = 1
pmres@gsMarkerSizeF = 0.005
pmres@gsMarkerColor = "black"

resP = True
resP@gsnFrame = False
resP@gsnPanelXWhiteSpacePercent = 8
resP@gsnPanelYWhiteSpacePercent = 5

; Load and plot data:
; -------------------

do iproc=1,nproc
   if (strategy.eq."common") then
      ; Load data
      data_nicas = addfile(datadir + "/" + prefix + "_nicas_" + mpicom + "_" + sprinti("%0.4i",nproc) + "-" + sprinti("%0.4i",iproc) + "_common_summary.nc","r")
      if (iproc.eq.1) then
         nl = data_nicas@nl
         levs = new(nl,integer)
         do il=0,nl-1
            levs(il) = toint(str_get_field(data_nicas@levs,il+1,":"))
         end do
         lon = data_nicas->lon
         lat = data_nicas->lat
         c1_to_c0 = data_nicas->c1_to_c0
         l1_to_l0 = data_nicas->l1_to_l0
         nc1 = dimsizes(c1_to_c0)
         nl1 = dimsizes(l1_to_l0)
      end if

      lcheck = data_nicas->lcheck
      lcheck@_FillValue = -999.0

      do il1=0,nl1-1
         ; Subgrid at a given level
         il0 = l1_to_l0(il1)
         nc2 = num(.not.ismissing(lcheck(il1,:)))
         lon_plot = new(nc2,double)
         lat_plot = new(nc2,double)
         lcheck_plot = new(nc2,double)
         ic2 = 0
         do ic1=0,nc1-1
            if (.not.ismissing(lcheck(il1,ic1))) then
               lon_plot(ic2) = lon(c1_to_c0(ic1)-1)
               lat_plot(ic2) = lat(c1_to_c0(ic1)-1)
               lcheck_plot(ic2) = lcheck(il1,ic1)
               ic2 = ic2+1
            end if
         end do

         ; Plot data
         output = "../fig/" + prefix + "_halos_" + mpicom + "_" + sprinti("%0.4i",nproc) + "-" + sprinti("%0.4i",iproc) + "_lev" + sprinti("%0.3i",levs(il0-1))
         wks = gsn_open_wks("epsi",output)
         gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
         gsn_reverse_colormap(wks)
         dot = new(nc2,graphic)

         res@tiMainString = "Subgrid @ lev " + levs(il0-1) + " (" + nc2 + " points)"
         res@sfXArray = lon_plot
         res@sfYArray = lat_plot
         plot = gsn_csm_contour_map(wks,lcheck_plot,res)
         dot = gsn_add_polymarker(wks,plot,lon_plot,lat_plot,pmres)

         draw(plot)
         frame(wks)
         delete(wks)
         delete(dot)
         delete(lon_plot)
         delete(lat_plot)
         delete(lcheck_plot)
         delete(res@sfXArray)
         delete(res@sfYArray)
         system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
      end do
   end if
end do

end
