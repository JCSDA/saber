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
nc_plot = 30
ic_colors = toint(fspan(2,201,nc_plot))

; Resources:
; ----------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

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

plres = True
plres@gsLineThicknessF = 3.0

; Load and plot data:
; -------------------

data = addfile(datadir + "/" + prefix + "_sampling.nc","r")
nl = data@nl
levs = new(nl,integer)
do il=0,nl-1
   levs(il) = toint(str_get_field(data@levs,il+1,":"))
end do
nc1 = data@nc1
nl0 = nl
nc3 = data@nc3
ic1_plot = toint(random_uniform(0,nc1-1,nc_plot))
lat_sampling = data->lat(:,:,ic1_plot)
lon_sampling = data->lon(:,:,ic1_plot)
lat_c1 = data->lat(0,0,:)
lon_c1 = data->lon(0,0,:)
local_diag = data@local_diag
displ_diag = data@displ_diag
if (local_diag.eq.".true.") then
   c2_to_c1 = data->c2_to_c1
   nc2 = dimsizes(c2_to_c1)
end if
dot_c1 = new(nc1,graphic)
if (local_diag.eq.".true.") then
   dot_c2 = new(nc2,graphic)
end if

output = "../fig/" + prefix + "_hdiag_c1"
wks = gsn_open_wks("epsi",output)
gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

res@tiMainString = "Samplings C1 and C2"
plot = gsn_csm_map(wks,res)
if (local_diag.eq.".true.") then
   pmres@gsMarkerColor = 170
   pmres@gsMarkerSizeF = 0.01
   do ic2=0,nc2-1
      ic1 = c2_to_c1(ic2)-1
      dot_c2(ic2) = gsn_add_polymarker(wks,plot,lon_c1(ic1),lat_c1(ic1),pmres)
   end do
end if
pmres@gsMarkerColor = 2
pmres@gsMarkerSizeF = 0.005
do ic1=0,nc1-1
   dot_c1(ic1) = gsn_add_polymarker(wks,plot,lon_c1(ic1),lat_c1(ic1),pmres)
end do

draw(plot)
frame(wks)
delete(wks)
system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")

dot = new((/nc_plot,nc3/),graphic)
lines = new((/nc_plot,nc3-1/),graphic)

do il0=0,nl0-1
   output = "../fig/" + prefix + "_hdiag_sampling_lev" + sprinti("%0.3i",(il0))
   wks = gsn_open_wks("epsi",output)
   gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

   res@tiMainString = "Sampling @ lev " + levs(il0)
   plot = gsn_csm_map(wks,res)
   do ic1=0,nc_plot-1
      pmres@gsMarkerColor = ic_colors(ic1)
      plres@gsLineColor = ic_colors(ic1)
      do ic3=0,nc3-1
         dot(ic1,ic3) = gsn_add_polymarker(wks,plot,lon_sampling(il0,ic3,ic1),lat_sampling(il0,ic3,ic1),pmres)
         if (ic3.gt.0) then
            lines(ic1,ic3-1) = gsn_add_polyline(wks,plot,(/lon_sampling(il0,0,ic1),lon_sampling(il0,ic3,ic1)/),(/lat_sampling(il0,0,ic1),lat_sampling(il0,ic3,ic1)/),plres)
         end if
      end do
   end do

   draw(plot)
   frame(wks)
   delete(wks)
   system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
end do

if ((local_diag.eq."true").or.(displ_diag.eq."true")) then
   do il=0,nl-1
      if (isfilepresent(datadir + "/" + prefix + "_sampling_" + sprinti("%0.3i",(il)) + ".nc")) then
         data = addfile(datadir + "/" + prefix + "_sampling_" + sprinti("%0.3i",(il)) + ".nc","r")
         ic2_plot = toint(random_uniform(0,nc2-1,nc_plot))
         local_mask = data->local_mask(ic2_plot,:)
         displ_mask = data->local_mask(ic2_plot,:)

         if (local_diag.eq."true") then
            output = "../fig/" + prefix + "_hdiag_local_mask"
            wks = gsn_open_wks("epsi",output)
            gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

            res@tiMainString = "Local mask @ independent level " + il
            plot = gsn_csm_map(wks,res)
            pmres@gsMarkerSizeF = 0.01
            do ic2=0,nc_plot-1
               pmres@gsMarkerColor = ic_colors(ic2)
               do ic1=0,nc1-1
                  if (local_mask(ic2,ic1).eq.1) then
                     dot_c1(ic1) = gsn_add_polymarker(wks,plot,lon_c1(ic1),lat_c1(ic1),pmres)
                  end if
               end do
            end do

            draw(plot)
            frame(wks)
            delete(wks)
            system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
         end if

         if (displ_diag.eq."true") then
            output = "../fig/" + prefix + "_hdiag_displ_mask"
            wks = gsn_open_wks("epsi",output)
            gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

            res@tiMainString = "Displacement mask @ independent level " + il
            plot = gsn_csm_map(wks,res)
            pmres@gsMarkerSizeF = 0.01
            do ic2=0,nc_plot-1
               pmres@gsMarkerColor = ic_colors(ic2)
               do ic1=0,nc1-1
                  if (displ_mask(ic2,ic1).eq.1) then
                     dot_c1(ic1) = gsn_add_polymarker(wks,plot,lon_c1(ic1),lat_c1(ic1),pmres)
                  end if
               end do
            end do

            draw(plot)
            frame(wks)
            delete(wks)
            system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")
         end if
      end if
   end do
end if

end
