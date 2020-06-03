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

procs = fspan(1,nproc,nproc)

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
res@trXMinF = 0.5
res@trXMaxF = nproc+0.5
res@trYMinF = 0.5
res@trYMaxF = nproc+0.5
res@tmXBMode = "Explicit"
res@tmYLMode = "Explicit"
res@tmXBValues = procs
res@tmYLValues = procs
res@tmXBLabels = procs
res@tmYLLabels = procs
res@sfXArray = procs
res@sfYArray = procs

resP = True
resP@gsnFrame = False
resP@gsnPanelXWhiteSpacePercent = 8
resP@gsnPanelYWhiteSpacePercent = 5

; Load and plot data:
; -------------------

com_AB_halo = new((/nproc,nproc/),integer)
com_AB_excl = new((/nproc,nproc/),integer)
com_AC_halo = new((/nproc,nproc/),integer)
com_AC_excl = new((/nproc,nproc/),integer)
do iproc=1,nproc
   data_nicas = addfile(datadir + "/" + prefix + "_nicas_" + mpicom + "_" + sprinti("%0.4i",nproc) + "-" + sprinti("%0.4i",iproc) + "_common.nc","r")
   com_AB_halo(iproc-1,:) = data_nicas->com_AB_jhalocounts
   com_AB_excl(iproc-1,:) = data_nicas->com_AB_jexclcounts
   com_AC_halo(iproc-1,:) = data_nicas->com_AC_jhalocounts
   com_AC_excl(iproc-1,:) = data_nicas->com_AC_jexclcounts
end do

output = "../fig/" + prefix + "_halos_" + mpicom + "_" + sprinti("%0.4i",nproc)
wks = gsn_open_wks("epsi",output)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

plot = new(4,graphic)

res@tiMainString = "Halos A/B, halo"
plot(0) = gsn_csm_contour(wks,com_AB_halo,res)
res@tiMainString = "Halos A/B, excl"
plot(1) = gsn_csm_contour(wks,com_AB_excl,res)
res@tiMainString = "Halos A/C, halo"
plot(2) = gsn_csm_contour(wks,com_AC_halo,res)
res@tiMainString = "Halos A/C, excl"
plot(3) = gsn_csm_contour(wks,com_AC_excl,res)

gsn_panel(wks,plot,(/2,2/),resP)
frame(wks)
delete(wks)
system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")

end
