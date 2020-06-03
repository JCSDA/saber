load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Parameters:
; -----------

system("mkdir -p fig")
home = getenv("HOME")
datadir = home + "/data/bump/arp/877D"
prefix = "aro_7G0N"
nproc = 2

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
res@tiMainString = "Distribution for " + nproc + " tasks"

; Load and plot data:
; -------------------

data_distribution = addfile(datadir + "/" + prefix + "_distribution_" + sprinti("%0.4i",nproc) + ".nc","r")
mpicom = data_distribution@mpicom
strategy = data_distribution@strategy

lon = data_distribution->lon
lat = data_distribution->lat
c0_to_proc = tofloat(data_distribution->c0_to_proc)

output = "../fig/" + prefix + "_distribution_" + sprinti("%0.4i",nproc)
system("mkdir -p ../fig")
wks = gsn_open_wks("epsi",output)
gsn_define_colormap(wks,"BlWhRe")

res@sfXArray = lon
res@sfYArray = lat
plot = gsn_csm_contour_map(wks,c0_to_proc,res)

draw(plot)
frame(wks)
delete(wks)
system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")

end
