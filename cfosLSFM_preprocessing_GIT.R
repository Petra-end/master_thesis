#IMPORT PACKAGES
library('reticulate')
library('abind')
library('freesurferformats')
library('tiff')
library('moments')
library(ggrepel)
np = import("numpy")

#CREATE THE CHANNEL
#open data as binary 5 channels
ni_dat = readTIFF('/beegfs/home/pmatyskova/project/nicolas_newdata/p_val_colors_fed_Refed.tif',
                  all=TRUE,as.is = TRUE,native = FALSE)

#separate the rgb channels
red_channel <- array(rep(0, 270*270*320), dim=c(270, 270, 320))
blue_channel <- array(rep(0, 270*270*320), dim=c(270, 270, 320))
green_channel <- array(rep(0, 270*270*320), dim=c(270, 270, 320))
for (i in 1:320) {
  red_channel[,,i] = ni_dat[[i]][,,1]
  green_channel[,,i] = ni_dat[[i]][,,2]
  blue_channel[,,i] = ni_dat[[i]][,,3]
}

#set to binary
red_channel[red_channel>0] = 1
green_channel[green_channel>0] = 1
blue_channel[blue_channel>0] = 1

#create cy channels
yellow_channel = red_channel * green_channel
cyan_channel = blue_channel * green_channel

#set rgb channels to 0 if value in cy channel already
red_channel[yellow_channel>0] = 0
green_channel[yellow_channel>0] = 0
green_channel[cyan_channel>0] = 0
blue_channel[cyan_channel>0] = 1

#to combine red & yellow and green & blue
ry_channel = red_channel + yellow_channel
ry_channel[ry_channel>0] = 1

gb_channel = blue_channel + green_channel
gb_channel[gb_channel>0] = 1

#rotate
red_channel = rotate3D(red_channel, axis = 1L, degrees = 270L)
green_channel = rotate3D(green_channel, axis = 1L, degrees = 270L)
cyan_channel = rotate3D(cyan_channel, axis = 1L, degrees = 270L)
yellow_channel = rotate3D(yellow_channel, axis = 1L, degrees = 270L)

ry_channel = rotate3D(ry_channel, axis = 1L, degrees = 270L)
gb_channel = rotate3D(gb_channel, axis = 1L, degrees = 270L)

#PAD WITH 0s
V1before <- array(rep(0, 130*320*270), dim=c(130, 320, 270))
V1after <- array(rep(0, 128*320*270), dim=c(128, 320, 270))
red_h = abind(V1before, ry_channel, along = 1)
red_h = abind(red_h, V1after, along = 1)
dim(red_h)

V3before = array(rep(0, 528*320*90), dim=c(528, 320, 90))
V3after = array(rep(0, 528*320*96), dim=c(528, 320, 96))
red_f = abind(V3before, red_h, along = 3)
red_f = abind(red_f, V3after, along = 3)
dim(red_f)
np$save("/beegfs/home/pmatyskova/project/nicnew_fedref_rychannelpad.npy", r_to_py(red_f))