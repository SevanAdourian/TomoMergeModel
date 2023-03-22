#!/bin/bash

depth=0050
input_file='../ascii/semucb/global_50.xyz'
output_file='./global_'$depth'km.ps'
range_global="-180/180/-90/90"
# range_global="0/360/-90/90"
# range_plot="-R160/270/30/85"
# range_regional="-150/-100/20/50"
# frame_plot="W20c"
frame_plot="A110/25/15c"
bor=90/45

# cat<<EOF>$dir_home"mask.xyz"
# 235 31.5
# 235 43
# 246 43
# 246 31.5
# 235 31.5
# EOF

cat<<EOF>$dir_home"mask.xyz"
-125 31.5
-125 43
-114 43
-114 31.5
-125 31.5
EOF

cat<<EOF>"./panoply_white.cpt"
# GMT palette Blues_09.cpt
#
# This product includes color specifications and designs
# developed by Cynthia Brewer (http://colorbrewer.org/).
#
# Converted to the cpt format by J.J.Green
# Sequential palette with 9 colours
#
# COLOR_MODEL = RGB
0       4 14 216        0.0625  4 14 216
0.0625  32 80 255       0.125   32 80 255
0.125   65 150 255      0.1875  65 150 255
0.1875  109 193 255     0.25    109 193 255
0.25    134 217 255     0.3125  134 217 255
0.3125  156 238 255     0.375   156 238 255
0.375   175 245 255     0.4375  175 245 255
0.4375  206 255 255     0.45    206 255 255
0.45    233 255 255     0.475   244 255 255
0.475   255 255 255     0.525   255 255 255
0.525   255 255 200     0.55    255 255 200
0.55    255 254 71      0.5625  255 254 71
0.5625  255 235 0       0.625   255 235 0
0.625   255 196 0       0.6875  255 196 0
0.6875  255 144 0       0.75    255 144 0
0.75    255 72 0        0.8125  255 72 0
0.8125  red     0.875   red
0.875   213 0 0 0.9375  213 0 0
0.9375  158 0 0 1       158 0 0
B       0 0 121
F       87 0 0
N       200
EOF

gmt makecpt -I -T-0.1/0.1/0.001 -D -Z -C"panoply_white.cpt" > "tomo.cpt"

# gmt surface -R$range_global -H1 -I0.5 $input_file -Goutgrid.grd
gmt surface -R$range_global -I0.5 $input_file -Goutgrid.grd

gmt grdimage -R -J$frame_plot outgrid.grd -Bp20WSen:."depth ${depth} km": -C$dir_home"tomo.cpt" -Y5c -X3c -K > $output_file
gmt pscoast -R -J -Dl -W0.03c,black -A5000 -O -K >> $output_file
gmt psxy -R -J "./mask.xyz" -W0.03c,purple,-. -O -K >> $output_file
gmt psscale -D10/-1c/5c/0.25ch -C"./tomo.cpt" -B0.1:"dlnVs (%)": -O >> $output_file

rm outgrid.grd tomo.cpt panoply_white.cpt gmt.history mask.xyz
