#!/usr/bin/python

import sys
import Image
import ImageEnhance

#fbase="/media/1A8F-0227/Landsat/LT50150282011282EDC00/LT50150282011282EDC00"
#fbase="/cygdrive/d/Landsat/LT50190272011278EDC00/LT50190272011278EDC00"
#fbase="/media/Acer/cygwin/home/Petey/data/landsat/LE70160282011169EDC00"
#outfile="image2.pgm"

fbase=sys.argv[1]
outfile=sys.argv[2]

factor=1

fblue=fbase+"_B1.TIF"
fgreen=fbase+"_B2.TIF"
fred=fbase+"_B3.TIF"

iblue=Image.open(fblue)
contr=ImageEnhance.Contrast(iblue)
iblue=contr.enhance(factor)
igreen=Image.open(fgreen)
contr=ImageEnhance.Contrast(igreen)
igreen=contr.enhance(factor)
ired=Image.open(fred)
contr=ImageEnhance.Contrast(ired)
ired=contr.enhance(factor)

irgb=Image.merge("RGB", (ired, igreen, iblue))

#color=ImageEnhance.Color(irgb)
#final=color.enhance(5)

#final.save(outfile)
irgb.save(outfile)

