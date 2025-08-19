#!/bin/csh -f

if ($#argv != 1) then
  echo ""
  echo "grd180.csh grid file"
  echo ""
  echo "    shift grid from 0 360 to -180 180"
  echo ""
  exit 1
endif

set lon1 = `gmt grdinfo $1 -C | awk '{printf("%.12f", $2-360.0)}'`
set lon2 = `gmt grdinfo $1 -C | awk '{printf("%.12f", $3-360.0)}'`
set lat1 = `gmt grdinfo $1 -C | awk '{print $4}'`
set lat2 = `gmt grdinfo $1 -C | awk '{print $5}'`

gmt grdedit $1 -R$lon1/$lon2/$lat1/$lat2

