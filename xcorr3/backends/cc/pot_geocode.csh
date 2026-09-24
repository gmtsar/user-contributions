#!/bin/csh -f
#
#  pot_geocode.csh  -  Pixel Offset Tracking: freq_xcorr.dat → geocoded offset grids
#
#  Based on make_xcorr_plot.csh, with added:
#    - Search boundary filtering
#    - Configurable SNR threshold
#    - proj_ra2ll.csh geocoding
#
#  Usage:
#    pot_geocode.csh nx ny master.PRM trans.dat [snr] [xsearch] [ysearch] [ri]
#
#  nx/ny must match the values used in xcorr_cc / xcorr
#
# ============================================================================

if ($#argv < 4) then
  echo ""
  echo "Usage: pot_geocode.csh nx ny master.PRM trans.dat [snr] [xsearch] [ysearch] [ri]"
  echo ""
  echo "  nx ny       - MUST match the -nx -ny used in xcorr_cc"
  echo "  master.PRM  - reference image parameter file"
  echo "  trans.dat   - coordinate transformation table (from topo/)"
  echo "  snr         - correlation threshold (default: 10)"
  echo "  xsearch     - search half-window range (default: 64)"
  echo "  ysearch     - search half-window azimuth (default: 64)"
  echo "  ri          - range interpolation factor (default: 2)"
  echo ""
  exit 1
endif

set nx      = $1
set ny      = $2
set master  = $3
set trans   = $4
set SNR     = 10
set xsearch = 64
set ysearch = 64
set ri      = 2
set PSNR    = 5
if ($#argv >= 5) set SNR     = $5
if ($#argv >= 6) set xsearch = $6
if ($#argv >= 7) set ysearch = $7
if ($#argv >= 8) set ri      = $8
if ($#argv >= 9) set PSNR    = $9

if (! -f freq_xcorr.dat) then
  echo "ERROR: freq_xcorr.dat not found in current directory"
  exit 1
endif

# ---- Compute pixel sizes from PRM (same method as original make_xcorr_plot.csh) ----
set PRF            = `grep PRF $master            | awk -F"=" '{print $2}'`
set SC_vel         = `grep SC_vel $master         | awk -F"=" '{print $2}'`
set earth_radius   = `grep earth_radius $master   | awk -F"=" '{print $2}'`
set SC_height      = `grep SC_height $master      | awk -F"=" '{print $2}'`
set rng_samp_rate  = `grep rng_samp_rate $master  | head -1 | awk '{print $3}'`

set ground_vel = `echo $SC_vel $earth_radius $SC_height | awk '{print $1/sqrt(1+$3/$2)}'`
set azi_size   = `echo $ground_vel $PRF | awk '{printf "%.6f", $1/$2}'`
set rng_size   = `echo $rng_samp_rate   | awk '{printf "%.6f", 299792458.0/$1/2}'`

# ---- Compute search boundary thresholds ----
set max_rng = `echo $xsearch $ri | awk '{printf "%.1f", $1/$2 - 2}'`
set max_azi = `echo $ysearch     | awk '{printf "%.1f", $1 - 2}'`

echo "============================================"
echo " pot_geocode.csh"
echo "============================================"
echo " azi pixel size:  $azi_size m"
echo " rng pixel size:  $rng_size m"
echo " nx=$nx  ny=$ny"
echo " SNR threshold:   $SNR"
echo " peak_snr thr:    $PSNR (col 6; skipped for legacy 5-col files)"
echo " xsearch=$xsearch  ysearch=$ysearch  ri=$ri"
echo " Max allowed offset: |rng| < $max_rng px, |azi| < $max_azi px"
echo "============================================"

# ---- Step 1: Filter freq_xcorr.dat ----
# peak_snr (col 6) rejects spurious/gross peaks; (NF < 6) keeps legacy 5-col files working.
echo "Filtering ..."
set nraw = `wc -l < freq_xcorr.dat`

awk '{ \
  if ($5 > '$SNR' && (NF < 6 || $6 > '$PSNR') && $2 > -'$max_rng' && $2 < '$max_rng' && $4 > -'$max_azi' && $4 < '$max_azi') \
    print $0 \
}' freq_xcorr.dat > pot_filtered.dat

set nfilt = `wc -l < pot_filtered.dat`
echo " Passed filter: $nfilt / $nraw"

if ($nfilt == 0) then
  echo "ERROR: no points passed the filter."
  exit 1
endif

# ---- Step 2: Extract azi / rng components (same column order as original) ----
awk '{print $1, $3, $4, $5}' pot_filtered.dat > azi.dat
awk '{print $1, $3, $2, $5}' pot_filtered.dat > rng.dat

# ---- Step 3: Grid parameters (same as original make_xcorr_plot.csh) ----
set xmin = `gmt gmtinfo azi.dat -C | awk '{print $1}'`
set xmax = `gmt gmtinfo azi.dat -C | awk '{print $2}'`
set ymin = `gmt gmtinfo azi.dat -C | awk '{print $3}'`
set ymax = `gmt gmtinfo azi.dat -C | awk '{print $4}'`

set xinc = `echo $xmax $xmin $nx | awk '{printf "%.12f", ($1-$2)/($3-1)}'`
set yinc = `echo $ymax $ymin $ny | awk '{printf "%.12f", ($1-$2)/($3-1)}'`

echo " Region: -R$xmin/$xmax/$ymin/$ymax"
echo " xinc: $xinc  yinc: $yinc"

# ---- Step 4: blockmedian + xyz2grd (identical to original) ----
echo "Gridding ..."
gmt blockmedian azi.dat -R$xmin/$xmax/$ymin/$ymax -I$xinc/$yinc -Wi | awk '{print $1, $2, $3}' > azi_b.dat
gmt blockmedian rng.dat -R$xmin/$xmax/$ymin/$ymax -I$xinc/$yinc -Wi | awk '{print $1, $2, $3}' > rng_b.dat

gmt xyz2grd azi_b.dat -R$xmin/$xmax/$ymin/$ymax -I$xinc/$yinc -Gaoff.grd
gmt xyz2grd rng_b.dat -R$xmin/$xmax/$ymin/$ymax -I$xinc/$yinc -Groff.grd

# ---- Step 5: Convert pixel offsets to meters (identical to original) ----
echo "Converting to meters ..."
gmt grdmath aoff.grd $azi_size MUL = azi_offset.grd
gmt grdmath roff.grd -$rng_size MUL = rng_offset.grd

# ---- Step 6: Geocode (radar coords → lat/lon) ----
echo "Geocoding azi_offset ..."
proj_ra2ll.csh $trans azi_offset.grd azi_offset_ll.grd
echo "Geocoding rng_offset ..."
proj_ra2ll.csh $trans rng_offset.grd rng_offset_ll.grd

# ---- Cleanup ----
rm -f azi.dat rng.dat azi_b.dat rng_b.dat pot_filtered.dat
rm -f rap llp llpb

# ---- Report ----
echo ""
echo "============================================"
echo " Output (radar coord): aoff.grd roff.grd azi_offset.grd rng_offset.grd"
echo " Output (geographic):  azi_offset_ll.grd rng_offset_ll.grd"
echo "============================================"
gmt grdinfo azi_offset_ll.grd -C | awk '{printf " azi: lon %.4f~%.4f lat %.4f~%.4f z %.2f~%.2f m\n",$2,$3,$4,$5,$6,$7}'
gmt grdinfo rng_offset_ll.grd -C | awk '{printf " rng: lon %.4f~%.4f lat %.4f~%.4f z %.2f~%.2f m\n",$2,$3,$4,$5,$6,$7}'
echo "Done."
