#!/bin/csh -f
#
#  pot_merge.csh  -  Merge north/south geocoded offset grids
#
#  1. Computes bias from overlap via XYZ (robust, no grid alignment issues)
#  2. Applies bias correction to south
#  3. Combines via blockmedian (overlap = median of bias-corrected values)
#
#  Usage: pot_merge.csh north_dir south_dir output_dir
#
# ============================================================================

if ($#argv < 3) then
  echo "Usage: pot_merge.csh north_dir south_dir output_dir"
  exit 1
endif

set north_dir  = $1
set south_dir  = $2
set output_dir = $3

mkdir -p $output_dir
set tmpdir = $output_dir/tmp_merge
mkdir -p $tmpdir

echo "============================================"
echo " pot_merge.csh"
echo "============================================"

# ---- Grid parameters ----
set info_n = `gmt grdinfo $north_dir/azi_offset_ll.grd -C`
set info_s = `gmt grdinfo $south_dir/azi_offset_ll.grd -C`

set n_latmin = `echo $info_n | awk '{print $4}'`
set n_latmax = `echo $info_n | awk '{print $5}'`
set s_latmin = `echo $info_s | awk '{print $4}'`
set s_latmax = `echo $info_s | awk '{print $5}'`
set dx = `echo $info_n | awk '{print $8}'`
set dy = `echo $info_n | awk '{print $9}'`

# Overlap
set ov_latmin = `echo $n_latmin $s_latmin | awk '{printf "%.10f", ($1>$2)?$1:$2}'`
set ov_latmax = `echo $n_latmax $s_latmax | awk '{printf "%.10f", ($1<$2)?$1:$2}'`
set ov_span = `echo $ov_latmax $ov_latmin | awk '{printf "%.4f",$1-$2}'`

echo " North lat: $n_latmin ~ $n_latmax"
echo " South lat: $s_latmin ~ $s_latmax"
echo " Overlap: $ov_latmin ~ $ov_latmax ($ov_span deg)"
echo " dx=$dx  dy=$dy"

# Clean merged -R (floor/ceil to grid nodes)
set n_lonmin = `echo $info_n | awk '{print $2}'`
set n_lonmax = `echo $info_n | awk '{print $3}'`
set s_lonmin = `echo $info_s | awk '{print $2}'`
set s_lonmax = `echo $info_s | awk '{print $3}'`

set m_xmin = `echo $n_lonmin $s_lonmin $dx | awk '{v=($1<$2)?$1:$2; printf "%.10f", int(v/$3)*$3}'`
set m_xmax = `echo $n_lonmax $s_lonmax $dx | awk '{v=($1>$2)?$1:$2; printf "%.10f", (int(v/$3)+1)*$3}'`
set m_ymin = `echo $n_latmin $s_latmin $dy | awk '{v=($1<$2)?$1:$2; printf "%.10f", int(v/$3)*$3}'`
set m_ymax = `echo $n_latmax $s_latmax $dy | awk '{v=($1>$2)?$1:$2; printf "%.10f", (int(v/$3)+1)*$3}'`
set merged_R = "-R$m_xmin/$m_xmax/$m_ymin/$m_ymax"
echo " Merged: $merged_R"
echo ""

# ============================================================================
foreach comp (azi_offset rng_offset snr)

  set n_grd = $north_dir/${comp}_ll.grd
  set s_grd = $south_dir/${comp}_ll.grd
  set o_grd = $output_dir/${comp}_merged.grd

  if (! -f $n_grd || ! -f $s_grd) then
    echo "  Skipping $comp (not found)"
    continue
  endif

  echo "---- $comp ----"

  # ---- Bias from overlap (XYZ-based, fast) ----
  if ($comp != "snr") then
    gmt grd2xyz $n_grd -s | awk -v y0=$ov_latmin -v y1=$ov_latmax '{if($2>=y0 && $2<=y1) print $3}' > $tmpdir/n_ov.txt
    gmt grd2xyz $s_grd -s | awk -v y0=$ov_latmin -v y1=$ov_latmax '{if($2>=y0 && $2<=y1) print $3}' > $tmpdir/s_ov.txt
    set n_med = `sort -g $tmpdir/n_ov.txt | awk '{a[NR]=$1}END{print a[int(NR/2)+1]}'`
    set s_med = `sort -g $tmpdir/s_ov.txt | awk '{a[NR]=$1}END{print a[int(NR/2)+1]}'`
    set bias = `echo $n_med $s_med | awk '{printf "%.6f",$1-$2}'`
    echo "  Bias (N-S): $bias  (N median=$n_med, S median=$s_med)"
  else
    set bias = "0"
  endif

  # ---- Combine XYZ: north as-is + south with bias correction ----
  echo "  Extracting and combining XYZ ..."
  gmt grd2xyz $n_grd -s > $tmpdir/combined.xyz
  gmt grd2xyz $s_grd -s | awk -v b=$bias '{printf "%.10f %.10f %.10f\n", $1, $2, $3+b}' >> $tmpdir/combined.xyz

  set npts = `wc -l < $tmpdir/combined.xyz`
  echo "  Combined: $npts points"

  # ---- blockmedian + xyz2grd ----
  echo "  Gridding ..."
  gmt blockmedian $tmpdir/combined.xyz $merged_R -I$dx/$dy -r > $tmpdir/median.xyz
  set ngrid = `wc -l < $tmpdir/median.xyz`
  echo "  Grid cells: $ngrid"
  gmt xyz2grd $tmpdir/median.xyz $merged_R -I$dx/$dy -r -G$o_grd

  echo "  -> $o_grd"
  gmt grdinfo $o_grd -C | awk '{printf "     lon %.4f~%.4f  lat %.4f~%.4f  z %.2f~%.2f\n",$2,$3,$4,$5,$6,$7}'
  echo ""

end

# ---- Cleanup ----
rm -rf $tmpdir

echo "============================================"
echo " Done. Output in $output_dir/"
ls -lh $output_dir/*.grd
echo "============================================"
