#!/bin/csh -f
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
# Shared Python dispatcher for the six numbered GSLC entrypoints.
set script_dir = "$0:h"
if ("$script_dir" == "$0") set script_dir = "."
if ($?NISAR_GSLC_PYTHON) then
    set candidates = ("$NISAR_GSLC_PYTHON")
else
    set candidates = (python python3)
endif
set selected = ""
foreach candidate ($candidates:q)
    "$candidate" -c 'import numpy, h5py, scipy, netCDF4, pyproj' >& /dev/null
    if ($status == 0) then
        set selected = "$candidate"
        break
    endif
end
if ("$selected" == "") then
    echo "Set NISAR_GSLC_PYTHON to a Python with numpy/h5py/scipy/netCDF4/pyproj."
    exit 1
endif
if ($#argv == 0) then
    exec "$selected" "$script_dir/gslc_workflow.py" --help
endif
exec "$selected" "$script_dir/gslc_workflow.py" $argv:q
