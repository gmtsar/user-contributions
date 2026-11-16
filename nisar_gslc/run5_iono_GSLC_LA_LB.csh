#!/bin/csh -f
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
# GSLC run5: estimate ionosphere; see WORKFLOW_GSLC.md.
set script_dir = "$0:h"
if ("$script_dir" == "$0") set script_dir = "."
exec "$script_dir/gslc_workflow.csh" iono $argv:q
