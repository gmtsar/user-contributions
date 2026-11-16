#!/bin/csh -f
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
# Combined GSLC run1-6 entry; preserve each stage's full log.
set script_dir = "$0:h"
if ("$script_dir" == "$0") set script_dir = "."
if ("$script_dir" !~ /*) set script_dir = "$cwd/$script_dir"
set project_dir = "$cwd"
set pdf_args = ()
set first_step = 1
set reuse_run4 = 0

while ($#argv > 0)
    switch ("$argv[1]")
    case --project:
        if ($#argv < 2) then
            echo "--project requires a project directory."
            exit 2
        endif
        set project_dir = "$argv[2]"
        shift argv
        shift argv
        breaksw
    case --no-pdf:
        set pdf_args = (--no-pdf)
        shift argv
        breaksw
    case --from-step:
        if ($#argv < 2) then
            echo "--from-step requires a number from 1 to 6."
            exit 2
        endif
        if ("$argv[2]" !~ [1-6]) then
            echo "--from-step must be a number from 1 to 6."
            exit 2
        endif
        set first_step = "$argv[2]"
        shift argv
        shift argv
        breaksw
    case --resume:
        set reuse_run4 = 1
        shift argv
        breaksw
    case -h:
    case --help:
        echo "Usage: $0 [--project DIRECTORY] [--no-pdf] [--from-step 1..6] [--resume]"
        echo "Read two GSLC .h5 inputs from project/raw/ and run stages 1 through 6."
        echo "Default project: current directory. Default output: GRD and GMT PDF."
        echo "Run4-6 also export longitude/latitude GRD and plots; LOS uses direct wavelength conversion."
        echo "Logs: project/logs/all_UTC_TIMESTAMP_PID/run1.log through run6.log."
        echo "Each stage must succeed before the next starts; existing outputs are protected."
        echo "--from-step 4 --resume: retain completed run4 unwraps, retry geographic exports, then run5 and run6."
        exit 0
    default:
        echo "Unknown argument: $argv[1]; use --help."
        exit 2
    endsw
end

set stamp = `date -u +%Y%m%dT%H%M%SZ`
set log_dir = "$project_dir/logs/all_${stamp}_$$"
mkdir -p "$log_dir"
if ($status != 0) then
    echo "Cannot create log directory: $log_dir"
    exit 1
endif

set steps = (run1_prepare_GSLC_LA_LB.csh run2_process_GSLC_LA_LB.csh run3_presnaphu_GSLC_LA_LB.csh run4_unwrap_GSLC_LA_LB.csh run5_iono_GSLC_LA_LB.csh run6_corrected_GSLC_LA.csh)
echo "GSLC run${first_step} -> run6; project: $project_dir"
echo "Stage logs: $log_dir"
@ step = 1
foreach script ($steps:q)
    if ($step < $first_step) then
        @ step++
        continue
    endif
    set stage_log = "$log_dir/run${step}.log"
    echo "START run${step}: $script; log: $stage_log"
    set stage_args = ($pdf_args:q)
    if ($step == 4 && $reuse_run4 == 1) set stage_args = ($stage_args:q --resume)
    "$script_dir/$script" --project "$project_dir" $stage_args:q >& "$stage_log"
    set stage_status = $status
    if ($stage_status != 0) then
        echo "FAILED run${step} (exit $stage_status). See $stage_log"
        exit $stage_status
    endif
    echo "FINISHED run${step}"
    @ step++
end
echo "COMPLETE run${first_step} -> run6; corrected A: $project_dir/LA_iono_corrected/"
echo "Corrected LOS: $project_dir/LA_iono_corrected/los.grd and lonlat/los_ll.grd"
