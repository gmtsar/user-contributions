#!/usr/bin/env python3
"""Compare complete xcorr output files, retaining every differing row.

Exit zero if coordinates/order and accepted-point tolerances pass. Byte equality
is reported separately. This does not certify performance or downstream fitting.
Example: python3 test/compare_outputs.py stock.dat parallel.dat --rows 1000
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

def compare(reference, candidate, rows, threshold=18.0, offset_tolerance=0.0625,
            correlation_tolerance=0.01):
    a=np.loadtxt(reference,ndmin=2);b=np.loadtxt(candidate,ndmin=2)
    if a.shape!=(rows,5) or b.shape!=(rows,5):
        raise ValueError(f'Expected {rows} rows and 5 columns; got {a.shape}, {b.shape}')
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Nonfinite output values')
    coordinates=bool(np.array_equal(a[:,[0,2]],b[:,[0,2]]))
    ma=a[:,4]>=threshold;mb=b[:,4]>=threshold;accepted=ma|mb
    delta=np.abs(a-b);changed=np.any(delta!=0,axis=1)
    offsets=float(delta[accepted][:,[1,3]].max(initial=0))
    correlation=float(delta[accepted,4].max(initial=0))
    mask_equal=bool(np.array_equal(ma,mb))
    result={
        'reference_sha256':hashlib.sha256(Path(reference).read_bytes()).hexdigest(),
        'candidate_sha256':hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),
        'rows':rows,'identical_coordinates_and_order':coordinates,
        'accepted_reference':int(ma.sum()),'accepted_candidate':int(mb.sum()),
        'identical_acceptance_mask':mask_equal,'threshold':threshold,
        'different_rows':int(changed.sum()),
        'different_accepted_rows':int((changed&accepted).sum()),
        'max_accepted_offset_difference_px':offsets,
        'max_accepted_correlation_difference':correlation,
        'max_all_offset_difference_px':float(delta[:,[1,3]].max(initial=0)),
        'pass':bool(coordinates and mask_equal and offsets<=offset_tolerance+1e-8
                    and correlation<=correlation_tolerance+1e-8),
    }
    result['byte_equal']=result['reference_sha256']==result['candidate_sha256']
    differences=np.column_stack([np.flatnonzero(changed)+1,a[changed],b[changed]])
    return result,differences

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference',type=Path);p.add_argument('candidate',type=Path)
    p.add_argument('--rows',required=True,type=int)
    p.add_argument('--threshold',type=float,default=18.0)
    p.add_argument('--offset-tolerance',type=float,default=0.0625)
    p.add_argument('--correlation-tolerance',type=float,default=0.01)
    p.add_argument('--differences',type=Path)
    args=p.parse_args()
    result,differences=compare(args.reference,args.candidate,args.rows,args.threshold,
                               args.offset_tolerance,args.correlation_tolerance)
    if args.differences:
        np.savetxt(args.differences,differences,
                   header='row reference[x dx y dy corr] candidate[x dx y dy corr]')
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['pass'] else 1)
