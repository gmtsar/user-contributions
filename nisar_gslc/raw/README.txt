Place exactly two NISAR GSLC .h5 files (two acquisition dates) in this folder.
Symbolic links to the original .h5 files are also supported.
Each input HDF5 contains both frequency A and frequency B; do not split them.

From the project directory, run:
  ./run1_prepare_GSLC_LA_LB.csh

To execute all six stages sequentially from the project directory, run:
  ./run1-6_gslc_all.csh

Run1 reads raw/*.h5, sorts the pair by acquisition metadata, and creates the
LA/LB configurations and links. Input HDF5 files are opened read-only.
When using --project /another/project, place the two inputs in that project's
raw/ directory instead. This README is ignored during input discovery.
