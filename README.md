1. Two Fortran scripts for calculating CFODD are in ./src. 
   Input data are binary files with direct access (per record size: lonxlat).
   Run either of them by typing ./run.sh.
2. The Fortran results will be dumped into ./CFODD_output/...
3. An NCL script to plot CFODD from the above outputs is in ./plot.
4. For MIROC5 model (T85), a land_sea mask file is kept in ./utils. 

