# required packages:
gdal  
os  
rasterio  
xml.etree  
numpy  
datetime  
math  
subprocess  
shutil  
glob  

# sample usage:

This is a first version of the script for testing. In order to try it, the .py file should be place in the same folder of the SAOCOM .xemt files and unzipped data folders:

![1](1.png)

This script uses a class saocom_sm_slc and a function read_saocom that are used to read SAOCOM-1 Stripmap Data in SLC format to GMTSAR.  To test it, it can be run in interactive mode from the linux terminal.

    >>> python 3 -i saocom_sm_slc.py
    >> read_saocom()
    
This will automatically iterate through all the .xemt files and data folders and will create a directory for each polarization channel with the corresponding PRM,LED and SLC files. It is very important that the original information (with no modifications) is placed within the same folder in order to get this code to work properly.

![2](2.png)





