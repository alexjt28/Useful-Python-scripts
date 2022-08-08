#!/usr/bin/python

#----------------------------------------------------------------------------------------------------------------------
# This script reads in rainfall from a time series netCDF file of climate model output and determines the northernmost
# latitudes reached each month/year by the West African monsoon.
#
# To define the northern limit of the West African Monsoon (WAM), we use a threshold of 2 mm/day,
# adapted from 60 mm/month, averaged over the zonal region 15°W–20°E (Pausata et al., 2016).
# With monthly model output, we calculate the northernmost latitude that exceeds
# 2 mm/day for each month over the zonal region. We find the northernmost latitude for each year,
# taken as the maximum of the year’s monthly values, and create a distribution of northernmost latitudes
# by year for each simulation. The median of this distribution is calculated and printed at the end.
#
# Author: Alex Thompson
# Date: 8/8/2022
#
# Example shown here is for 100 years of monthly files from a case named "b.e12.B1850C5.f19_g16.iPI.01"
#----------------------------------------------------------------------------------------------------------------------

#*******************************************
# Import packages
#*******************************************

import numpy as np
import xarray as xr
import Ngl

#*******************************************
# Read in and define data
#*******************************************

# Open netCDF file containing input data
fname = "[/filepath/to/case/]b.e12.B1850C5.f19_g16.iPI.01.cam.h0.timeseries_of_all_vars.0801-0900.nc"
data = xr.open_dataset(fname)
#print(data.variables) # prints ALL variables in file, lots of text output

# Dimensions from netCDF file
time = data.time
lev  = data.lev
lat  = data.lat
lon  = data.lon
lonWE = lon.where(lon < 180.,lon-360)
#print(time)
#print(lev)
#print(lat)
#print(lon)
#print(lonWE)

# Create year and month dimensions
MON    = np.array([0,1,2,3,4,5,6,7,8,9,10,11])    # array of each month, integer
numyrs = np.array(len(time)/12).astype(dtype=int) # number of years, integer
nummon = np.array(len(MON)).astype(dtype=int)     # number of months, integer

# Create precipitation variable
prect = (data.PRECC + data.PRECL)*86400000. # prect=precc+precl, units: m/s -> mm/day
prect.attrs["units"] = "mm/day"             # define new units as mm/day
#print(prect)

#*************************************************
# Define bounded region to perform calculation on 
#*************************************************

westlon  = 345.  # western longitude for domain 
eastlon  = 20.   # eastern longitude for domain
southlat = 10.   # southern latitude for domain
northlat = 36.   # northern latitude for domain

lonw = list(lon.values).index(lon.sel(lon=westlon, method='nearest').lon) # get index value 
lone = list(lon.values).index(lon.sel(lon=eastlon, method='nearest').lon)
lats = list(lat.values).index(lat.sel(lat=southlat, method='nearest').lat)
latn = list(lat.values).index(lat.sel(lat=northlat, method='nearest').lat)
print("index = %d and coord = %d° (- = W, + = E)" %(lonw,lonWE.isel(lon=lonw)))
print("index = %d and coord = %d° (- = W, + = E)" %(lone,lonWE.isel(lon=lone)))
print("index = %d and coord = %d° (- = S, + = N)" %(lats,lat.isel(lat=lats)))
print("index = %d and coord = %d° (- = S, + = N)" %(latn,lat.isel(lat=latn)))

#*****************************************
# Prepare for performing zonal average 
#*****************************************

# If lon traverses prime meridian, need to create separate array that wraps end of array
if lonw > lone:
 lonarray_west = np.arange(lonw,list(lon.values).index(lon.sel(lon=lon[-1]).lon)+1,1) # from index lonw to end
 lonarray_east = np.arange(0,list(lon.values).index(lon.sel(lon=lon[lone]).lon)+1,1)  # from index 0 to lone
 lonarray = np.append(lonarray_west,lonarray_east)
 #print(lonarray)

# Index specific lat/lons, if/else specifies which lon indices to use
if lonw < lone:
 apply_coord_indices = prect[:,lats:latn,lonw:lone]
else:
 apply_coord_indices = prect[:,lats:latn,lonarray]

# Perform zonal average [time,lats:latn]
monsoonvar = apply_coord_indices.mean(axis=2) # averaging along longitudes
#print(monsoonvar)

#**************************
# Initialize variables
#**************************

# Create variables for # of years and latitude cells
times     = monsoonvar[:,0].size
numlats   = monsoonvar[0,:].size
latcoords = lat[lats-1:latn]
#print(latcoords)

# Initialize variables to run with loop, start northward_lat at 0 so it can climb
rain_per_lat = np.zeros(numlats)
northward_lat = np.zeros(times,dtype=int)

#******************************************************************************
# Run a loop to determine highest latitude reaching a specified rainfall limit
#******************************************************************************

# Set maximum rainfall value
max_rain = 2.0 # mm/day

# Run loop
for n in range(0,times):
 for i in range(0,numlats):
  rain_per_lat[:] = monsoonvar[n,:]
  if rain_per_lat[i] >= max_rain:
   northward_lat[n] = northward_lat[n]+1   # counts how many latitude bands (s->nn) exceed max_rain

# northward_lat now contains an index value for every month of the farthest north the max_rain value reached

#*************************************************************************
# Find max latitude (as index) that reached max_rain for each time period 
#*************************************************************************

# Initialize variables
timeperiod             = MON.size
numgroups              = times // timeperiod              # // means that 'numgroups' will be an integer
northlat_per_timegroup = np.zeros(numgroups,dtype=int)

# Run loop (p=year, a=beg_index_of_year, b=end_index_of_year)
for p in range(1,numgroups+1,1):
    b = (p*timeperiod)-1
    a = b - (timeperiod-1)
    northlat_per_timegroup[p-1] = max(northward_lat[a:b])

#**********************************************
# Convert indices into latitude coordinates
#**********************************************

# Initialize variable
northlat_vals = np.zeros(numgroups)

# Run loop
for c in range(0,len(northlat_per_timegroup)):
 northlat_vals[c] = latcoords[northlat_per_timegroup[c]]

# Print northernmost latitude value for each timegroup (aka year)
print("Northernmost latitudes per timegroup:")
print(northlat_vals)

#**************************
# Find median of values
#**************************

MEDIAN_MAXLAT = np.median(northlat_vals)
print("Median value = %d" %(MEDIAN_MAXLAT))
