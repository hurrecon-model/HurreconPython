# Copyright (C) President and Fellows of Harvard College 2020

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public
#   License along with this program.  If not, see
#   <http://www.gnu.org/licenses/>.

# The HURRECON Model estimates wind speed, wind direction, enhanced Fujita 
# scale wind damage, and duration of gale and hurricane winds as a function
# of hurricane location and maximum wind speed.

# Emery R. Boose
# March 2022

# Python version 3.7.11

# Note: Pandas datetime functions are currently limited to the years 
# 1678-2262 and so are not used here.

# Note: If setting PROJ_LIB as an environmental variable at the OS level
# causes problems with other programs, try setting it as below substituting
# the correct path for your system.


### MODULES ###############################################

import os

# set PROJ_LIB as Python environmental variable
os.environ['PROJ_LIB'] = "C:/Anaconda3/envs/hf/Library/share/proj"

import sys
import math
import time
import fiona
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors
import rasterio as rio
from rasterio.plot import show
from rasterio.mask import mask
from descartes import PolygonPatch


### INTERNAL FUNCTIONS ####################################

# get_fujita_wind_speeds returns a list containing the minimum 3-second 
# gust speed (meters/second) for each enhanced Fujita class.
# returns a list of 3-second gust speeds

def get_fujita_wind_speeds():
	ef0 = 29.1 # meters/second
	ef1 = 38.0 # meters/second
	ef2 = 49.2 # meters/second
	ef3 = 60.4 # meters/second
	ef4 = 73.8 # meters/second
	ef5 = 89.4 # meters/second

	return [ef0, ef1, ef2, ef3, ef4, ef5]

# get_fujita_colors returns a vector containing the color for each 
# enhanced Fujita class.
# returns a vector of colors

def get_fujita_colors():
	ef0_col = "purple"
	ef1_col = "blue"
	ef2_col = "green"
	ef3_col = "yellow"
	ef4_col = "orange"
	ef5_col = "red"
	efx_col = "grey"

	return [ef0_col, ef1_col, ef2_col, ef3_col, ef4_col, ef5_col, efx_col]

# format_time_difference_hms returns a time difference formatted as
# hours:minutes:seconds.
#   start_time - start time
#   end_time - end time
# returns a time difference formatted as hh:mm:ss.

def format_time_difference_hms(start_time, end_time):
	dsec = end_time - start_time

	hours = int(dsec//3600)
	minutes = int((dsec - 3600*hours)//60)
	seconds = round(dsec - 3600*hours - 60*minutes, 1)

	hh = str(hours)
	mm = str(minutes)
	ss = str(seconds)

	t_diff = hh.zfill(2) + ":" + mm.zfill(2) + ":" + ss.zfill(4) 

	return t_diff

# format_time_difference_ms returns a time difference in milliseconds.
#   start_time - start time
#   end_time - end time
# returns a time difference in milliseconds

def format_time_difference_ms(start_time, end_time):
	t_diff = str(round(1000*(end_time - start_time))).zfill(3)

	return t_diff

# check_file_exists displays an error message and stops execution if
# the specified file does not exist.
#   file_name - name of file
# no return value

def check_file_exists(file_name):
	if os.path.exists(file_name) == False or os.path.isfile(file_name) == False:
		sys.exit("File not found: " + file_name)

# check_legend_location checks if the specified legend location is valid.
#   loc - legend location
# returns True or False

def check_legend_location(loc):
    locations = ["upper left", "upper right", "lower left", "lower right"]

    return loc in locations

# read_site_file reads a site file and returns a list containing the
# latitude (degrees), longitude (degrees), and cover type (water=1, land=2) 
# for the specified site.
#   site_name - name of site
# returns a list of latitude, longitude, and cover type

def read_site_file(site_name):
	cwd = os.getcwd()
	site_file = cwd + "/input/sites.csv"
	check_file_exists(site_file)
	ss = pd.read_csv(site_file)

	# get site location & cover type
	index = np.where(ss.site_name == site_name)[0].tolist()
   
	if len(index) == 0:
		sys.exit("Site not found")

	site_latitude = ss.latitude[index[0]]
	site_longitude = ss.longitude[index[0]]
	cover_type = ss.cover_type[index[0]]

	return [site_latitude, site_longitude, cover_type]

# read_parameter_file reads a parameter file and returns a list containing
# the radius of maximum wind (rmw) (kilometers) and scaling parameter (profile 
# exponent) (s_par). If width is True, parameters are returned for the specified 
# hurricane, if available; otherwise parameters for ALL are returned.
#   hur_id - hurricane id
#   width - whether to use width parameters for the specified hurricane
# returns a list of rmw and s_par

def read_parameter_file(hur_id, width):
	cwd = os.getcwd()
	par_file = cwd + "/input/parameters.csv"
	check_file_exists(par_file)
	pp = pd.read_csv(par_file)

	# get rmw & s_par parameters
	if width:
		index = np.where(pp.hur_id == hur_id)[0].tolist()
		if len(index) == 0:
			sys.exit("Parameters missing for: " + hur_id)
	
	else:
		index = np.where(pp.hur_id == "ALL")[0].tolist()
		if len(index) == 0:
			sys.exit("Parameter file must contain an entry for ALL")

	rmw = pp.rmw[index[0]]
	s_par = pp.s_par[index[0]]

	return [rmw, s_par]

# get_fixed_model_parameters returns a list of fixed model parameters,
# including asymmetry factor, inflow angle, friction factor, and gust factor.
#   cover_type - cover type (1=water, 2=land)
# returns a list of fixed model parameters

def get_fixed_model_parameters(cover_type):
	asymmetry_factor = 1.0

	# water
	if cover_type == 1:
		inflow_angle = 20 # degrees
		friction_factor = 1.0
		gust_factor = 1.2

	#land
	elif cover_type == 2:
		inflow_angle = 40 # degrees
		friction_factor = 0.8
		gust_factor = 1.5

	else:
		sys.exit("Cover type must be 1 (water) or 2 (land)")

	return [asymmetry_factor, inflow_angle, friction_factor, gust_factor]

# get_time_step calculates the time step (minutes) for regional modeling, 
# assuming a maximum hurricane forward speed of 20 meters per second (1200 
# meters per minute). Values are rounded to the nearest 1, 2, 3, 5, 10, 15, 
# 30, or 60 minutes.
# returns a time step in minutes

def get_time_step():
	# read land-water file
	cwd = os.getcwd()
	land_water_file = cwd + "/input/land_water.tif"
	check_file_exists(land_water_file)
	land_water = rio.open(land_water_file)

	# get cell height in meters (at latitude 45 degrees)
	nrows = land_water.height
	lat_min = land_water.bounds.bottom
	lat_max = land_water.bounds.top
	cell_y = 111132 * (lat_max - lat_min)/nrows

	# close land-water file
	land_water.close()

	# calculate time step
	ts = round(cell_y/1200)

	if ts <= 1:
		time_step = 1
	elif ts <= 2:
		time_step = 2
	elif ts <= 4:
		time_step = 3
	elif ts <= 7:
		time_step = 5
	elif ts <= 12:
		time_step = 10
	elif ts <= 22:
		time_step = 15
	elif ts <= 45:
		time_step = 30
	else:
		time_step = 60

	return time_step

# calculate_julian_day calculates the Julian day and fraction from
# the specified year, month, day, hour, and minute.
#   year - year
#   month - month
#   day - day
#   hour - hour
#   minute - minute
# returns the Julian day and fraction

def calculate_julian_day(year, month, day, hour, minute):
	# get Julian day
	if month == 1:
		jd = 0
	elif month == 2:
		jd = 31
	elif month == 3:
		jd = 59
	elif month == 4:
		jd = 90
	elif month == 5:
		jd = 120
	elif month == 6:
		jd = 151
	elif month == 7:
		jd = 181
	elif month == 8:
		jd = 212
	elif month == 9:
		jd = 243
	elif month == 10:
		jd = 273
	elif month == 11:
		jd = 304
	elif month == 12:
		jd = 334

	jd = jd + day

	# correct for leap year
	if year % 4 == 0 and ((year % 100 != 0) or (year % 400 == 0)):
		if month >= 3:
			jd = jd + 1

	# get fraction of day
	jd = jd + hour/24 + minute/1440

	return jd

# calculate_standard_datetime returns a standard datetime in the format
# YYYY-MM-DDThh:mm from the specified year and Julian day with fraction.
#   year - year
#   jd - Julian day with fraction
# returns a standard datetime

def calculate_standard_datetime(year, jd):
	# separate day & fraction
	jd_int = int(math.floor(jd))
	jd_frac = jd - math.floor(jd)

	# correct for leap year
	if year % 4 == 0 and ((year % 100 != 0) or (year % 400 == 0)):
		leap_year = True
	else:
		leap_year = False

	if (jd_int > 60 and leap_year == True):
		jd_int = jd_int - 1

	# get month
	if jd_int <= 31:
		mon = 1
		day = jd_int
	elif jd_int <= 59:
		mon = 2
		day = jd_int - 31
	elif jd_int <= 90:
		mon = 3
		day = jd_int - 59
	elif jd_int <= 120:
		mon = 4
		day = jd_int - 90
	elif jd_int <= 151:
		mon = 5
		day = jd_int - 120
	elif jd_int <= 181:
		mon = 6
		day = jd_int - 151
	elif jd_int <= 212:
		mon = 7
		day = jd_int - 181
	elif jd_int <= 243:
		mon = 8
		day = jd_int - 212
	elif jd_int <= 273:
		mon = 9
		day = jd_int - 243
	elif jd_int <= 304:
		mon = 10
		day = jd_int - 273
	elif jd_int <= 334:
		mon = 11
		day = jd_int - 304
	else:
		mon = 12
		day = jd_int - 334
    
	hour = int(math.floor(jd_frac*24))
	minute = int(round(jd_frac*1440 - hour*60))

	if minute == 60:
		minute = 0
		hour = hour + 1

	if hour == 24:
		hour = 0
		day = day + 1

	if (day == 32) or (day == 31 and mon in (4, 6, 9, 11)):
		day = 1
		mon = mon + 1

	if (mon == 2) and ((day == 29 and leap_year == False) or (day == 30)):
		day = 1
		mon = mon + 1

	if mon == 13:
		mon = 1
		year = year + 1

	datetime = str(year) + "-" + str(mon).zfill(2) + "-" + str(day).zfill(2) +  "T" + str(hour).zfill(2) + ":" + str(minute).zfill(2)

	return datetime

# read_hurricane_track_file reads a hurricane track file and returns
# a data frame of track data for the specfied hurricane.
#   hur_id - hurricane id
# returns a data frame of track data

def read_hurricane_track_file(hur_id):
	# read hurricane track file
	cwd = os.getcwd()
	track_file = cwd + "/input/tracks.csv"
	check_file_exists(track_file)
	zz = pd.read_csv(track_file)

	# subset by hurricane name
	xx = zz.loc[zz.hur_id == hur_id]

	if len(xx) == 0:
		sys.exit("Hurricane not in track file")

	# create data frame
	hur_id_list = list(xx.hur_id)
	name_list = list(xx.name)
	date_time_list = list(xx.date_time)
	jd_list = list(xx.jd)
	status_list = list(xx.status)
	latitude_list = list (xx.latitude)
	longitude_list = list(xx.longitude)
	wind_max_list = list(xx.wind_max)

	tt_cols = ['hur_id', 'name', 'date_time', 'jd', 'status', 'latitude', 'longitude', 'wind_max']

	tt = pd.DataFrame(data=list(zip(hur_id_list, name_list, date_time_list, jd_list, status_list,
		latitude_list, longitude_list, wind_max_list)), columns=tt_cols)

	return tt 

# interpolate_hurricane_location_max_wind performs a linear interpolation
# of hurricane latitude & longitude (degrees) and maximum sustained wind
# speed (meters/second) using data from a hurricane track file and a 
# specified time step.
#   tt - data frame of track data
#   time_step - time step (minutes)
# returns a list containing lists of year, Julian day, latitude, longitude, 
# and max wind speed

def interpolate_hurricane_location_max_wind(tt, time_step):
	tt_rows = len(tt)

	# initialize lists
	yr_list = []
	jd_list = []
	lat_list = []
	lon_list = []
	wmax_list = []

	# interpolate values for each segment of track
	for i in range(0, tt_rows-1):
		new_rows = int(round(1440*(tt.jd[i+1] - tt.jd[i])/(time_step))) + 1

		jd = np.linspace(start=tt.jd[i], stop=tt.jd[i+1], num=new_rows).tolist()
		lat = np.linspace(start=tt.latitude[i], stop=tt.latitude[i+1], num=new_rows).tolist()
		lon = np.linspace(start=tt.longitude[i], stop=tt.longitude[i+1], num=new_rows).tolist()
		wmax = np.linspace(start=tt.wind_max[i], stop=tt.wind_max[i+1], num=new_rows).tolist()

		# remove last element to avoid duplication
		del jd[-1]
		del lat[-1]
		del lon[-1]
		del wmax[-1]

		yr = [int(tt.date_time[0][0:4])] * (new_rows - 1)

		yr_list.extend(yr)
		jd_list.extend(jd)
		lat_list.extend(lat)
		lon_list.extend(lon)
		wmax_list.extend(wmax)

	# add final row
	yr_list.append(int(tt.date_time[tt_rows-1][0:4]))
	jd_list.append(tt.jd[tt_rows-1])
	lat_list.append(tt.latitude[tt_rows-1])
	lon_list.append(tt.longitude[tt_rows-1])
	wmax_list.append(tt.wind_max[tt_rows-1])

	return [yr_list, jd_list, lat_list, lon_list, wmax_list]

# estimate_range uses the Pythagorean equation to estimate the range 
# (kilometers) from one point to another based on the latitude & longitude 
# of each point. Note: overestimates range unless on both points are on
# the same meridian.
#   lat1 - latitude of first point (degrees)
#   lon1 - longitude of first point (degrees)
#   lat2 - latitude of second point (degrees)
#   lon2 - longitude of second point (degrees)
# returns range in kilometers

def estimate_range(lat1, lon1, lat2, lon2):
	R = 6367 # radius of earth in kilometers (at latitude 45 degrees)
	d2r = 0.017453292519943295  # pi / 180

	lat_avg = d2r*(lat1 + lat2)/2
	x = d2r*(lon1 - lon2)*math.cos(d2r*lat_avg)
	y = d2r*(lat1 - lat2)
	range_est = R * math.sqrt(x**2 + y**2)

	return range_est

# calculate_range uses the Haversine formula to calculate the range 
# (kilometers) from one point to another based on the latitude & longitude
# of each point.
#   lat1 - latitude of first point (degrees)
#   lon1 - longitude of first point (degrees)
#   lat2 - latitude of second point (degrees)
#   lon2 - longitude of second point (degrees)
# returns range in kilometers

def calculate_range(lat1, lon1, lat2, lon2):
	R = 6367 # radius of earth in kilometers (at latitude 45 degrees)
	d2r = 0.017453292519943295  # pi / 180

	# nearly same point
	if abs(lat2 - lat1) < 0.000001 and abs(lon2 - lon1) < 0.000001:
		rang = 0

	else:
		# date line
		if lon1 > 90 and lon2 < -90:
			lon2 <- lon2 + 360
		if lon1 < -90 and lon2 > 90:
			lon1 <- lon1 + 360

		# convert degrees to radians
		rlat1 = d2r*lat1
		rlat2 = d2r*lat2
		rlon1 = d2r*lon1
		rlon2 = d2r*lon2

		A = (math.sin((rlat2-rlat1)/2))**2 + math.cos(rlat1)*math.cos(rlat2)*(math.sin((rlon2-rlon1)/2))**2
		C = 2 * math.atan2(math.sqrt(A), math.sqrt(1-A))
		rang = R * C

	return rang

# calculate_bearing uses the Haversine formula to calculate the bearing 
# (degrees) from one point to another based on the latitude & longitude 
# of each point.
#   lat1 - latitude of first point (degrees)
#   lon1 - longitude of first point (degrees)
#   lat2 - latitude of second point (degrees)
#   lon2 - longitude of second point (degrees)
# returns bearing in degrees

def calculate_bearing(lat1, lon1, lat2, lon2):
	d2r = 0.017453292519943295  # pi / 180
	r2d = 57.29577951308232  # 180 / pi

	# nearly same point
	if abs(lat2 - lat1) < 0.000001 and abs(lon2 - lon1) < 0.000001:
		bear = 0

	else:
		# date line
		if lon1 > 90 and lon2 < -90:
			lon2 <- lon2 + 360
		if lon1 < -90 and lon2 > 90:
			lon1 <- lon1 + 360

		# same longitude
		if lon1 == lon2:
			if lat1 > lat2:
				bear = 180
			else:
				bear = 0

		# different longitude
		else:
			# convert degrees to radians
			rlat1 = d2r*lat1
			rlat2 = d2r*lat2
			rlon1 = d2r*lon1
			rlon2 = d2r*lon2

			B2 = math.atan2(math.sin(rlon2-rlon1)*math.cos(rlat2), math.cos(rlat1)*math.sin(rlat2) - 
				math.sin(rlat1)*math.cos(rlat2)*math.cos(rlon2-rlon1))

			# convert radians to degrees
			B = r2d*B2

			if lon1 < lon2:
				# quadrants I, IV
				bear = B
			else:  
				# quadrants II, III
				bear = 360 + B

	return bear

# get_maximum_wind_speed returns the maximum sustained wind speed for
# the specified hurricane.
#   hur_id - hurricane id
# returns maximum sustained wind speed (meters/second)

def get_maximum_wind_speed(hur_id):
	# read hurricane track file
	cwd = os.getcwd()
	track_file = cwd + "/input/tracks.csv"
	check_file_exists(track_file)
	zz = pd.read_csv(track_file)

	# subset by hurricane name
	xx = zz.loc[zz.hur_id == hur_id]

	# get maximum wind speed
	wmax = max(xx.wind_max)

	return wmax

# get_maximum_range estimates the range (kilometers) at which sustained wind
# speeds are less than gale (17.5 meters/second).
#   wmax - maximum sustained wind speed (meters/second)
#   rmw - radius of maximum winds (kilometers)
#   s_par - scaling parameter
# returns range in kilometers

def get_maximum_range(wmax, rmw, s_par):
	rang = rmw
	wspd = 100

	while wspd >= 17.5:
		rang = rang + 10
		x = (rmw/rang)**s_par
		wspd = wmax * math.sqrt(x * math.exp(1-x))

	return rang

# interpolate_hurricane_speed_bearing performs a linear interpolation of hurricane
# speed (meters/second) and bearing (degrees) along a hurricane track based on
# mid-segment values.
#   tt - data frame of track values
#   jd_list - list of Julian day values
# returns a list containing lists of hurricane speed & bearing

def interpolate_hurricane_speed_bearing(tt, jd_list):
	tt_rows = len(tt)

	# intialize lists
	vv_jd = []
	vv_spd = []
	vv_bear = []

	# calculate mid-segment hurricane speed & bearing
	for i in range(0, tt_rows-1):
		hur_range = calculate_range(tt.latitude[i], tt.longitude[i], 
			tt.latitude[i+1], tt.longitude[i+1])

		hur_bear = calculate_bearing(tt.latitude[i], tt.longitude[i], 
			tt.latitude[i+1], tt.longitude[i+1])

		interval_sec = (tt.jd[i+1] - tt.jd[i]) * 1440 * 60

		vv_jd.append(tt.jd[i] + (tt.jd[i+1] - tt.jd[i])/2)
		vv_spd.append(1000*hur_range/interval_sec)
		vv_bear.append(hur_bear)

	# get number of rows
	vv_rows = len(vv_jd)

	# intialize lists
	bear_list = []
	spd_list = []

 	# interpolate hurricane speed & bearing for each segment
	for i in range(0, vv_rows+1):
		# before mid-point of 1st segment
		if i == 0:
			index = np.where(jd_list <= vv_jd[0])[0]
			new_rows = len(index)

			bear = [vv_bear[0]] * new_rows
			spd = [vv_spd[0]] * new_rows

			bear_list.extend(bear)
			spd_list.extend(spd)

		# interpolate between mid-points
		elif i <= vv_rows-1:
			index = np.where((jd_list > vv_jd[i-1]) & (jd_list <= vv_jd[i]))[0]
			new_rows = len(index)

			# bearing
			b1 = vv_bear[i-1]
			b2 = vv_bear[i]

			if b2 - b1 > 180:
				b1 = b1 + 360
			elif b1 - b2 > 180:
				b2 = b2 + 360

			bear = np.linspace(start=b1, stop=b2, num=new_rows).tolist()
			
			# speed
			spd = np.linspace(start=vv_spd[i-1], stop=vv_spd[i], num=new_rows).tolist()

			bear_list.extend(bear)
			spd_list.extend(spd)

		# after mid-point of last segment
		else:
			index = np.where(jd_list > vv_jd[vv_rows-1])[0]
			new_rows = len(index)

			bear = [vv_bear[vv_rows-1]] * new_rows
			spd = [vv_spd[vv_rows-1]] * new_rows

			bear_list.extend(bear)
			spd_list.extend(spd)
			
	# adjust bearing as needed
	for i in range(0, len(bear_list)):
		if bear_list[i] < 0:
			bear_list[i] = bear_list[i] + 360

		if bear_list[i] > 360:
			bear_list[i] = bear_list[i] - 360

	return [spd_list, bear_list]

# calculate_site_range_bearing calculates the range (kilometers) and bearing
# (degrees) from a site to the hurricane center.
#   lat_list - list of hurricane latitudes (degrees)
#	lon_list - list of hurricane longitudes (degrees)
#   site_latitude - latitude of site (degrees)
#   site_longitude - longitude of site (degrees)
# returns a list containing lists of site range & bearing

def calculate_site_range_bearing(lat_list, lon_list, site_latitude, site_longitude):
	num = len(lat_list)

	site_range = [calculate_range(site_latitude, site_longitude, 
		lat_list[i], lon_list[i]) for i in range(0, num)]

	site_bear = [calculate_bearing(site_latitude, site_longitude, 
		lat_list[i], lon_list[i]) for i in range(0, num)]

	return [site_range, site_bear]

# calculate_wind_direction calculates the wind direction (degrees) at the
# specified site.
#   hurr_lat - latitude of hurricane (degrees)
#   site_bear - bearing from site to hurricane center (degrees)
#   inflow_angle - cross-isobar inflow angle (degrees)
# returns a calculated wind direction in degrees

def calculate_wind_direction (hurr_lat, site_bear, inflow_angle):
	# northern hemisphere: tangent minus inflow angle
	if hurr_lat > 0:
		wind_dir = site_bear - 90 - inflow_angle
		if wind_dir < 0:
			wind_dir = wind_dir + 360
  
	# southern hemisphere: tangent plus inflow angle
	else:
		wind_dir = site_bear + 90 + inflow_angle
		if wind_dir > 360:
			wind_dir = wind_dir - 360

	return wind_dir

# calculate_wind_speed calculates the sustained wind speed (meters/second) at
# the specified site.
#   site_bear - bearing from site to hurricane center (degrees)
#   site_range - range from site to hurricance center (kilometers)
#   hur_lat - latitude of hurricane (degrees)
#   hur_bear - hurricane bearing (degrees)
#   hur_spd - hurricane speed (meters/second)
#   wind_max - maximum sustained wind speed (meters/second)
#   rmw - radius of maximum winds (kilometers)
#   s_par - scaling parameter
#   asymmetry_factor - asymmetry factor
#   friction_factor - friction factor
# returns a calculated sustained wind speed (meters/second)

def calculate_wind_speed(site_bear, site_range, hur_lat, hur_bear, hur_spd, 
 	wind_max, rmw, s_par, asymmetry_factor, friction_factor):
	
	# hurricane eye (avoid division by zero)
	if site_range == 0:
		wind_spd = 0 

	else:
		# northern hemisphere: clockwise angle from path
		if hur_lat > 0:
			T = site_bear - hur_bear + 180

		# southern hemisphere: counterclockwise angle from path
		else:
			T = site_bear - hur_bear

		X = (rmw / site_range)**s_par

		# sustained wind speed at radius of maximum wind (rmw)
		Z = wind_max - hur_spd * asymmetry_factor * (1 - math.sin(T * math.pi/180))/2

		if Z < 0:
			Z = 0
      
		# sustained wind speed at site
		wind_spd = Z * math.sqrt(X * math.exp(1 - X))

		# adjust for land or water
		wind_spd = wind_spd * friction_factor

	return wind_spd

# calculate_wind_gust calculates the wind gust speed (meters/second) from 
# the sustained wind speed (meters/second) and the gust factor.
#   wind_spd - sustained wind speed (meters/second)
#   gust_factor - gust factor
# returns wind gust speed (meters/second)

def calculate_wind_gust(wind_spd, gust_factor):
	gust_spd = gust_factor * wind_spd

	return gust_spd

# calculate_enhanced_fujita_scale returns the enhanced Fujita scale value
# based on the wind gust speed (meters/second).
#   gust_spd - wind gust speed (meters/second)
# returns the enhanced Fujita scale value

def calculate_enhanced_fujita_scale(gust_spd):
	# get enhanced Fujita wind speeds
	ef = get_fujita_wind_speeds()

	if gust_spd < ef[0]:
		ef_sca = -1
	elif gust_spd < ef[1]:
		ef_sca = 0
	elif gust_spd < ef[2]:
		ef_sca = 1
	elif gust_spd < ef[3]:
		ef_sca = 2
	elif gust_spd < ef[4]:
		ef_sca = 3
	elif gust_spd < ef[5]:
		ef_sca = 4
	else:
		ef_sca = 5

	return ef_sca

# calculate_wind_speed_direction calculates the wind speed, gust speed, wind
# direction, and enhanced Fujita scale wind damage at a site.
#   sbear_list - list of site bearings (degrees)
#	srange_list - list of site ranges (kilometers)
#	lat_list - list of hurricane latitudes (degrees)
#	bear_list - list of hurricane bearings (degrees)
#	spd_list - list of hurricane forward speeds (meters/second)
#	wmax_list - list of maximum sustained wind speeds (meters/second)
#   inflow_angle - cross-isobar inflow angle (degrees)
#   rmw - radius of maximum winds (kilometers)
#   s_par - scaling parameter
#   asymmetry_factor - asymmetry factor
#   friction_factor - friction factor
#   gust_factor - gust factor
# returns a list containing lists of wind speed, gust speed, wind direction, 
# and enhanced Fujita scale

def calculate_wind_speed_direction(sbear_list, srange_list, lat_list, bear_list,
	spd_list, wmax_list, inflow_angle, rmw, s_par, asymmetry_factor, friction_factor, 
	gust_factor):

	num = len(sbear_list)

	# wind speed
	wspd_list = [calculate_wind_speed(sbear_list[i], srange_list[i], lat_list[i], 
		bear_list[i], spd_list[i], wmax_list[i], rmw, s_par, asymmetry_factor, 
		friction_factor) for i in range(0, num)]

	# gust speed
	gspd_list = [calculate_wind_gust(wspd_list[i], gust_factor) for i in range(0, num)]

	# wind direction
	wdir_list = [calculate_wind_direction(lat_list[i], sbear_list[i], inflow_angle)
		for i in range(0, num)]

	# enhanced Fujita scale
	ef_list = [calculate_enhanced_fujita_scale(gspd_list[i]) for i in range(0, num)]

	return [wspd_list, gspd_list, wdir_list, ef_list]

# get_standard_date_time creates a list of standard datetime values in the format
# YYYY-MM-DDThh:mm.
#   yr_list - list of years
#	jd_list - list of Julian days
# returns a list of datetime values

def get_standard_date_time(yr_list, jd_list):
	num = len(yr_list)

	dt_list = [calculate_standard_datetime(yr_list[i], jd_list[i])
		for i in range(0, num)]
	
	return dt_list

# get_peak_values returns a data frame of peak values for a given
# hurricane and site.
#   hur_id - hurricane id
#   site_name - name of site
#   mm - data frame of modeled values
# returns a data frame of peak values

def get_peak_values(hur_id, site_name, mm):
	
	# get time step in minutes
	h1 = int(mm.date_time[0][12:13])
	m1 = int(mm.date_time[0][15:16])
	t1 = h1 * 60 + m1

	h2 = int(mm.date_time[1][12:13])
	m2 = int(mm.date_time[1][15:16])
	t2 = h2 * 60 + m2

	time_step = t2 - t1

	# get peak wind
	pk = mm.loc[mm.wind_spd == max(mm.wind_spd), ]

	date_time = pk.date_time.values[0]
	wind_dir = pk.wind_dir.values[0]
	wind_spd = pk.wind_spd.values[0]
	gust_spd = pk.gust_spd.values[0]
	ef_sca = pk.ef_sca.values[0]

	# get wind duration in hours
	ef0_obs = mm.loc[mm.ef_sca >= 0, ]
	ef1_obs = mm.loc[mm.ef_sca >= 1, ]
	ef2_obs = mm.loc[mm.ef_sca >= 2, ]
	ef3_obs = mm.loc[mm.ef_sca >= 3, ]
	ef4_obs = mm.loc[mm.ef_sca >= 4, ]
	ef5_obs = mm.loc[mm.ef_sca >= 5, ]

	ef0 = len(ef0_obs) * time_step/60
	ef1 = len(ef1_obs) * time_step/60
	ef2 = len(ef2_obs) * time_step/60
	ef3 = len(ef3_obs) * time_step/60
	ef4 = len(ef4_obs) * time_step/60
	ef5 = len(ef5_obs) * time_step/60

	# create data fame of peak values
	kk_cols = ['site_name', 'hur_id', 'date_time', 'wind_dir', 'wind_spd',
		'gust_spd', 'ef_sca', 'ef0', 'ef1', 'ef2', 'ef3', 'ef4', 'ef5']

	kk = pd.DataFrame(columns=kk_cols)

	kk.loc[0] = [site_name, hur_id, date_time, wind_dir, wind_spd, gust_spd, ef_sca, 
		ef0, ef1, ef2, ef3, ef4, ef5]

	return kk

# get_regional_peak_wind calculates peak values for wind speed (meters/second), 
# enhanced Fujita scale, wind direction (degrees), cardinal wind direction (1-8),
# and duration of EF0, EF1, EF2, EF3, EF4, and EF5 winds (minutes) for a given 
# hurricane over a region. Results are returned in a list.
#   hur_id - hurricane id
#	lat_list - list of hurricane latitudes (degrees)
#	lon_list - list of hurricane longitudes (degrees)
#	wmax_list - list of maximum sustained wind speeds (meters/second)
#	bear_list - list of hurricane bearings (degrees)
#	spd_list - list of hurricane forward speeds (meters/second)
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   water - whether to calculate values over water
#   console - whether to display messages in console
# returns a list containing 10 raster arrays

def get_regional_peak_wind(hur_id, lat_list, lon_list, wmax_list, bear_list, 
	spd_list, width, time_step, water, console):

	# get number of positions
	num = len(lat_list)
	
	# read land-water file
	cwd = os.getcwd()
	land_water_file = cwd + "/input/land_water.tif"
	check_file_exists(land_water_file)
	land_water = rio.open(land_water_file)

	# get regional values
	nrows = land_water.height
	ncols = land_water.width

	lat_min = land_water.bounds.bottom
	lat_max = land_water.bounds.top
	cell_y = (lat_max - lat_min)/nrows 

	lon_min = land_water.bounds.left
	lon_max = land_water.bounds.right
	cell_x = (lon_max - lon_min)/ncols
 
	# create arrays for peak wind speed & direction
	ss = np.zeros((nrows, ncols), dtype=np.int16)   # wind speed (m/s)
	ff = np.zeros((nrows, ncols), dtype=np.int16)   # enhanced Fujita scale
	dd = np.zeros((nrows, ncols), dtype=np.int16)   # wind direction (degrees)
	cc = np.zeros((nrows, ncols), dtype=np.int16)   # cardinal wind direction (1-8)
	f0 = np.zeros((nrows, ncols), dtype=np.int16)   # duration of EF0 winds (minutes)
	f1 = np.zeros((nrows, ncols), dtype=np.int16)   # duration of EF1 winds (minutes)
	f2 = np.zeros((nrows, ncols), dtype=np.int16)   # duration of EF2 winds (minutes)
	f3 = np.zeros((nrows, ncols), dtype=np.int16)   # duration of EF3 winds (minutes)
	f4 = np.zeros((nrows, ncols), dtype=np.int16)   # duration of EF4 winds (minutes)
	f5 = np.zeros((nrows, ncols), dtype=np.int16)   # duration of EF5 winds (minutes)

	xx = np.zeros((nrows, ncols), dtype=np.float64) # floating point wind speed (m/s)

	# create array from raster
	land_water_array = land_water.read(1)
	land_water.close()

	# read parameters file
	pars = read_parameter_file(hur_id, width)
	rmw = pars[0]
	s_par = pars[1]

	# get fixed model parameters by cover type (water=1, land=2)
	asymmetry = [get_fixed_model_parameters(1)[0], get_fixed_model_parameters(2)[0]]
	inflow = [get_fixed_model_parameters(1)[1], get_fixed_model_parameters(2)[1]]
	friction = [get_fixed_model_parameters(1)[2], get_fixed_model_parameters(2)[2]]
	gust = [get_fixed_model_parameters(1)[3], get_fixed_model_parameters(2)[3]]

	# get maximum wind speed over track
	wmax_track = get_maximum_wind_speed(hur_id)

	# get maximum range for gale winds
	range_maximum = get_maximum_range(wmax_track, rmw, s_par)

	# record total elasped time
	start_time = time.time()

	# calculate peak wind speed & direction and gale & hurricane duration for each location
	for i in range(0, nrows):
		for j in range(0, ncols):
			# get cover type from land_water layer
			cover_type = int(land_water_array[nrows-i-1][j])

			if cover_type == 2 or water == True:
				# get site latitude & longitude
				site_latitude = lat_min + (i + 0.5)*cell_y
				site_longitude = lon_min + (j + 0.5)*cell_x

				asymmetry_factor = asymmetry[cover_type-1]
				inflow_angle = inflow[cover_type-1]
				friction_factor = friction[cover_type-1]
				gust_factor = gust[cover_type-1]

				for k in range(0, num):
					hur_latitude  = lat_list[k]
					hur_longitude = lon_list[k]

					# site range
					site_range = calculate_range(site_latitude, site_longitude,
						hur_latitude, hur_longitude)

					# skip if too far away
					if site_range < range_maximum:
						# site bearing
						site_bear = calculate_bearing(site_latitude, site_longitude,
							hur_latitude, hur_longitude)

						# wind speed (m/s)
						wspd = calculate_wind_speed(site_bear, site_range, hur_latitude, 
							bear_list[k], spd_list[k], wmax_list[k], rmw, s_par, 
							asymmetry_factor, friction_factor)

						# update values if gale or higher
						if wspd >= 17.5:
							# enhanced Fujita scale
							gspd = calculate_wind_gust(wspd, gust_factor)
							fsca = calculate_enhanced_fujita_scale(gspd)

                            # update duration (minutes)
							if fsca >= 0:
								f0[nrows-i-1][j] = f0[nrows-i-1][j] + time_step
							if fsca >= 1:
								f1[nrows-i-1][j] = f1[nrows-i-1][j] + time_step
							if fsca >= 2:
								f2[nrows-i-1][j] = f2[nrows-i-1][j] + time_step
							if fsca >= 3:
								f3[nrows-i-1][j] = f3[nrows-i-1][j] + time_step
							if fsca >= 4:
								f4[nrows-i-1][j] = f4[nrows-i-1][j] + time_step
							if fsca >= 5:
								f5[nrows-i-1][j] = f5[nrows-i-1][j] + time_step

							# peak wind speed
							if xx[nrows-i-1][j] < wspd:							
								xx[nrows-i-1][j] = wspd

								ss[nrows-i-1][j] = int(round(wspd))

								# peak wind direction			
								wdir = calculate_wind_direction(hur_latitude, site_bear, inflow_angle)
								dd[nrows-i-1][j] = int(round(wdir))

		# report progress
		if console == True:
			x = round(i*100/nrows)
			if x % 10 == 0:
				print("          ", end="")
				print("\r", x, "%", end="")

	# calculate other values
	for i in range(0, nrows):
		for j in range(0, ncols):
			# get cover type from land_water layer
			cover_type = int(land_water_array[nrows-i-1][j])

			if cover_type == 2 or water == True:
				wspd = xx[nrows-i-1][j]

				# update values if gale or higher
				if wspd >= 17.5:
					asymmetry_factor = asymmetry[cover_type-1]
					inflow_angle = inflow[cover_type-1]
					friction_factor = friction[cover_type-1]
					gust_factor = gust[cover_type-1]

					# enhanced Fujita scale
					gspd = calculate_wind_gust(wspd, gust_factor)
					fsca = calculate_enhanced_fujita_scale(gspd)
					ff[nrows-i-1][j] = fsca + 2

					# cardinal wind direction (1 = north, 2 = northeast, etc)
					wdir = dd[nrows-i-1][j]
					cdir = math.floor((wdir+22.5)/45) + 1
					if cdir > 8:
						cdir = 1
					cc[nrows-i-1][j] = cdir

	# report elapsed time
	if console == True:
		elapsed_time = format_time_difference_hms(start_time, time.time())
		print("\r", elapsed_time, "\n", end="")

	return [ss, ff, dd, cc, f0, f1, f2, f3, f4, f5]

# get_regional_datetime calculates wind speed (meters/second), enhanced 
# Fujita scale, wind direction (degrees), and cardinal wind direction for 
# a given hurricane over a region at a specified datetime. Results are 
# returned as a list of 4 raster arrays.
#   hur_id - hurricane id
#   lat - hurricane latitude (degrees)
#   lon - hurricane longitude (degrees)
#   wmax - maximum sustained wind speed (meters/second)
#   bear - hurricane bearing (degrees)
#   spd - hurricane forward speed (meters/second)
#   width - whether to use width parameters for the specified hurricane
#   water - whether to calculate values over water
#   console - whether to display progress in console
# returns a list containing 4 raster layers

def get_regional_datetime(hur_id, lat, lon, wmax, bear, spd, width, 
  water, console):

	# read land-water file
	cwd = os.getcwd()
	land_water_file = cwd + "/input/land_water.tif"
	check_file_exists(land_water_file)
	land_water = rio.open(land_water_file)

	# get regional values
	nrows = land_water.height
	ncols = land_water.width

	lat_min = land_water.bounds.bottom
	lat_max = land_water.bounds.top
	cell_y = (lat_max - lat_min)/nrows 

	lon_min = land_water.bounds.left
	lon_max = land_water.bounds.right
	cell_x = (lon_max - lon_min)/ncols
 
	# create arrays for peak wind speed & direction
	ss = np.zeros((nrows, ncols), dtype=np.int16)  # wind speed (m/s)
	ff = np.zeros((nrows, ncols), dtype=np.int16)  # enhanced Fujita scale
	dd = np.zeros((nrows, ncols), dtype=np.int16)  # wind direction (degrees)
	cc = np.zeros((nrows, ncols), dtype=np.int16)  # cardinal wind direction (1-8)

	# create array from raster
	land_water_array = land_water.read(1)
	land_water.close()

	# read parameters file
	pars = read_parameter_file(hur_id, width)
	rmw = pars[0]
	s_par = pars[1]

	# get fixed model parameters by cover type (water=1, land=2)
	asymmetry = [get_fixed_model_parameters(1)[0], get_fixed_model_parameters(2)[0]]
	inflow = [get_fixed_model_parameters(1)[1], get_fixed_model_parameters(2)[1]]
	friction = [get_fixed_model_parameters(1)[2], get_fixed_model_parameters(2)[2]]
	gust = [get_fixed_model_parameters(1)[3], get_fixed_model_parameters(2)[3]]

	# get maximum wind speed over track
	wmax_track = get_maximum_wind_speed(hur_id)

	# get maximum range for gale winds
	range_maximum = get_maximum_range(wmax_track, rmw, s_par)

	# record total elasped time
	start_time = time.time()

	# calculate peak wind speed & direction and gale & hurricane duration for each location
	for i in range(0, nrows):
		for j in range(0, ncols):
			# get cover type from land_water layer
			cover_type = int(land_water_array[nrows-i-1][j])

			if cover_type == 2 or water == True:
				# get site latitude & longitude
				site_latitude = lat_min + (i + 0.5)*cell_y
				site_longitude = lon_min + (j + 0.5)*cell_x

				asymmetry_factor = asymmetry[cover_type-1]
				inflow_angle = inflow[cover_type-1]
				friction_factor = friction[cover_type-1]
				gust_factor = gust[cover_type-1]

				hur_latitude  = lat
				hur_longitude = lon

				# site range
				site_range = calculate_range(site_latitude, site_longitude,
					hur_latitude, hur_longitude)

				# skip if too far away
				if site_range < range_maximum:
					# site bearing
					site_bear = calculate_bearing(site_latitude, site_longitude,
						hur_latitude, hur_longitude)

					# wind speed (m/s)
					wspd = calculate_wind_speed(site_bear, site_range, hur_latitude, bear, 
						spd, wmax, rmw, s_par, asymmetry_factor, friction_factor)

					# update values if gale or higher
					if wspd >= 17.5:
						# wind speed (m/s)
						ss[nrows-i-1][j] = int(round(wspd))

						# wind direction (degrees)			
						wdir = calculate_wind_direction(hur_latitude, site_bear, inflow_angle)
						dd[nrows-i-1][j] = int(round(wdir))

		# report progress
		if console == True:
			x = round(i*100/nrows)
			if x % 10 == 0:
				print("          ", end="")
				print("\r", x, "%", end="")

	# calculate other values
	for i in range(0, nrows):
		for j in range(0, ncols):
			# get cover type from land_water layer
			cover_type = int(land_water_array[nrows-i-1][j])

			if cover_type == 2 or water == True:
				wspd = ss[nrows-i-1][j]

				# update values if gale or higher
				if wspd >= 17.5:
					asymmetry_factor = asymmetry[cover_type-1]
					inflow_angle = inflow[cover_type-1]
					friction_factor = friction[cover_type-1]
					gust_factor = gust[cover_type-1]

					# enhanced Fujita scale
					gspd = calculate_wind_gust(wspd, gust_factor)
					fsca = calculate_enhanced_fujita_scale(gspd)
					ff[nrows-i-1][j] = fsca + 2

					# cardinal wind direction (1 = north, 2 = northeast, etc)
					wdir = dd[nrows-i-1][j]
					cdir = math.floor((wdir+22.5)/45) + 1
					if cdir > 8:
						cdir = 1
					cc[nrows-i-1][j] = cdir

	# report elapsed time
	if console == True:
		elapsed_time = format_time_difference_hms(start_time, time.time())
		print("\r", elapsed_time, "\n", end="")

	return [ss, ff, dd, cc]

# get_regional_summary_csv compiles regional results for all hurricanes.
# Results are returned as a data frame of hurricane ids and maximum enhanced
# Fujita scale values.
# returns a data frame of summary values

def get_regional_summary_csv():
	# get current working directory
	cwd = os.getcwd()

	# read ids file
	ids_file = cwd + "/input/ids.csv"
	check_file_exists(ids_file)
	ii = pd.read_csv(ids_file)
	ii_rows = len(ii)

	# create lists for peak Fujita value across region
	hur_id_list = []
	efmax_list = []

	# record values for each hurricane
	for i in range(0, ii_rows):
		# get hurricane id
		hur_id = ii.hur_id[i]

		# read regional hurricane file in Geotiff format
		hur_tif_file = cwd + "/region-all/" + hur_id + ".tif"
		check_file_exists(hur_tif_file)
		hur_tif = rio.open(hur_tif_file)

		# get enhanced Fujita scale layer
		ff_array = hur_tif.read(2) # enhanced Fujita scale

		# update peak Fujita value
		efmax = np.amax(ff_array) - 2

		hur_id_list.append(hur_id)
		efmax_list.append(efmax)

	# create data frame
	kk_cols = ['hur_id', 'efmax']
	kk = pd.DataFrame(data=list(zip(hur_id_list, efmax_list)), columns=kk_cols)

	return kk

# get_regional_summary_tif compiles regional results for all hurricanes.
# Results are returned as a list of 7 raster arrays representing the maximum
# Fujita value and the number of storms for each Fujita value.
# returns a list of 7 raster arrays

def get_regional_summary_tif():
	# get current working directory
	cwd = os.getcwd()

	# read ids file
	ids_file = cwd + "/input/ids.csv"
	check_file_exists(ids_file)
	ii = pd.read_csv(ids_file)
	ii_rows = len(ii)

	# read land-water file
	land_water_file = cwd + "/input/land_water.tif"
	check_file_exists(land_water_file)
	land_water = rio.open(land_water_file)

	# get regional values
	nrows = land_water.height
	ncols = land_water.width

	# create arrays for peak wind speed & direction
	efm = np.zeros((nrows, ncols), dtype=np.int16)
	ef0 = np.zeros((nrows, ncols), dtype=np.int16)
	ef1 = np.zeros((nrows, ncols), dtype=np.int16)
	ef2 = np.zeros((nrows, ncols), dtype=np.int16)
	ef3 = np.zeros((nrows, ncols), dtype=np.int16)
	ef4 = np.zeros((nrows, ncols), dtype=np.int16)
	ef5 = np.zeros((nrows, ncols), dtype=np.int16)

	# record values for each hurricane
	for i in range(0, ii_rows):
		# get hurricane id
		hur_id = ii.hur_id[i]

		# read regional hurricane file in Geotiff format
		hur_tif_file = cwd + "/region-all/" + hur_id + ".tif"
		check_file_exists(hur_tif_file)
		hur_tif = rio.open(hur_tif_file)

		# get enhanced Fujita scale layer
		ff_array = hur_tif.read(2) # enhanced Fujita scale

		# update enhanced Fujita scale
		for j in range(0, nrows):
			for k in range(0, ncols):
				val = ff_array[j, k]

				if val > 0 and efm[j, k] < val:
					efm[j, k] = val
 
				if val == 2:
					ef0[j, k] = ef0[j, k] + 1
      
				elif val == 3:
					ef0[j, k] = ef0[j, k] + 1
					ef1[j, k] = ef1[j, k] + 1
      
				elif val == 4:
					ef0[j, k] = ef0[j, k] + 1
					ef1[j, k] = ef1[j, k] + 1
					ef2[j, k] = ef2[j, k] + 1

				elif val == 5:
					ef0[j, k] = ef0[j, k] + 1
					ef1[j, k] = ef1[j, k] + 1
					ef2[j, k] = ef2[j, k] + 1
					ef3[j, k] = ef3[j, k] + 1
        
				elif val == 6:
					ef0[j, k] = ef0[j, k] + 1
					ef1[j, k] = ef1[j, k] + 1
					ef2[j, k] = ef2[j, k] + 1
					ef3[j, k] = ef3[j, k] + 1
					ef4[j, k] = ef4[j, k] + 1          

				elif val == 7:
					ef0[j, k] = ef0[j, k] + 1
					ef1[j, k] = ef1[j, k] + 1
					ef2[j, k] = ef2[j, k] + 1
					ef3[j, k] = ef3[j, k] + 1
					ef4[j, k] = ef4[j, k] + 1          
					ef5[j, k] = ef5[j, k] + 1          

	return[efm, ef0, ef1, ef2, ef3, ef4, ef5]

# get_values_at_datetime returns a data frame of modeled values for the
# specified datetime. The data frame includes hurricane latitude (degrees)
# longitude (degrees), maximum sustained wind speed (meters/second), forward 
# speed (meters/second), and bearing (degrees).
#   hur_id - hurricane id
#   tt - data frame of track data
#   dt - datetime in the format YYYY-MM-DDThh:mm
# returns a data frame of modeled values

def get_values_at_datetime(hur_id, tt, dt):
	# interpolate hurricane location & max wind speed
	mm = interpolate_hurricane_location_max_wind(tt, time_step=1)	
	yr_list = mm[0]
	jd_list = mm[1]
	lat_list = mm[2]
	lon_list = mm[3]
	wmax_list = mm[4]
	dt_list = get_standard_date_time(yr_list, jd_list)
	
	# abort if no match
	if dt not in dt_list:
		sys.exit("Datetime not found")

	# interpolate hurricane speed & bearing
	mm = interpolate_hurricane_speed_bearing(tt, jd_list)
	spd_list = mm[0]
	bear_list = mm[1]

	# get values for specified datetime
	index = dt_list.index(dt)
  
	# extract values
	lat = lat_list[index]
	lon = lon_list[index]
	wmax = wmax_list[index]
	spd = spd_list[index]
	bear = bear_list[index]

	# create data fame of peak values
	pp_cols = ['lat', 'lon', 'wmax', 'spd', 'bear']

	pp = pd.DataFrame(columns=pp_cols)

	pp.loc[0] = [lat, lon, wmax, spd, bear]

	return pp

# get_track_lat_lon returns a data frame of track data for the specified hurricane
# if the maximum enhanced Fujita value exceeds a specified value.
#   hur_id - hurricane id
#   fuj_min - minimum enhanced Fujita value
#   tt - data frame of track data (all hurricanes)
#   kk - data frame of summary data
# returns a data frame of track data (this hurricane)

def get_track_lat_lon(hur_id, fuj_min, tt, kk):
	xx = kk.loc[kk.hur_id == hur_id, "efmax"]

	if len(xx) == 0:
		return ""

	efmax = int(kk.loc[kk.hur_id == hur_id, "efmax"])

	if efmax < fuj_min:
		return ""
  
	else:
		xx = tt.loc[tt.hur_id == hur_id, ]
		return xx 


### UTILITY FUNCTIONS #####################################

# hurrecon_set_path sets the path for the current set of model runs.
#   hur_path - path for current model runs
#   console - whether to display messages in console
# no return value

def hurrecon_set_path(hur_path, console=True):
	if hur_path == "":
		sys.exit("Need to enter a path")

	elif os.path.exists(hur_path) == False:
		sys.exit("Path does not exist")

	os.chdir(hur_path)

	if console == True:
		print("Path set to", hur_path)

# hurrecon_create_land_water creates a land-water raster file in GeoTiff 
# format from vector boundary files in shapefile format. The land-water file
# (land_water.tif) is assumed to be aligned with lines of latitude and 
# longitude.  Boundary files are assumed to be named boundary.* on the vector 
# subdirectory. The land-water file is created on the input subdirectory.
#   nrows - number of rows
#   ncols - number of columns
#   xmn - minimum longitude (degrees)
#   xmx - maximum longitude (degrees)
#   ymn - minimum latitude (degrees)
#   ymx - maximum latitude (degrees)
#   console - whether to display messages in console
# no return value

def hurrecon_create_land_water(nrows, ncols, xmn, xmx, ymn, ymx, console=True):
	# get current working directory
	cwd = os.getcwd()

	# open boundaries file
	boundaries_file = cwd + "/vector/boundaries.shp"
	shapefile = fiona.open(boundaries_file, "r")
	shapes = [feature["geometry"] for feature in shapefile]
	shapefile.close()

	# create GeoTiff file
	lw_tif_file = cwd + "/input/land_water.tif"

	transform = rio.transform.from_bounds(west=xmn, south=ymn, east=xmx, north=ymx, 
		width=ncols, height=nrows)

	# fill with zeros & save to file
	aa = np.zeros((nrows, ncols), dtype=np.int16)

	lw_tif = rio.open(lw_tif_file, 'w', driver='GTiff', height=nrows, width=ncols, 
		count=1, dtype='int16', crs='+proj=latlong', transform=transform, 
		nodata=None)

	lw_tif.write(aa, 1)
	lw_tif.close()

	# mask with boundaries file (water = -99)
	lw_tif = rio.open(lw_tif_file, 'r', driver='GTiff', height=nrows, width=ncols, 
		nt=1, dtype='int16', crs='+proj=latlong', transform=transform, nodata=None)

	out_image, out_transform = mask(lw_tif, shapes, nodata=-99, crop=False)
	out_meta = lw_tif.meta
	
	out_meta.update({"driver": "GTiff",	"height": nrows, "width": ncols, 
		"transform": out_transform})

	lw_tif.close()

	# reclasssify (water = 1, land = 2)
	for i in range(0, nrows):
		for j in range(0, ncols):
			if out_image[0][i][j] == -99:
				out_image[0][i][j] = 1
			else:
				out_image[0][i][j] = 2

	# save to file
	dest = rio.open(lw_tif_file, "w", **out_meta)
	dest.write(out_image)
	dest.close()

	# calculate cell dimensions in kilometers
	lat_avg = (ymx + ymn)/2

	cell_height = 111*(ymx-ymn)/nrows
	cell_width = 111*(xmx-xmn)*math.cos(lat_avg*math.pi/180)/ncols

	if console == True:
		print("Cell height =", round(cell_height) , "kilometers")
		print("Cell width  =", round(cell_width), "kilometers\n")

# hurrecon_reformat_hurdat2 reformats a HURDAT2 file from the National 
# Hurricane Center for use with the HURRECON model. The input file is assumed
# to be in space-delimited text format. The output file (hurdat2_tracks.csv)
# contains full track information for each hurricane plus columns for standard 
# datetime and Julian day with fraction.
#   hurdat2_file - name of HURDAT2 file
#   path - optional path for input & output files
#   console - whether to display messages in console
# no return value

def hurrecon_reformat_hurdat2(hurdat2_file, path="", console=True):
	# output files
	track_file = "hurdat2_tracks.csv"

	if path != "":
		if path[len(path)-1] != "/":
			path = path + "/"
		hurdat2_file = path + hurdat2_file
		track_file = path + track_file

	# read hurdat2 file
	file_in = open(hurdat2_file, 'r') 
	hurdat = file_in.readlines() 
	nlines = len(hurdat)

	# close hurdat2 file
	file_in.close()

	# initialize lists
	tracks_hur_id_list = []
	tracks_name_list = []
	tracks_date_time_list = []
	tracks_jd_list = []
	tracks_status_list = []
	tracks_latitude_list = []
	tracks_longitude_list = []
	tracks_wind_max_list = []

	# current line number in hurdat
	line_num = -1

	# process file
	while line_num < nlines - 1:
		# get hurricane id, name, and number of positions
		line_num = line_num + 1
		row = hurdat[line_num].split(",") 
		hur_id = row[0].strip()
		name = row[1].strip()
		positions = int(row[2].strip())

		# process observations
		for i in range(0, positions):
			line_num = line_num + 1
			row = hurdat[line_num].split(",")
			date = row[0].strip()
			time = row[1].strip()
			status = row[3].strip()

			lat = row[4].strip()
			latitude = float(lat[0:len(lat)-1])

			lon = row[5].strip()
			longitude = -float(lon[0:len(lon)-1])

			wind_max = float(row[6].strip())

			# convert knots to meters per second
			wind_max = round(0.514444 * wind_max, 1)

			# get peak wind
			if i == 0:
				wind_peak = wind_max
			else:
				if wind_peak < wind_max:
					wind_peak = wind_max

			# get datetime
			year = date[0:4]
			month = date[4:6]
			day = date[6:8]
			hour = int(time[0:2])
			minute = int(time[2:4])

			date_time = year + "-" + month + "-" + day + "T" + str(hour).zfill(2) + ":" + str(minute).zfill(2)

			# get Julian date
			jd = float(calculate_julian_day(int(year), int(month), int(day), hour, minute))

			# store values in tracks lists
			tracks_hur_id_list.append(hur_id)
			tracks_name_list.append(name)
			tracks_date_time_list.append(date_time)
			tracks_jd_list.append(jd)
			tracks_status_list.append(status)
			tracks_latitude_list.append(latitude)
			tracks_longitude_list.append(longitude)
			tracks_wind_max_list.append(wind_max)
	
		# report progress
		if console == True:
			x = round(line_num*100/nlines)
			if x % 10 == 0:
				print("          ", end="")
				print("\r", x, "%", end="")

	# create data frame
	tracks_cols = ['hur_id', 'name', 'date_time', 'jd', 'status', 'latitude', 'longitude', 'wind_max']

	tracks= pd.DataFrame(data=list(zip(tracks_hur_id_list, tracks_name_list, tracks_date_time_list, 
		tracks_jd_list, tracks_status_list, tracks_latitude_list, tracks_longitude_list, tracks_wind_max_list)), 
		columns=tracks_cols)

	# save to file
	tracks.to_csv(track_file, index=False)

	# get number of storms
	ii = set(tracks_hur_id_list)

	# display number of storms
	if console == True:
		print("\nNumber of storms = ", len(ii))
		print("Number of observations = ", len(tracks))

# hurrecon_extract_tracks extracts track data from an input track file
# (input_tracks.csv) created from HURDAT2 using hurrecon_reformat_hurdat2
# or created from other sources using the same file structure. The geographic 
# window used to select hurricanes is set by the land-water file and optionally
# extended by the margin parameter. Selection begins by identifying all positions
# in the window where winds reach or exceed hurricane speed (33 meters/second). 
# If at least one such position exists, the track is extended to include one 
# position before and one position after the first and last hurricane position, 
# if possible. If the resulting track contains at least two positions and the 
# maximum sustained wind speed equals or exceeds wind_min, the track is included.
# For included storms, summary data are written to ids.csv, track data are written 
# to tracks.csv, and track data for all positions are written to tracks_all.csv.
#   margin - optional extension of the geographic window set by the
#     land-water file (degrees)
#   wind_min - minimum value of maximum sustained wind speed 
#     (meters/second)
#   status - whether to limit search to storms with hurricane status
#   console - whether to display messages in console
# no return value

def hurrecon_extract_tracks(margin=0, wind_min=33, status=True, console=True):
	# get current working directory
	cwd = os.getcwd()

	# output files
	ids_file = cwd + "/input/ids.csv"
	track_file = cwd + "/input/tracks.csv"
	track_all_file = cwd + "/input/tracks_all.csv"

	# read input tracks file
	input_track_file = cwd + "/input/input_tracks.csv"
	check_file_exists(input_track_file)
	tt = pd.read_csv(input_track_file)
	tt_rows = len(tt)

	# get ids
	xx = tt[["hur_id", "name"]]
	zz = xx.drop_duplicates()
	ii_hur_id = list(zz.hur_id)
	ii_name = list(zz.name)
	ii_rows = len(ii_hur_id)
	
	# read land-water file
	land_water_file = cwd + "/input/land_water.tif"
	check_file_exists(land_water_file)
	land_water = rio.open(land_water_file)

	lat_min = land_water.bounds.bottom - margin
	lat_max = land_water.bounds.top + margin
	lon_min = land_water.bounds.left - margin
	lon_max = land_water.bounds.right + margin

	land_water.close()

	# initialize lists
	ids_hur_id_list = []
	ids_name_list = []
	ids_positions_list = []
	ids_wind_peak_list = []

	tracks_hur_id_list = []
	tracks_name_list = []
	tracks_date_time_list = []
	tracks_jd_list = []
	tracks_status_list = []
	tracks_latitude_list = []
	tracks_longitude_list = []
	tracks_wind_max_list = []

	tracks_all_hur_id_list = []
	tracks_all_date_time_list = []
	tracks_all_latitude_list = []
	tracks_all_longitude_list = []

	# process files
	for i in range(0, ii_rows):
		# get hurricane id & name
		hur_id = ii_hur_id[i]
		name = ii_name[i]

		if status == True:
			index = np.where((tt.hur_id == hur_id) & (tt.latitude >= lat_min) & (tt.latitude <= lat_max) & (tt.longitude >= lon_min) & (tt.longitude <= lon_max) & (tt.wind_max >= 33) & (tt.status == "HU"))[0].tolist()
		else:
			index = np.where((tt.hur_id == hur_id) & (tt.latitude >= lat_min) & (tt.latitude <= lat_max) & (tt.longitude >= lon_min) & (tt.longitude <= lon_max) & (tt.wind_max >= 33))[0].tolist()

		index_all = np.where(tt.hur_id == hur_id)[0].tolist()

		# get start & end position
		if len(index) > 0:
			start_index = index[0]

			if start_index > 0:
				if tt.hur_id[start_index - 1] == hur_id:
					start_index = start_index - 1
 
			end_index = index[len(index) - 1]

			if end_index < tt_rows - 1:
				if tt.hur_id[end_index + 1] == hur_id:
					end_index = end_index + 1

			# subset by hurricane name and bounding region
			xx = tt.loc[(tt.hur_id == hur_id) & (tt.index >= start_index) & (tt.index <= end_index)]
			positions = len(xx)
			wind_peak = max(xx.wind_max)
		
			zz = tt.loc[index_all, ]

			# store id & tracks if at least 2 positions & exceeds minimum wind speed
			if positions > 1 and wind_peak >= wind_min:
				ids_hur_id_list.append(hur_id)
				ids_name_list.append(name)
				ids_positions_list.append(positions)
				ids_wind_peak_list.append(wind_peak)

				tracks_hur_id_list.extend(xx.hur_id)
				tracks_name_list.extend(xx.name)
				tracks_date_time_list.extend(xx.date_time)
				tracks_jd_list.extend(xx.jd)
				tracks_status_list.extend(xx.status)
				tracks_latitude_list.extend(xx.latitude)
				tracks_longitude_list.extend(xx.longitude)
				tracks_wind_max_list.extend(xx.wind_max)

				tracks_all_hur_id_list.extend(zz.hur_id)
				tracks_all_date_time_list.extend(zz.date_time)
				tracks_all_latitude_list.extend(zz.latitude)
				tracks_all_longitude_list.extend(zz.longitude)

		# report progress
		if console == True:
			x = round(i*100/ii_rows)
			if x % 10 == 0:
				print("          ", end="")
				print("\r", x, "%", end="")

	# create data frames
	ids_cols = ['hur_id', 'name', 'positions', 'wind_peak']
	ids = pd.DataFrame(data=list(zip(ids_hur_id_list, ids_name_list, ids_positions_list, ids_wind_peak_list)),
		columns=ids_cols)

	tracks_cols = ['hur_id', 'name', 'date_time', 'jd', 'status', 'latitude', 'longitude', 'wind_max']

	tracks = pd.DataFrame(data=list(zip(tracks_hur_id_list, tracks_name_list, tracks_date_time_list, 
		tracks_jd_list, tracks_status_list, tracks_latitude_list, tracks_longitude_list, tracks_wind_max_list)), 
		columns=tracks_cols)

	tracks_all_cols = ['hur_id', 'date_time', 'latitude', 'longitude']

	tracks_all = pd.DataFrame(data=list(zip(tracks_all_hur_id_list, tracks_all_date_time_list, 
		tracks_all_latitude_list, tracks_all_longitude_list)), columns=tracks_all_cols)

	# save to file
	ids.to_csv(ids_file, index=False)
	tracks.to_csv(track_file, index=False)
	tracks_all.to_csv(track_all_file, index=False)

	# display number of storms
	if console == True:
		print("\nNumber of storms = ", len(ids))
		print("Number of observations = ", len(tracks))


### MODELING FUNCTIONS ####################################

# hurrecon_model_site calculates wind speed (meters/second), gust speed 
# (meters/second), wind direction (degrees), and enhanced Fujita scale wind 
# damage for a given hurricane and site. If width is True, the radius of 
# maximum wind (rmw) and scaling parameter (s_par) for this hurricane are 
# used; otherwise values for ALL are used. If save is True, results are 
# saved to a CSV file on the site subdirectory; otherwise results are returned
# as a data frame.
#   hur_id - hurricane id
#   site_name - name of site
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   save - whether to save results to a CSV file
#   console - whether to display messages in console
# returns a data frame of results if save is False

def hurrecon_model_site(hur_id, site_name, width=False, time_step=1, save=True, 
	console=True):

	# record total elapsed time
	start_time = time.time()

    # get current working directory
	cwd = os.getcwd()

	# read sites file
	sites = read_site_file(site_name)
	site_latitude = sites[0]
	site_longitude = sites[1]
	cover_type = sites[2]

	# read parameters file
	pars = read_parameter_file(hur_id, width)
	rmw = pars[0]
	s_par = pars[1]

	# get fixed parameters
	fixed = get_fixed_model_parameters(cover_type)
	asymmetry_factor = fixed[0]
	inflow_angle = fixed[1]
	friction_factor = fixed[2]
	gust_factor = fixed[3]

	# read hurricane track file
	tt = read_hurricane_track_file(hur_id)

	# interpolate hurricane location & max wind speed
	mm = interpolate_hurricane_location_max_wind(tt, time_step)
	yr_list = mm[0]
	jd_list = mm[1]
	lat_list = mm[2]
	lon_list = mm[3]
	wmax_list = mm[4]

	# get number of rows
	mm_rows = len(mm[0])

	# interpolate hurricane speed & bearing
	mm = interpolate_hurricane_speed_bearing(tt, jd_list)
	spd_list = mm[0]
	bear_list = mm[1]

	# calculate range & bearing from site to hurricane center
	mm = calculate_site_range_bearing(lat_list, lon_list, site_latitude, site_longitude)
	srange_list = mm[0]
	sbear_list = mm[1]

	# calculate wind speed, wind direction & enhanced Fujita scale at site
	mm = calculate_wind_speed_direction(sbear_list, srange_list, lat_list, bear_list, 
		spd_list, wmax_list, inflow_angle, rmw, s_par, asymmetry_factor, friction_factor, 
		gust_factor)
	wspd_list = mm[0]
	gspd_list = mm[1]
	wdir_list = mm[2]
	ef_list   = mm[3]

	# get standard date & time
	dt_list = get_standard_date_time(yr_list, jd_list)

	# get constant parameters
	rmw_list = [rmw] * mm_rows
	spar_list = [s_par] * mm_rows

	# create data frame
	mm_cols = ["date_time", "year", "jd", "latitude", "longitude", "wind_max", "hur_bear",
		"hur_spd", "site_bear", "site_range", "rmw", "s_par", "wind_dir", "wind_spd",
		"gust_spd", "ef_sca"]

	mm = pd.DataFrame(data=list(zip(dt_list, yr_list, jd_list, lat_list, lon_list, wmax_list, 
		bear_list, spd_list, sbear_list, srange_list, rmw_list, spar_list, wdir_list, wspd_list, 
		gspd_list, ef_list)), columns=mm_cols)

	# display total elapsed time
	if console == True:
		print(format_time_difference_ms(start_time, time.time()), " ms")

	# output
	if save == True:
		# save modeled data to CSV file
		modeled_file = cwd + "/site/" + hur_id + " " + site_name + ".csv"
		mm.to_csv(modeled_file, index=False)
	
		if console == True:
			print("\rSaving to", modeled_file)
	else:
		# return modeled data as data frame
		return mm

# hurrecon_model_site_all creates a table of peak values for all hurricanes
# for a given site. If width is True, the radius of maximum wind (rmw) and 
# scaling parameter (s_par) for the given hurricane are used; otherwise values
# for ALL are used. If save is True, results are saved to a CSV file on the 
# site-all subdirectory; otherwise results are returned as a data frame.
#   site_name - name of site
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   save - whether to save results to a CSV file
#   console - whether to display messages in console
# returns a data frame of results if save is False

def hurrecon_model_site_all(site_name, width=False, time_step=1, save=True, 
	console=True):

	# get current working directory
	cwd = os.getcwd()

	# read ids file
	ids_file = cwd + "/input/ids.csv"
	check_file_exists(ids_file)
	ii = pd.read_csv(ids_file)
	ii_rows = len(ii)

	# initialize lists
	snam_list = []
	hid_list = []
	dt_list = []
	wdir_list = []
	wspd_list = []
	gspd_list = []
	ef_list = []
	ef0_list = []
	ef1_list = []
	ef2_list = []
	ef3_list = []
	ef4_list = []
	ef5_list = []

	# record total elasped time
	start_time = time.time()

	# get peak values for each hurricane
	for i in range(0, ii_rows):
		# get hurricane name
		hur_id = ii.hur_id[i]

		# get modeled output
		mm = hurrecon_model_site(hur_id, site_name, width, time_step, save=False,
			console=False)

		# get peak values
		pk = get_peak_values(hur_id, site_name, mm)
		
		snam_list.append(pk.site_name[0])
		hid_list.append(pk.hur_id[0])
		dt_list.append(pk.date_time[0])
		wdir_list.append(pk.wind_dir[0])
		wspd_list.append(pk.wind_spd[0])
		gspd_list.append(pk.gust_spd[0])
		ef_list.append(pk.ef_sca[0])
		ef0_list.append(pk.ef0[0])
		ef1_list.append(pk.ef1[0])
		ef2_list.append(pk.ef2[0])
		ef3_list.append(pk.ef3[0])
		ef4_list.append(pk.ef4[0])
		ef5_list.append(pk.ef5[0])

		# report progress
		if console == True:
			x = round(i*100/ii_rows)
			if x % 10 == 0:
				print("          ", end="")
				print("\r", x, "%", end="")

	if console == True:
		elapsed_time = format_time_difference_hms(start_time, time.time())
		print("\r", elapsed_time, "\n", end="")

	# create data frame for peak values
	peak_values_cols = ['site_name', 'hur_id', 'date_time', 'wind_dir', 'wind_spd',
		'gust_spd', "ef_sca", "ef0", "ef1", "ef2", "ef3", "ef4", "ef5"]

	peak_values = pd.DataFrame(data=list(zip(snam_list, hid_list, dt_list, wdir_list, wspd_list, 
		gspd_list, ef_list, ef0_list, ef1_list, ef2_list, ef3_list, ef4_list, ef5_list)), 
		columns=peak_values_cols)

	# output
	if (save == True):
		# save modeled data to CSV file
		site_peak_file = cwd + "/site-all/" + site_name + " Peak Values.csv"
		peak_values.to_csv(site_peak_file, index=False)
	
		if console == True:
			print("Saving to", site_peak_file)
	else:
		# return modeled data as data frame
		return peak_values

# hurrecon_model_region calculates peak wind speed (meters/second), peak 
# enhanced Fujita scale, peak wind direction (degrees), peak cardinal wind 
# direction, and duration of EF0, EF1, EF2, EF3, EF4, and EF5 winds (minutes)
# for a given hurricane over a region. If width is TRUE, the radius of maximum 
# wind (rmw) and scaling parameter (s_par) for the given hurricane are used; 
# otherwise values for ALL are used. If time_step is NULL, the time step is 
# calculated. If water is FALSE, results are calculated for land areas only. 
# If save is TRUE, results are saved as a GeoTiff file on the region subdirectory.
#   hur_id - hurricane id
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   water - whether to caculate results over water
#   save - whether to save results to a GeoTiff file
#   console - whether to display messages in console
# returns a list of 10 raster arrays if save is False

def hurrecon_model_region(hur_id, width=False, time_step="", water=False, save=True, 
	console=True):
	
	# get current working directory
	cwd = os.getcwd()

	# get time step if necessary
	if time_step == "":
		time_step = get_time_step()

	if console == True:
		print("Time step =", time_step, "minutes")
 
 	# read hurricane track file
	tt = read_hurricane_track_file(hur_id)

	# interpolate hurricane location & max wind speed
	mm = interpolate_hurricane_location_max_wind(tt, time_step)
	jd_list = mm[1]
	lat_list = mm[2]
	lon_list = mm[3]
	wmax_list = mm[4]

	# interpolate hurricane speed & bearing
	mm = interpolate_hurricane_speed_bearing(tt, jd_list)
	spd_list = mm[0]
	bear_list = mm[1]

	# get modeled values over region
	peak_list = get_regional_peak_wind(hur_id, lat_list, lon_list, wmax_list,
		bear_list, spd_list, width, time_step, water, console)

	# output
	if (save == True):
		# read land-water file
		land_water_file = cwd + "/input/land_water.tif"
		check_file_exists(land_water_file)
		land_water = rio.open(land_water_file)

		# open GeoTiff file
		hur_tif_file = cwd + "/region/" + hur_id + ".tif"

		profile = land_water.profile
		profile.update(dtype='int16', nodata=-9999, count=10)
        
		hur_tif = rio.open(hur_tif_file, 'w', **profile)

		hur_tif.write(peak_list[0], 1)
		hur_tif.write(peak_list[1], 2)
		hur_tif.write(peak_list[2], 3)
		hur_tif.write(peak_list[3], 4)
		hur_tif.write(peak_list[4], 5)
		hur_tif.write(peak_list[5], 6)
		hur_tif.write(peak_list[6], 7)
		hur_tif.write(peak_list[7], 8)
		hur_tif.write(peak_list[8], 9)
		hur_tif.write(peak_list[9], 10)

		hur_tif.close()

		if console == True:
			print("Saving to", hur_tif_file)

	else:
		# return modeled values as a list of 6 raster arrays
		return peak_list

# hurrecon_model_region_dt calculates wind speed (meters/second), enhanced
# Fujita scale, wind direction (degrees), and cardinal wind direction for a
# given hurricane over a region at a specified datetime. If width is
# True, the radius of maximum wind (rmw) and scaling parameter (s_par) for 
# this hurricane are used; otherwise values for ALL are used . If water is 
# False, results are calculated for land areas only. If save is True, results
# are saved as a GeoTiff file on the region-dt subdirectory; otherwise results
# are returned as a list of 4 raster arrays.
#   hur_id - hurricane id
#   dt - datetime in the format YYYY-MM-DDThh:mm
#   width - whether to use width parameters for the specified hurricane
#   water - whether to caculate results over water
#   save - whether to save results to a GeoTiff file
#   console - whether to display messages in console
# returns a list of 4 rasters if save is False

def hurrecon_model_region_dt(hur_id, dt, width=False, water=False, save=True, 
	console=True):

	# get current working directory
	cwd = os.getcwd()

 	# read hurricane track file
	tt = read_hurricane_track_file(hur_id)

	# get values for specified datetime
	pp = get_values_at_datetime(hur_id, tt, dt)

	# get modeled values over region
	datetime_list = get_regional_datetime(hur_id, pp.lat[0], pp.lon[0], pp.wmax[0], 
		pp.bear[0], pp.spd[0], width, water, console)

	# output
	if (save == True):
		# read land-water file
		land_water_file = cwd + "/input/land_water.tif"
		check_file_exists(land_water_file)
		land_water = rio.open(land_water_file)

		# open GeoTiff file
		dt2 = dt.replace(":", "")
		hur_tif_file = cwd + "/region-dt/" + hur_id + " " + dt2 + ".tif"

		profile = land_water.profile
		profile.update(dtype='int16', nodata=-9999, count=6)
        
		hur_tif = rio.open(hur_tif_file, 'w', **profile)

		hur_tif.write(datetime_list[0], 1)
		hur_tif.write(datetime_list[1], 2)
		hur_tif.write(datetime_list[2], 3)
		hur_tif.write(datetime_list[3], 4)

		hur_tif.close()

		if console == True:
			print("Saving to", hur_tif_file)

	else:
		# return modeled values as a list of 4 raster arrays
		return datetime_list

# hurrecon_model_region_all calculates peak wind speed (meters/second), 
# enhanced Fujita scale, wind direction (degrees), cardinal wind direction, 
# duration of gale winds (minutes), and duration of hurricane winds (minutes) 
# over a region for all hurricanes. If width is True, the radius of maximum 
# wind (rmw) and scaling parameter (s_par) for the given hurricane are used;
# otherwise values for ALL are used. If no value is provided for time step, 
# the time step is calculated. If water is False, results are calculated for 
# land areas only. Results for each hurricane are saved in a GeoTiff file on 
# the region-all subdirectory. Summary results for all hurricanes (summary.tif,
# summary.csv) are also calculated and saved to the region-all subdirectory.
# If returns is True, summary values are returned.
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   water - whether to calculate results over water
#   console - whether to display messages in console
#   returns - whether to return summary values
# returns a list containing a data frame and a list of raster arrays

def hurrecon_model_region_all(width=False, time_step="", water=False, 
	console=True, returns=False):
	
	# get current working directory
	cwd = os.getcwd()

	# get time step if necessary
	if time_step == "":
		time_step = get_time_step()
	
	if console == True:
		print("Time step =", time_step, "minutes")

	# read ids file
	ids_file = cwd + "/input/ids.csv"
	check_file_exists(ids_file)
	ii = pd.read_csv(ids_file)
	ii_rows = len(ii)

	# read land-water file
	land_water_file = cwd + "/input/land_water.tif"
	check_file_exists(land_water_file)
	land_water = rio.open(land_water_file)

	profile = land_water.profile
	profile.update(dtype='int16', nodata=-9999, count=7)

	# record total elasped time
	start_time = time.time()

	# get regional estimate for each hurricane
	for i in range(0, ii_rows):
		# get hurricane name
		hur_id = ii.hur_id[i]

		# report progress
		if console == True:
			x = round(i*100/ii_rows)
			print("          ", end="")
			print("\r", x, "%", end="")

		# generate & save regional results
		peak_list = hurrecon_model_region(hur_id, width, time_step, water, 
			save=False, console=False)

		# open GeoTiff file
		hur_tif_file = cwd + "/region-all/" + hur_id + ".tif"
		hur_tif = rio.open(hur_tif_file, 'w', **profile)

		hur_tif.write(peak_list[0], 1)
		hur_tif.write(peak_list[1], 2)
		hur_tif.write(peak_list[2], 3)
		hur_tif.write(peak_list[3], 4)
		hur_tif.write(peak_list[4], 5)
		hur_tif.write(peak_list[5], 6)

		hur_tif.close()

	# get & save summary.csv file
	kk = get_regional_summary_csv()
	peak_file = cwd + "/region-all/summary.csv"
	kk.to_csv(peak_file, index=False)

	# get & save summary.tif file
	sum_list = get_regional_summary_tif()
	sum_tif_file = cwd + "/region-all/summary.tif"
	sum_tif = rio.open(sum_tif_file, 'w', **profile)

	sum_tif.write(sum_list[0], 1)
	sum_tif.write(sum_list[1], 2)
	sum_tif.write(sum_list[2], 3)
	sum_tif.write(sum_list[3], 4)
	sum_tif.write(sum_list[4], 5)
	sum_tif.write(sum_list[5], 6)
	sum_tif.write(sum_list[6], 7)

	sum_tif.close()

	if console == True:
		# display total elapsed time
		elapsed_time = format_time_difference_hms(start_time, time.time())
		print("\r", elapsed_time, "\n", end="")

		# display where results are saved
		reg_all_dir = cwd + "/region-all/"
		print("Saving to", reg_all_dir)

	# return a list of summary results
	if returns == True:
		return [kk, sum_list]

### SUMMARIZING FUNCTIONS #################################

# hurrecon_summarize_land_water displays information about the current 
# land-water file (land_water.tif) in the console.
#   console - whether to display messages in console
# returns a string containing summary information if console is False

def hurrecon_summarize_land_water(console=True):
	# get current working directory
	cwd = os.getcwd()

	# read land-water file
	land_water_file = cwd + "/input/land_water.tif"
	check_file_exists(land_water_file)
	land_water = rio.open(land_water_file)

	# get regional values
	nrows = land_water.height
	ncols = land_water.width

	ymn = land_water.bounds.bottom
	ymx = land_water.bounds.top

	xmn = land_water.bounds.left
	xmx = land_water.bounds.right

	# calculate cell dimensions in kilometers
	lat_avg = (ymx + ymn)/2

	cell_height = 111*(ymx-ymn)/nrows
	cell_width = 111*(xmx-xmn)*math.cos(lat_avg*math.pi/180)/ncols

	# get default time step

	time_step = get_time_step()

	st = "Rows: " + str(nrows) + "  Columns: " + str(ncols) + "\n"
	st = st + "Latitude: " + str(ymn) + " to " + str(ymx) + " degrees" + "\n"
	st = st + "Longitude: " + str(xmn) + " to " + str(xmx) + " degrees" + "\n"
	st = st + "Cell height: " + str(round(cell_height)) + " kilometers" + "\n"
	st = st + "Cell width: " + str(round(cell_width)) + " kilometers" + "\n"
	st = st + "Time Step: " + str(time_step) + " minutes" + "\n"

	if console == True:
		print(st)
	else:
		return st

# hurrecon_summarize_tracks displays information about the current ids 
# file (ids.csv) in the console.
#   console - whether to display messages in console
# returns a string containing summary information if console is False

def hurrecon_summarize_tracks(console=True):
	# get current working directory
	cwd = os.getcwd()

	# read ids file
	ids_file = cwd + "/input/ids.csv"
	check_file_exists(ids_file)
	ii = pd.read_csv(ids_file)
	ii_rows = len(ii)

	positions_total = sum(ii.positions)

	wind_peak_min = min(ii.wind_peak)
	wind_peak_max = max(ii.wind_peak)

	st = "Number of storms = " + str(ii_rows) + "\n"
	st = st + "Number of positions = " + str(positions_total) + "\n"
	st = st + "Minimum peak wind = " + str(wind_peak_min) + " m/s" + "\n"
	st = st + "Maximum peak wind = " + str(wind_peak_max) + " m/s" + "\n"

	if console == True:
		print(st)
	else:
		return st

# hurrecon_summarize_site displays peak values for a given hurricane
# and site in the console.
#   hur_id - hurricane id
#   site_name - name of site
#   console - whether to display messages in console
# returns a string containing summary information if console is False

def hurrecon_summarize_site(hur_id, site_name, console=True):
	# get current working directory
	cwd = os.getcwd()

	# read data
	modeled_name = hur_id + " " + site_name
	modeled_file = cwd + "/site/" + modeled_name + ".csv"
	check_file_exists(modeled_file)
	mm = pd.read_csv(modeled_file)

	# get peak values
	pk = get_peak_values(hur_id, site_name, mm)

	# print peak values
	st = modeled_name + "\n"

	st = st + "PEAK: " + str(pk.date_time[0]) + "\n"
	st = st + "Wind dir: " + str(int(round(pk.wind_dir[0]))) + " deg" + "\n"
	st = st + "Wind spd: " + str(int(round(pk.wind_spd[0]))) + " m/s" + "\n"
	st = st + "Gust spd: " + str(int(round(pk.gust_spd[0]))) + " m/s" + "\n"

	if pk.ef0[0] > 0:
		st = st + "EF0: " + str(round(pk.ef0[0], 1)) + " hours" + "\n"
	if pk.ef1[0] > 0:
		st = st + "EF1: " + str(round(pk.ef1[0], 1)) + " hours" + "\n"
	if pk.ef2[0] > 0:
		st = st + "EF2: " + str(round(pk.ef2[0], 1)) + " hours" + "\n"
	if pk.ef3[0] > 0:
		st = st + "EF3: " + str(round(pk.ef3[0], 1)) + " hours" + "\n"
	if pk.ef4[0] > 0:
		st = st + "EF4: " + str(round(pk.ef4[0], 1)) + " hours" + "\n"
	if pk.ef5[0] > 0:
		st = st + "EF5: " + str(round(pk.ef5[0], 1)) + " hours" + "\n"

	if console == True:
		print(st)
	else:
		return st

### PLOTTING FUNCTIONS ####################################

# hurrecon_plot_site creates a time-series plot (wind speed, gust 
# speed, or wind direction as a function of datetime) or a scatter 
# plot (wind speed or gust speed as a function of wind direction) 
# for a given hurricane and site. Optional start and end datetimes 
# may be specified. X-variables: datetime or wind_direction. 
# Y-variables: wind_speed, gust_speed, or wind_direction.
#   hur_id - hurricane id
#   site_name - name of site
#   start_datetime - optional start datetime in format YYYY-MM-DD hh:mm
#   end_datetime - optional end datetime in format YYYY-MM-DD hh:mm
#   xvar - dependent variable
#   yvar - independent variable
#   adjust - whether to subtract 360 degrees from wind directions
#      greater than 180 degrees in scatter plot
#   legend_loc - legend location
#   main_title - optional title
# no return value

def hurrecon_plot_site(hur_id, site_name, start_datetime='', end_datetime='', 
	xvar="datetime", yvar="wind_speed", adjust=False, legend_loc="upper right",
	main_title=""):
	
	# register matplotlib converters
	from pandas.plotting import register_matplotlib_converters
	register_matplotlib_converters()

	# get current working directory
	cwd = os.getcwd()

	# get enhanced Fujita wind speeds & colors
	ef = get_fujita_wind_speeds()
	ef_col = get_fujita_colors()

	# check legend location
	if (not(check_legend_location(legend_loc))):
		legend_loc = "upper right"

	# read data
	modeled_file = cwd + "/site/" + hur_id + " " + site_name + ".csv"
	check_file_exists(modeled_file)
	mm = pd.read_csv(modeled_file)
	mm_rows = len(mm)

	# add datetime
	dt_list = []

	for i in range(0, mm_rows):
		year = int(mm.date_time[i][0:4])
		month = int(mm.date_time[i][5:7])
		day = int(mm.date_time[i][8:10])
		hour = int(mm.date_time[i][11:13])
		minute = int(mm.date_time[i][14:16])

		dt_list.append(dt.datetime(year, month, day, hour, minute))

	mm['dt'] = dt_list

	# x variable
	if xvar == "datetime":
		plot_type = "time_series"
		x_var = "dt"
		x_label = "Datetime (UTC)"

	elif xvar == "wind_direction":
		plot_type = "scatter_plot"
		x_var = "wind_dir"
		x_label = "Wind Direction (deg)"

	else:
		sys.exit("xvar must be datetime or wind_direction")

	# y variable
	if yvar == "wind_speed":
		y_var = "wind_spd"
		y_label = "Wind Speed (m/s)"

	elif yvar == "gust_speed":
		y_var = "gust_spd"
		y_label = "Gust Speed (m/s)"

	elif yvar == "wind_direction":
		y_var = "wind_dir"
		y_label = "Wind Direction (deg)"

	else:
		sys.exit("yvar must be wind_speed, gust_speed, or wind_direction")

	# adjust wind direction
	if plot_type == "scatter_plot" and adjust == True:
		wdir_list = []

		for i in range(0, mm_rows):
			wdir = mm.wind_dir[i]
			if wdir > 180:
				wdir = wdir - 360

			wdir_list.append(wdir)

		mm['wind_dir2'] = wdir_list
		x_var = "wind_dir2"

	# subset by data range
	if start_datetime != "":
		sdate = start_datetime
	else:
		sdate = mm.dt[0]

	if end_datetime != "":
		edate = end_datetime
	else:
		edate = mm.dt[mm_rows-1]

	mm_plot = mm.loc[(mm.dt >= sdate) & (mm.dt <= edate), ]

	# subset by enhanced Fujita scale
	mm_plot_efx = mm_plot.loc[(mm_plot.gust_spd < ef[0]), : ]
	mm_plot_ef0 = mm_plot.loc[(mm_plot.gust_spd >= ef[0]) & (mm_plot.gust_spd < ef[1]), : ]
	mm_plot_ef1 = mm_plot.loc[(mm_plot.gust_spd >= ef[1]) & (mm_plot.gust_spd < ef[2]), : ]
	mm_plot_ef2 = mm_plot.loc[(mm_plot.gust_spd >= ef[2]) & (mm_plot.gust_spd < ef[3]), : ]
	mm_plot_ef3 = mm_plot.loc[(mm_plot.gust_spd >= ef[3]) & (mm_plot.gust_spd < ef[4]), : ]
	mm_plot_ef4 = mm_plot.loc[(mm_plot.gust_spd >= ef[4]) & (mm_plot.gust_spd < ef[5]), : ]
	mm_plot_ef5 = mm_plot.loc[(mm_plot.gust_spd >= ef[5]), : ]

	gust_max = max(mm.gust_spd)

	# get title
	if main_title == "":
		main_title = hur_id + " " + site_name

	# create plot
	plt.figure(figsize=(10,10))

	if plot_type == "time_series":
		plt.xlim(min(mm_plot.dt), max(mm_plot.dt))

	plt.scatter(mm_plot_efx[x_var], mm_plot_efx[y_var], label="No damage", color=ef_col[6])
	if gust_max >= ef[0]:
		plt.scatter(mm_plot_ef0[x_var], mm_plot_ef0[y_var], label="EF0 damage", color=ef_col[0])
	if gust_max >= ef[1]:
		plt.scatter(mm_plot_ef1[x_var], mm_plot_ef1[y_var], label="EF1 damage", color=ef_col[1])
	if gust_max >= ef[2]:
		plt.scatter(mm_plot_ef2[x_var], mm_plot_ef2[y_var], label="EF2 damage", color=ef_col[2])
	if gust_max >= ef[3]:
		plt.scatter(mm_plot_ef3[x_var], mm_plot_ef3[y_var], label="EF3 damage", color=ef_col[3])
	if gust_max >= ef[4]:
		plt.scatter(mm_plot_ef4[x_var], mm_plot_ef4[y_var], label="EF4 damage", color=ef_col[4])
	if gust_max >= ef[5]:
		plt.scatter(mm_plot_ef5[x_var], mm_plot_ef5[y_var], label="EF5 damage", color=ef_col[5])

	plt.title(main_title, fontsize=20)
	plt.xlabel(x_label, fontsize=18)
	plt.ylabel(y_label, fontsize=18)
	plt.legend(loc=legend_loc)
	plt.show()
	plt.clf()

# hurrecon_plot_site_all creates a time-series plot of peak values 
# for all hurricanes for a given site. Optional start and end years
# may be specified. Variables to plot: wind_speed, gust_speed, or
# wind_direction.
#   site_name - name of site
#   start_year - optional start year
#   end_year - optional end year
#   var - variable to plot
#   legend_loc - legend location
#   main_title - optional title
# no return value

def hurrecon_plot_site_all(site_name, start_year='', end_year='', 
	var="wind_speed", legend_loc="upper right", main_title=""):

	# get current working directory
	cwd = os.getcwd()

	# get enhanced Fujita wind speeds & colors
	ef = get_fujita_wind_speeds()
	ef_col = get_fujita_colors()

	# check legend location
	if (not(check_legend_location(legend_loc))):
		legend_loc = "upper right"

	# read data
	peak_file = cwd + "/site-all/" + site_name + " Peak Values.csv"
	check_file_exists(peak_file)
	kk = pd.read_csv(peak_file)
	kk_rows = len(kk)

	# get axis labels
	x_var = "year"
	x_label = "Year"

	if var == "wind_speed":
		y_var = "wind_spd"
		y_label = "Wind Speed (m/s)"
	elif var == "gust_speed":
		y_var = "gust_spd"
		y_label = "Gust Speed (m/s)"
	elif var == "wind_direction":
		y_var = "wind_dir"
		y_label = "Wind Direction (deg)"
	else:
		sys.exit("var must be wind_speed, gust_speed, or wind_direction")

	# subset by year
	yr_list = []

	for i in range(0, kk_rows):
		yr_list.append(int(kk.date_time[i][0:4]))

	kk['year'] = yr_list

	if start_year != "":
		syear = start_year
	else:
		syear = kk.year[0]

	if end_year != "":
		eyear = end_year
	else:
		eyear = kk.year[kk_rows-1]

	kk_plot = kk.loc[(kk.year >= syear) & (kk.year <= eyear), ]

	# subset by enhanced Fujita scale
	kk_plot_efx = kk_plot.loc[(kk_plot.gust_spd < ef[0]), : ]
	kk_plot_ef0 = kk_plot.loc[(kk_plot.gust_spd >= ef[0]) & (kk_plot.gust_spd < ef[1]), : ]
	kk_plot_ef1 = kk_plot.loc[(kk_plot.gust_spd >= ef[1]) & (kk_plot.gust_spd < ef[2]), : ]
	kk_plot_ef2 = kk_plot.loc[(kk_plot.gust_spd >= ef[2]) & (kk_plot.gust_spd < ef[3]), : ]
	kk_plot_ef3 = kk_plot.loc[(kk_plot.gust_spd >= ef[3]) & (kk_plot.gust_spd < ef[4]), : ]
	kk_plot_ef4 = kk_plot.loc[(kk_plot.gust_spd >= ef[4]) & (kk_plot.gust_spd < ef[5]), : ]
	kk_plot_ef5 = kk_plot.loc[(kk_plot.gust_spd >= ef[5]), : ]

	gust_max = max(kk.gust_spd)

	if main_title == "":
		main_title = site_name

	# create plot
	plt.figure(figsize=(10,10))

	plt.xlim(min(kk_plot.year - 10), max(kk_plot.year + 10))

	plt.scatter(kk_plot_efx[x_var], kk_plot_efx[y_var], label="No damage", s=70, color=ef_col[6])
	if gust_max >= ef[0]:
		plt.scatter(kk_plot_ef0[x_var], kk_plot_ef0[y_var], label="EF0 damage", s=70, color=ef_col[0])
	if gust_max >= ef[1]:
		plt.scatter(kk_plot_ef1[x_var], kk_plot_ef1[y_var], label="EF1 damage", s=70, color=ef_col[1])
	if gust_max >= ef[2]:
		plt.scatter(kk_plot_ef2[x_var], kk_plot_ef2[y_var], label="EF2 damage", s=70, color=ef_col[2])
	if gust_max >= ef[3]:
		plt.scatter(kk_plot_ef3[x_var], kk_plot_ef3[y_var], label="EF3 damage", s=70, color=ef_col[3])
	if gust_max >= ef[4]:
		plt.scatter(kk_plot_ef4[x_var], kk_plot_ef4[y_var], label="EF4 damage", s=70, color=ef_col[4])
	if gust_max >= ef[5]:
		plt.scatter(kk_plot_ef5[x_var], kk_plot_ef5[y_var], label="EF5 damage", s=70, color=ef_col[5])

	plt.title(main_title, fontsize=20)
	plt.xlabel(x_label, fontsize=18)
	plt.ylabel(y_label, fontsize=18)
	plt.legend(loc=legend_loc)
	plt.show()
	plt.clf()

# hurrecon_plot_tracks creates a regional plot of the land-water file
# and selected hurricane tracks.
#   select - show all positions (all), only positions used as
#     model input (model), or none (none)
#   wind_min - the minimum value of maximum sustained wind speed 
#    (meters/second)
#   main_title - optional title
#   colormap - color palette
# no return value

def hurrecon_plot_tracks(select="all", wind_min=33, main_title="", 
	colormap="default"):
	
	# get current working directory
	cwd = os.getcwd()

	# read land-water file
	land_water_file = cwd + "/input/land_water.tif"
	check_file_exists(land_water_file)
	land_water = rio.open(land_water_file)

	lat_min = land_water.bounds.bottom
	lat_max = land_water.bounds.top
	lon_min = land_water.bounds.left
	lon_max = land_water.bounds.right

	# read boundaries file
	boundaries_file = cwd + "/vector/boundaries.shp"
	shapefile = fiona.open(boundaries_file, "r")
	features = [feature["geometry"] for feature in shapefile]
	shapefile.close()

	# get hurricane tracks
	ids_file = cwd + "/input/ids.csv"
	check_file_exists(ids_file)
	ii = pd.read_csv(ids_file)

	track_file = cwd + "/input/tracks.csv"
	check_file_exists(track_file)
	tt = pd.read_csv(track_file)

	track_all_file = cwd + "/input/tracks_all.csv"
	check_file_exists(track_all_file)
	tt_all = pd.read_csv(track_all_file)
	
	# color palettes
	if colormap == "default":
		cmap = plt.get_cmap('Greens')
	
	else:
		cmap = plt.get_cmap(colormap)

	# get title
	if main_title == "":
		main_title = "Hurricane Tracks (" + str(wind_min) + " m/s)"

	# create plot
	vals = [0, 1, 2]
	labs = ["", "water", "land"]

	fig, ax = plt.subplots(figsize=(15, 15))
	plt.axis([lon_min, lon_max, lat_min, lat_max])
	plt.xlabel('Longitude (degrees)')
	plt.ylabel('Latitude (degrees)')
	patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
	ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
	plt.title(main_title)
	img = land_water.read(1)		
	plt.imshow(img, cmap=cmap, vmin=0.9)
	cbar = plt.colorbar(cmap=cmap, shrink=0.3)
	cbar.set_ticks(vals)
	cbar.set_ticklabels(labs)
	cbar.ax.set_title("   cover")

	if select == "all":
		for i in range(0, len(ii)):
			if ii.wind_peak[i] >= wind_min:
				hur_id = ii.loc[i, "hur_id"]
				xx = tt_all.loc[tt_all.hur_id == hur_id, ]
				x_coord = list(xx.longitude)
				y_coord = list(xx.latitude)
				plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)

	elif select == "model":
		for i in range(0, len(ii)):
			if ii.wind_peak[i] >= wind_min:
				hur_id = ii.loc[i, "hur_id"]
				xx = tt.loc[tt.hur_id == hur_id, ]
				x_coord = list(xx.longitude)
				y_coord = list(xx.latitude)
				plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)

	show((land_water, 1), cmap=cmap, vmin=0.9)
	plt.clf()

	# close land-water file
	land_water.close()

# hurrecon_plot_region creates regional plots of peak wind speed, peak 
# enhanced Fujita scale, peak wind direction, peak cardinal wind direction,
# and duration of EF0, EF1, EF2, EF3, EF4, and EF5 winds for a given hurricane.
# Variables to plot: wind_speed, fujita_scale, wind_direction, wind_compass, 
# ef0_duration, ef1_duration, ef2_duration, ef3_duration, ef4_duration,
# and ef5_duration.
#   hur_id - hurricane id
#   var - variable to plot
#   positions - whether to plot original positions
#   main_title - optional title
#   colormap - color palette
# no return value

def hurrecon_plot_region(hur_id, var="fujita_scale", positions=False, 
	main_title="", colormap="default"):
	
	# get current working directory
	cwd = os.getcwd()

	# read GeoTiff file
	hur_tif_file = cwd + "/region/" + hur_id + ".tif"
	check_file_exists(hur_tif_file)
	hur_tif = rio.open(hur_tif_file)

	# get extent
	lat_min = hur_tif.bounds.bottom
	lat_max = hur_tif.bounds.top
	lon_min = hur_tif.bounds.left
	lon_max = hur_tif.bounds.right

	# read boundaries file
	boundaries_file = cwd + "/vector/boundaries.shp"
	shapefile = fiona.open(boundaries_file, "r")
	features = [feature["geometry"] for feature in shapefile]
	shapefile.close()

	# get hurricane track
	track_all_file = cwd + "/input/tracks_all.csv"
	check_file_exists(track_all_file)
	tt_all = pd.read_csv(track_all_file)
	index = np.where(tt_all.hur_id == hur_id)[0].tolist()

	lat_list = [tt_all.latitude[index[i]] for i in range(0, len(index))]
	lon_list = [tt_all.longitude[index[i]] for i in range(0, len(index))]

	# get colormaps with white background
	viridis = plt.get_cmap('viridis')
	viridis.set_under('white') 	

	rainbow = plt.get_cmap('rainbow')
	rainbow.set_under('white')

	ef_col = get_fujita_colors()
  
	ff_all_vals = [0, 1, 2, 3, 4, 5, 6, 7]
	ff_all_labs = ['', 'None', 'EF0', 'EF1', 'EF2', 'EF3', 'EF4', 'EF5']
	ff_all_cols = [ef_col[6], ef_col[0], ef_col[1], ef_col[2], ef_col[3], ef_col[4], ef_col[5]]

	ff_max = hur_tif.read(2).max()
	ff_cols = []

	for i in range(0, ff_max):
		ff_cols.append(ff_all_cols[i])

	fscale = matplotlib.colors.ListedColormap(ff_cols)
	fscale.set_under('white')

	# color palettes
	if colormap == "default":
		if var == "fujita_scale":
			cmap = fscale

		elif var == "wind_compass":
			cmap = rainbow

		else:
			cmap = viridis
	else:
		cmap = plt.get_cmap(colormap)
		cmap.set_under('white')

	# create plot
	if var == "wind_speed":
		if np.amax(hur_tif.read(1)) > 0:
			if main_title == "":
				main_title = hur_id + " Peak Wind Speed"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(1)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   m/sec")
			show((hur_tif, 1), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No wind speed\n")


	elif var == "fujita_scale":
		if np.amax(hur_tif.read(2)) > 0:
			if main_title == "":
				main_title = hur_id + " Peak Fujita Scale"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(2)
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.set_ticks(ff_all_vals)
			cbar.set_ticklabels(ff_all_labs)
			cbar.ax.set_title("   EF Scale")
			show((hur_tif, 2), cmap=cmap, vmin=0.9)
			plt.clf()
		else:
			print("No Fujita values\n")


	elif var == "wind_direction":
		if np.amax(hur_tif.read(3)) > 0:
			if main_title == "":
				main_title = hur_id + " Peak Wind Direction"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(3)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   degrees")
			show((hur_tif, 3), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No wind direction\n")

	elif var == "wind_compass":
		if np.amax(hur_tif.read(4)) > 0:
			if main_title == "":
				main_title = hur_id + " Peak Wind Direction"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(4)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
			cbar.set_ticklabels(['', 'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
			cbar.ax.set_title("   direction")
			show((hur_tif, 4), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No wind compass\n")

	elif var == "ef0_duration":
		if np.amax(hur_tif.read(5)) > 0:
			if main_title == "":
				main_title = hur_id + " EF0 Winds"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(5)/60
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   hours")
			show((hur_tif, 5), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No EF0 winds\n")

	elif var == "ef1_duration":
		if np.amax(hur_tif.read(6)) > 0:
			if main_title == "":
				main_title = hur_id + " EF1 Winds"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(6)/60
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   hours")
			show((hur_tif, 6), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No EF1 winds\n")

	elif var == "ef2_duration":
		if np.amax(hur_tif.read(7)) > 0:
			if main_title == "":
				main_title = hur_id + " EF2 Winds"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(7)/60
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   hours")
			show((hur_tif, 7), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No EF2 winds\n")

	elif var == "ef3_duration":
		if np.amax(hur_tif.read(8)) > 0:
			if main_title == "":
				main_title = hur_id + " EF3 Winds"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(8)/60
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   hours")
			show((hur_tif, 8), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No EF3 winds\n")

	elif var == "ef4_duration":
		if np.amax(hur_tif.read(9)) > 0:
			if main_title == "":
				main_title = hur_id + " EF4 Winds"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(9)/60
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   hours")
			show((hur_tif, 9), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No EF4 winds\n")
	
	elif var == "ef5_duration":
		if np.amax(hur_tif.read(10)) > 0:
			if main_title == "":
				main_title = hur_id + " EF5 Winds"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.title(main_title)
			img = hur_tif.read(10)/60
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   hours")
			show((hur_tif, 10), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No EF5 winds\n")

	else:
		sys.exit("var must be wind_speed, fujita_scale, wind_direction, wind_compass, ef0_duration, ef1_duration, ef2_duration, ef3_duration, ef4_duration, or ef5_duration")

	# close GeoTiff file
	hur_tif.close()

# hurrecon_plot_region_dt creates regional plots of enhanced Fujita scale, 
# wind speed, wind direction, and cardinal wind direction for a given hurricane
# at a specified datetime. Variables to plot: wind_speed, fujita_scale, 
# wind_direction, or wind_compass.
#   hur_id - hurricane id
#   dt - datetime in the format YYYY-MM-DDThh:mm
#   var - variable to plot
#   positions - whether to plot original positions
#   main_title - optional title
#   colormap - color palette
# no return value

def hurrecon_plot_region_dt(hur_id, dt, var="fujita_scale", positions=False,
	main_title="", colormap="default"):

	# get current working directory
	cwd = os.getcwd()

	# read GeoTiff file
	dt2 = dt.replace(":", "")
	hur_tif_file = cwd + "/region-dt/" + hur_id + " " + dt2 + ".tif"
	check_file_exists(hur_tif_file)
	hur_tif = rio.open(hur_tif_file)

	# get extent
	lat_min = hur_tif.bounds.bottom
	lat_max = hur_tif.bounds.top
	lon_min = hur_tif.bounds.left
	lon_max = hur_tif.bounds.right

	# read boundaries file
	boundaries_file = cwd + "/vector/boundaries.shp"
	shapefile = fiona.open(boundaries_file, "r")
	features = [feature["geometry"] for feature in shapefile]
	shapefile.close()

	# get current location
	tt = read_hurricane_track_file(hur_id)
	pp = get_values_at_datetime(hur_id, tt, dt)

	# get hurricane track
	track_all_file = cwd + "/input/tracks_all.csv"
	check_file_exists(track_all_file)
	tt_all = pd.read_csv(track_all_file)
	index = np.where(tt_all.hur_id == hur_id)[0].tolist()

	lat_list = [tt_all.latitude[index[i]] for i in range(0, len(index))]
	lon_list = [tt_all.longitude[index[i]] for i in range(0, len(index))]

	# get colormaps with white background
	viridis = plt.get_cmap('viridis')
	viridis.set_under('white') 	

	rainbow = plt.get_cmap('rainbow')
	rainbow.set_under('white') 	

	ef_col = get_fujita_colors()
  
	ff_all_vals = [0, 1, 2, 3, 4, 5, 6, 7]
	ff_all_labs = ['', 'None', 'EF0', 'EF1', 'EF2', 'EF3', 'EF4', 'EF5']
	ff_all_cols = [ef_col[6], ef_col[0], ef_col[1], ef_col[2], ef_col[3], ef_col[4], ef_col[5]]

	ff_max = hur_tif.read(2).max()
	ff_cols = []

	for i in range(0, ff_max):
		ff_cols.append(ff_all_cols[i])

	fscale = matplotlib.colors.ListedColormap(ff_cols)
	fscale.set_under('white')

	# color palettes
	if colormap == "default":
		if var == "fujita_scale":
			cmap = fscale

		elif var == "wind_compass":
			cmap = rainbow

		else:
			cmap = viridis
	else:
		cmap = plt.get_cmap(colormap)
		cmap.set_under('white')

	# create plot
	if var == "fujita_scale":
		if np.amax(hur_tif.read(2)) > 0:
			if main_title == "":
				main_title = hur_id + ' Fujita Scale ' + dt
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.plot(pp.lon[0], pp.lat[0], 'o', color='brown', markersize=6)
			plt.title(main_title)
			img = hur_tif.read(2)
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.set_ticks(ff_all_vals)
			cbar.set_ticklabels(ff_all_labs)
			cbar.ax.set_title("   EF Scale")
			show((hur_tif, 2), cmap=cmap, vmin=0.9)
			plt.clf()
		else:
			print("No Fujita values\n")

	elif var == "wind_speed":
		if np.amax(hur_tif.read(1)) > 0:
			if main_title == "":
				main_title = hur_id + ' Wind Speed ' + dt
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.plot(pp.lon[0], pp.lat[0], 'o', color='brown', markersize=6)
			plt.title(main_title)
			img = hur_tif.read(1)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   m/sec")
			show((hur_tif, 1), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No wind speed\n")

	elif var == "wind_direction":
		if np.amax(hur_tif.read(3)) > 0:
			if main_title == "":
				main_title = hur_id + ' Wind Direction ' + dt
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.plot(pp.lon[0], pp.lat[0], 'o', color='brown', markersize=6)
			plt.title(main_title)
			img = hur_tif.read(3)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar =  plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   degrees")
			show((hur_tif, 3), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No wind direction\n")

	elif var == "wind_compass":
		if np.amax(hur_tif.read(4)) > 0:
			if main_title == "":
				main_title = hur_id + ' Wind Direction ' + dt
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='brown', linewidth=0.8)
			if (positions == True):
				plt.plot(lon_list, lat_list, 'o', color='brown', markersize=3)
			plt.plot(pp.lon[0], pp.lat[0], 'o', color='brown', markersize=6)
			plt.title(main_title)
			img = hur_tif.read(4)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
			cbar.set_ticklabels(['', 'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
			cbar.ax.set_title("   direction")
			show((hur_tif, 4), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No wind compass\n")

	else:
		sys.exit("var must be wind_speed, fujita_scale, wind_direction, wind_compass")

	# close GeoTiff file
	hur_tif.close()

# hurrecon_plot_region_all creates regional plots of maximum enhanced Fujita
# value and number of storms for each enhanced Fujita value for all hurricanes.
# Variables to plot: efmax, ef0, ef1, ef2, ef3, ef4, or ef5.
#   var - variable to plot
#   tracks - whether to also plot hurricane tracks
#   main_title - optional title
#   colormap - color palette
# no return value

def hurrecon_plot_region_all(var="efmax", tracks=False, main_title="",
	colormap="default"):
	
	# get current working directory
	cwd = os.getcwd()

	# read GeoTiff file
	sum_tif_file = cwd + "/region-all/" + "summary.tif"
	check_file_exists(sum_tif_file)
	sum_tif = rio.open(sum_tif_file)

	# get extent
	lat_min = sum_tif.bounds.bottom
	lat_max = sum_tif.bounds.top
	lon_min = sum_tif.bounds.left
	lon_max = sum_tif.bounds.right

	# read boundaries file
	boundaries_file = cwd + "/vector/boundaries.shp"
	shapefile = fiona.open(boundaries_file, "r")
	features = [feature["geometry"] for feature in shapefile]
	shapefile.close()

	# get colormaps with white background
	viridis = plt.get_cmap('viridis')
	viridis.set_under('white') 	

	ef_col = get_fujita_colors()
  
	ff_all_vals = [0, 1, 2, 3, 4, 5, 6, 7]
	ff_all_labs = ['', 'None', 'EF0', 'EF1', 'EF2', 'EF3', 'EF4', 'EF5']
	ff_all_cols = [ef_col[6], ef_col[0], ef_col[1], ef_col[2], ef_col[3], ef_col[4], ef_col[5]]

	ff_max = sum_tif.read(1).max()
	ff_cols = []

	for i in range(0, ff_max):
		ff_cols.append(ff_all_cols[i])

	fscale = matplotlib.colors.ListedColormap(ff_cols)
	fscale.set_under('white')

	# color palettes
	if colormap == "default":
		if var == "efm":
			cmap = fscale

		else:
			cmap = viridis
	else:
		cmap = plt.get_cmap(colormap)
		cmap.set_under('white')

	# get hurricane tracks
	if tracks == True:
		ids_file = cwd + "/input/ids.csv"
		check_file_exists(ids_file)
		ii = pd.read_csv(ids_file)

		track_all_file = cwd + "/input/tracks_all.csv"
		check_file_exists(track_all_file)
		tt = pd.read_csv(track_all_file)

		summary_file = cwd + "/region-all/summary.csv"
		check_file_exists(summary_file)
		kk = pd.read_csv(summary_file)

	# create plot
	if var == "efmax":
		if np.amax(sum_tif.read(1)) > 0:
			if main_title == "":
				main_title = "Peak Fujita Scale"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title(main_title)
			img = sum_tif.read(1)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.set_ticks(ff_all_vals)
			cbar.set_ticklabels(ff_all_labs)
			cbar.ax.set_title("   EF Scale")
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 0
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						x_coord = list(xx.longitude)
						y_coord = list(xx.latitude)
						plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)
			show((sum_tif, 1), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No Fujita values\n")

	elif var == "ef0":
		if np.amax(sum_tif.read(2)) > 0:
			if main_title == "":
				main_title = "Fujita Scale 0"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title(main_title)
			img = sum_tif.read(2)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   storms")
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 0
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						x_coord = list(xx.longitude)
						y_coord = list(xx.latitude)
						plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)
			show((sum_tif, 2), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No F0 values\n")

	elif var == "ef1":
		if np.amax(sum_tif.read(3)) > 0:
			if main_title == "":
				main_title = "Fujita Scale 1"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title(main_title)
			img = sum_tif.read(3)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   storms")
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 1
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						x_coord = list(xx.longitude)
						y_coord = list(xx.latitude)
						plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)
			show((sum_tif, 3), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No F1 values\n")

	elif var == "ef2":
		if np.amax(sum_tif.read(4)) > 0:
			if main_title == "":
				main_title = "Fujita Scale 2"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title(main_title)
			img = sum_tif.read(4)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   storms")
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 2
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						x_coord = list(xx.longitude)
						y_coord = list(xx.latitude)
						plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)
			show((sum_tif, 4), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No F2 values\n")

	elif var == "ef3":
		if np.amax(sum_tif.read(5)) > 0:
			if main_title == "":
				main_title = "Fujita Scale 3"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title(main_title)
			img = sum_tif.read(5)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   storms")
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 3
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						x_coord = list(xx.longitude)
						y_coord = list(xx.latitude)
						plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)
			show((sum_tif, 5), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No F3 values\n")
	
	elif var == "ef4":
		if np.amax(sum_tif.read(6)) > 0:
			if main_title == "":
				main_title = "Fujita Scale 4"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title(main_title)
			img = sum_tif.read(6)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   storms")
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 4
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						x_coord = list(xx.longitude)
						y_coord = list(xx.latitude)
						plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)
			show((sum_tif, 6), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No F4 values\n")

	elif var == "ef5":
		if np.amax(sum_tif.read(7)) > 0:
			if main_title == "":
				main_title = "Fujita Scale 5"
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title(main_title)
			img = sum_tif.read(7)		
			plt.imshow(img, cmap=cmap, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.ax.set_title("   storms")
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 5
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						x_coord = list(xx.longitude)
						y_coord = list(xx.latitude)
						plt.plot(x_coord, y_coord, color='grey', linewidth=0.8)
			show((sum_tif, 7), cmap=cmap, vmin=0.9)
			plt.clf()	
		else:
			print("No F5 values\n")

	else:
		sys.exit("var must be efmax, ef0, ef1, ef2, ef3, ef4, or ef5")

	# close GeoTiff file
	sum_tif.close()
