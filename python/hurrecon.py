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
# May 2020

# Python version 3.7.4.

# Note: Pandas datetime functions are currently limited to the years 
# 1678-2262 and so are not used here.


### MODULES ###############################################

import os
import sys
import math
import time
import fiona
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.plot import show
from rasterio.mask import mask
from descartes import PolygonPatch


### INTERNAL FUNCTIONS ####################################

# get_fujita_wind_speeds returns a list containing the minimum 3-second 
# gust speed (meters/second) for each enhanced Fujita class.
#   returns a list of 3-second gust speeds

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
#   returns a vector of colors

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
#   returns a time difference formatted as hh:mm:ss.

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
#   returns a time difference in milliseconds

def format_time_difference_ms(start_time, end_time):
	t_diff = str(round(1000*(end_time - start_time))).zfill(3)

	return t_diff

# check_file_exists displays an error message and stops execution if
# the specified file does not exist.
#   file_name - name of file
#   no return value

def check_file_exists(file_name):
	if os.path.exists(file_name) == False or os.path.isfile(file_name) == False:
		print("\nFile not found:", file_name, "\n")
		sys.exit()

# read_site_file reads a site file and returns a list containing the
# latitude (degrees), longitude (degrees), and cover type (water=1, land=2) 
# for the specified site.
#   site_name - name of site
#   returns a list of latitude, longitude, and cover type

def read_site_file(site_name):
	cwd = os.getcwd()
	site_file = cwd + "/input/sites.csv"
	check_file_exists(site_file)
	ss = pd.read_csv(site_file)

	# get site location & cover type
	index = np.where(ss.site_name == site_name)[0].tolist()
   
	if len(index) == 0:
		print("\nSite not found\n")
		sys.exit()

	site_latitude = ss.latitude[index[0]]
	site_longitude = ss.longitude[index[0]]
	cover_type = ss.cover_type[index[0]]

	return [site_latitude, site_longitude, cover_type]

# read_parameter_file reads a parameter file and returns a list containing
# the radius of maximum wind (rmw) (kilometers) and profile exponent (s_par). 
# If width is True, parameters are returned for the specified hurricane, 
# if available; otherwise parameters for ALL are returned.
#   hur_id - hurricane id
#   width - whether to use width parameters for the specified hurricane
#   returns a list of rmw and s_par

def read_parameter_file(hur_id, width):
	cwd = os.getcwd()
	par_file = cwd + "/input/parameters.csv"
	check_file_exists(par_file)
	pp = pd.read_csv(par_file)

	# get rmw & s_par parameters
	if width:
		index = np.where(pp.hur_id == hur_id)[0].tolist()
	else:
		index = np.where(pp.hur_id == "ALL")[0].tolist()

	if len(index) == 0:
		print("\nParameter file must contain an entry for ALL\n")
		sys.exit()

	rmw = pp.rmw[index[0]]
	s_par = pp.s_par[index[0]]

	return [rmw, s_par]

# get_fixed_model_parameters returns a list of fixed model parameters,
# including asymmetry factor, inflow angle, friction factor, and gust factor.
#   cover_type - cover type (1=water, 2=land)
#   returns a list of fixed model parameters

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
		print("\nCover type must be 1 (water) or 2 (land)\n")
		sys.exit()

	return [asymmetry_factor, inflow_angle, friction_factor, gust_factor]

# get_time_step calculates the time step (minutes) for regional modeling, 
# assuming a maximum hurricane forward speed of 20 meters per second (1200 
# meters per minute). Values are rounded to the nearest 1, 2, 3, 5, 10, 15, 
# 30, or 60 minutes.
#   returns a time step

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
#   returns the Julian day and fraction

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
#   returns a standard datetime

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
#   returns a data frame of track data

def read_hurricane_track_file(hur_id):
	# read hurricane track file
	cwd = os.getcwd()
	track_file = cwd + "/input/tracks.csv"
	check_file_exists(track_file)
	zz = pd.read_csv(track_file)

	# subset by hurricane name
	xx = zz.loc[zz.hur_id == hur_id]

	if len(xx) == 0:
		print("\nHurricane not in track file\n")
		sys.exit()

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
#   returns a data frame of interpolated data

def interpolate_hurricane_location_max_wind(tt, time_step):
	tt_rows = len(tt)

	# initialize lists
	jd_list = []
	dt_list = []
	yr_list = []
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

		dt = [''] * (new_rows - 1)
		yr = [int(tt.date_time[0][0:4])] * (new_rows - 1)

		jd_list.extend(jd)
		dt_list.extend(dt)
		yr_list.extend(yr)
		lat_list.extend(lat)
		lon_list.extend(lon)
		wmax_list.extend(wmax)

	# add final row
	jd_list.append(tt.jd[tt_rows-1])
	dt_list.append('')
	yr_list.append(int(tt.date_time[tt_rows-1][0:4]))
	lat_list.append(tt.latitude[tt_rows-1])
	lon_list.append(tt.longitude[tt_rows-1])
	wmax_list.append(tt.wind_max[tt_rows-1])

	# create data frame for modeled data
	mm_cols = ['date_time', 'year', 'jd', 'latitude', 'longitude', 'wind_max']
	
	mm = pd.DataFrame(data=list(zip(dt_list, yr_list, jd_list, lat_list, lon_list, wmax_list)), 
		columns=mm_cols)

	return mm

# calculate_range_bearing returns a list containing the range (kilometers)
# and bearing (degrees) from one point to another using the latitude & longitude
# of each point.
#   lat1 - latitude of first point
#   lon1 - longitude of first point
#   lat2 - latitude of second point
#   lon2 - longitude of second point
#   returns a list containing range & bearing

def calculate_range_bearing(lat1, lon1, lat2, lon2):
	R = 6367 # radius of earth in kilometers (at latitude 45 degrees)

	d2r = 0.017453292519943295  # pi / 180
	r2d = 57.29577951308232  # 180 / pi

	# nearly same point
	if abs(lat2 - lat1) < 0.000001 and abs(lon2 - lon1) < 0.000001:
		rang = 0
		bear = 0

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

		cos_rlat1 = math.cos(rlat1)
		cos_rlat2 = math.cos(rlat2)

		A = (math.sin((rlat2-rlat1)/2))**2 + cos_rlat1*cos_rlat2*(math.sin((rlon2-rlon1)/2))**2
		C = 2 * math.atan2(math.sqrt(A), math.sqrt(1-A))
		rang = R * C

		# same longitude
		if lon1 == lon2:
			if lat1 > lat2:
				bear = 180
			else:
				bear = 0

		# different longitude
		else:

			B2 = math.atan2(math.sin(rlon2-rlon1)*cos_rlat2, cos_rlat1*math.sin(rlat2) - 
				math.sin(rlat1)*cos_rlat2*math.cos(rlon2-rlon1))

			# convert radians to degrees
			B = r2d*B2

			if lon1 < lon2:
				# quadrants I, IV
				bear = B
			else:  
				# quadrants II, III
				bear = 360 + B

	return [rang, bear]

# interpolate_hurricane_speed_bearing performs a linear interpolation of hurricane
# speed (meters/second) and bearing (degrees) along a hurricane track based on
# mid-segment values.
#   tt - data frame of track values
#   mm - data frame of modeled values
#   returns a data frame of modeled values

def interpolate_hurricane_speed_bearing(tt, mm):
	tt_rows = len(tt)

	# intialize lists
	jd_list = []
	spd_list = []
	bear_list = []

	# calculate mid-segment hurricane speed & bearing
	for i in range(0, tt_rows-1):
		hur_range_bear = calculate_range_bearing(tt.latitude[i], tt.longitude[i], 
			tt.latitude[i+1], tt.longitude[i+1])

		interval_sec = (tt.jd[i+1] - tt.jd[i]) * 1440 * 60

		jd_list.append(tt.jd[i] + (tt.jd[i+1] - tt.jd[i])/2)
		spd_list.append(1000*hur_range_bear[0]/interval_sec)
		bear_list.append(hur_range_bear[1])

	# create data frame for mid-segment values
	vv_cols = ['jd', 'hur_bear', 'hur_spd']
	vv = pd.DataFrame(data=list(zip(jd_list, bear_list, spd_list)), columns=vv_cols)

	vv_rows = len(vv)

	# intialize lists
	bear_list = []
	spd_list = []

 	# interpolate hurricane speed & bearing for each segment
	for i in range(0, vv_rows+1):
		# before mid-point of 1st segment
		if i == 0:
			index = np.where(mm.jd <= vv.jd[0])[0].tolist()
			new_rows = len(index)

			bear = [vv.hur_bear[0]] * new_rows
			spd = [vv.hur_spd[0]] * new_rows

			bear_list.extend(bear)
			spd_list.extend(spd)

		# interpolate between mid-points
		elif i <= vv_rows-1:
			index = np.where((mm.jd > vv.jd[i-1]) & (mm.jd <= vv.jd[i]))[0].tolist()
			new_rows = len(index)

			# bearing
			b1 = vv.hur_bear[i-1]
			b2 = vv.hur_bear[i]

			if b2 - b1 > 180:
				b1 = b1 + 360
			elif b1 - b2 > 180:
				b2 = b2 + 360

			bear = np.linspace(start=b1, stop=b2, num=new_rows).tolist()
			
			# speed
			spd = np.linspace(start=vv.hur_spd[i-1], stop=vv.hur_spd[i], num=new_rows).tolist()

			bear_list.extend(bear)
			spd_list.extend(spd)

		# after mid-point of last segment
		else:
			index = np.where(mm.jd > vv.jd[vv_rows-1])[0].tolist()
			new_rows = len(index)

			bear = [vv.hur_bear[vv_rows-1]] * new_rows
			spd = [vv.hur_spd[vv_rows-1]] * new_rows

			bear_list.extend(bear)
			spd_list.extend(spd)
			
	# adjust bearing as needed
	for i in range(0, len(bear_list)):
		if bear_list[i] < 0:
			bear_list[i] = bear_list[i] + 360

		if bear_list[i] > 360:
			bear_list[i] = bear_list[i] - 360

	# add to modeled data frame
	mm['hur_bear'] = bear_list
	mm['hur_spd'] = spd_list

	return mm

# calculate_site_range_bearing calculates the range (kilometers) and bearing
# (degrees) from a site to the hurricane center.
#   mm - data frame of modeled values
#   site_latitude - latitude of site
#   site_longitude - longitude of site
#   returns a data frame of modeled values

def calculate_site_range_bearing(mm, site_latitude, site_longitude):
	mm_rows = len(mm)

	site_range_bear = [calculate_range_bearing(site_latitude, site_longitude, 
		mm.latitude[i], mm.longitude[i]) for i in range(0, mm_rows)]

	mm['site_bear']  = [site_range_bear[i][1] for i in range(0, mm_rows)]
	mm['site_range'] = [site_range_bear[i][0] for i in range(0, mm_rows)]

	return mm

# calculate_wind_direction calculates the wind direction (degrees) at the
# specified site.
#   hurr_lat - latitude of hurricane (degrees)
#   site_bear - bearing from site to hurricane center (degrees)
#   inflow_angle - cross-isobar inflow angle (degrees)
#   returns a calculated wind direction

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
#   s_par - profile exponent
#   asymmetry_factor - asymmetry factor
#   friction_factor - friction factor
#   returns a calculated sustained wind speed (meters/second)

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
#   wind_spd - sustained wind speed
#   gust_factor - gust factor
#   returns wind gust speed

def calculate_wind_gust(wind_spd, gust_factor):
	gust_spd = gust_factor * wind_spd

	return gust_spd

# calculate_enhanced_fujita_scale returns the enhanced Fujita scale value
# based on the wind gust speed (meters/second)
#   gust_spd - wind gust speed
#   returns the enhanced Fujita scale value

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
#   mm - data frame of modeled values
#   inflow_angle - cross-isobar inflow angle (degrees)
#   cover_type - cover type (1=water, 2=land)
#   rmw - radius of maximum winds (kilometers)
#   s_par - profile exponent
#   asymmetry_factor - asymmetry factor
#   friction_factor - friction factor
#   gust_factor - gust factor
#   returns a data frame of modeled values

def calculate_wind_speed_direction(mm, inflow_angle, cover_type, rmw, 
	s_par, asymmetry_factor, friction_factor, gust_factor):
	
	mm_rows = len(mm)

	# wind speed
	wspd = [calculate_wind_speed(mm.site_bear[i], mm.site_range[i], mm.latitude[i], 
			mm.hur_bear[i], mm.hur_spd[i], mm.wind_max[i], rmw, s_par, asymmetry_factor, 
			friction_factor) for i in range(0, mm_rows)]

	# gust speed
	gust = [calculate_wind_gust(wspd[i], gust_factor) for i in range(0, mm_rows)]

	# wind direction
	wdir = [calculate_wind_direction(mm.latitude[i], mm.site_bear[i], inflow_angle)
		for i in range(0, mm_rows)]

	# enhanced Fujita scale
	ef = [calculate_enhanced_fujita_scale(gust[i]) for i in range(0, mm_rows)]

	mm['rmw'] = rmw
	mm['s_par'] = s_par
	mm['wind_dir'] = wdir
	mm['wind_spd'] = wspd
	mm['gust_spd'] = gust
	mm['ef_sca'] = ef

	return mm

# add_standard_date_time adds a standard datetime column in the format
# YYYY-MM-DDThh:mm to a data frame of modeled values.
#   mm - data frame of modeled values
#   returns a data frame of modeled values

def add_standard_date_time(mm):
	mm_rows = len(mm)

	mm['date_time'] = [calculate_standard_datetime(mm.year[i], mm.jd[i])
		for i in range(0, mm_rows)]
	return mm

# get_peak_values returns a data frame of peak values for a given
# hurricane and site.
#   hur_id - hurricane id
#   site_name - name of site
#   mm - data frame of modeled values
#   returns a data frame of peak values

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
		'gust_spd', "ef_sca", "ef0", "ef1", "ef2", "ef3", "ef4", "ef5"]

	kk = pd.DataFrame(columns=kk_cols)

	kk.loc[0] = [site_name, hur_id, date_time, wind_dir, wind_spd, gust_spd, ef_sca, 
		ef0, ef1, ef2, ef3, ef4,ef5]

	return kk

# get_regional_peak_wind calculates the peak wind speed (meters/second), enhanced 
# Fujita scale, wind direction (degrees), cardinal wind direction, gale wind duration 
# (minutes), and hurricane wind duration (minutes) for a given hurricane over a region.
#   hur_id - hurricane id
#   mm - data frame of modeled values
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   water - whether to calculate values over water
#   timing - whether to report progress
#   returns a list containng 6 raster layers

def get_regional_peak_wind(hur_id, mm, width, time_step, water, timing):
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
	gg = np.zeros((nrows, ncols), dtype=np.int16)  # duration of gale force winds (minutes)
	hh = np.zeros((nrows, ncols), dtype=np.int16)  # duration of hurricane winds (minutes)

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

	# record total elasped time if timing is True
	if timing == True:
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

				# range & bearing from site to hurricane center
				site_range_bear = [calculate_range_bearing(site_latitude, site_longitude,
						mm.latitude[k], mm.longitude[k]) for k in range(0, len(mm))]

				# wind speed (m/s)
				wspd = [calculate_wind_speed(site_range_bear[k][1], site_range_bear[k][0], mm.latitude[k], 
					mm.hur_bear[k], mm.hur_spd[k], mm.wind_max[k], rmw, s_par, asymmetry_factor, friction_factor)
					for k in range(0, len(mm))]

				# update values if gale or higher
				wmax = max(wspd)
				if wmax >= 17.5:
					# duration of gale force winds (minutes)
					index = [index for index in wspd if index >= 17.5]
					gg[nrows-i-1][j] = gg[nrows-i-1][j] + len(index)*time_step

					# duration of hurricane force winds (minutes)
					index = [index for index in wspd if index >= 33]
					hh[nrows-i-1][j] = hh[nrows-i-1][j] + len(index)*time_step

					# peak wind speed
					index = wspd.index(wmax)
					ss[nrows-i-1][j] = int(round(wmax))

					# peak wind direction			
					wdir = calculate_wind_direction(mm.latitude[0], site_range_bear[index][1], inflow_angle)
					dd[nrows-i-1][j] = int(round(wdir))

		# report progress
		if timing == True:
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
	if timing == True:
		elapsed_time = format_time_difference_hms(start_time, time.time())
		print("\r", elapsed_time, "\n", end="")

	return [ss, ff, dd, cc, gg, hh]

# get_regional_summary compiles regional results for all hurricanes.
# Results are saved in a GeoTiff file (summary.tif) with 7 layers and in
# a CSV file of hurricane ids and maximum enhanced Fujita scale values
# (summary.csv) on the region subdirectory.
#   no return value

def get_regional_summary():
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

	# create lists for peak Fujita value across region
	hur_id_list = []
	efmax_list = []

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
		hur_tif_file = cwd + "/region/" + hur_id + ".tif"
		check_file_exists(hur_tif_file)
		hur_tif = rio.open(hur_tif_file)

		# get enhanced Fujita scale layer
		ff_array = hur_tif.read(2) # enhanced Fujita scale

		# update peak Fujita value
		efmax = np.amax(ff_array) - 2

		hur_id_list.append(hur_id)
		efmax_list.append(efmax)

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

	# save to csv file
	peak_cols = ['hur_id', 'efmax']
	peak = pd.DataFrame(data=list(zip(hur_id_list, efmax_list)), columns=peak_cols)

	peak_file = cwd + "/region/summary.csv"
	peak.to_csv(peak_file, index=False)

	# save to GeoTiff file
	sum_tif_file = cwd + "/region/summary.tif"

	profile = land_water.profile
	profile.update(dtype='int16', nodata=-9999, count=7)

	sum_tif = rio.open(sum_tif_file, 'w', **profile)

	sum_tif.write(efm, 1)
	sum_tif.write(ef0, 2)
	sum_tif.write(ef1, 3)
	sum_tif.write(ef2, 4)
	sum_tif.write(ef3, 5)
	sum_tif.write(ef4, 6)
	sum_tif.write(ef5, 7)

	sum_tif.close()

# get_track_lat_lon returns a data frame of track data for the specified hurricane
# if the maximum enhanced Fujita value exceeds a specified value.
#   hur_id - hurricane id
#   fuj_min - minimum enhanced Fujita value
#   tt - data frame of track data
#   kk - data frame of summary data
#   returns a data frame of track data

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

#' hurrecon_set_path sets the path for the current set of model runs.
#'   hur_path - path for current model runs
#'   no return value

def hurrecon_set_path(hur_path):
	if hur_path == "":
		print("\nNeed to enter a path\n")
		sys.exit()

	elif os.path.exists(hur_path) == False:
		print("\nPath does not exist\n")
		sys.exit()

	os.chdir(hur_path)
	print("Path set to", hur_path)

# hurrecon_create_land_water creates a land-water raster file in GeoTiff 
# format from vector boundary files in shapefile format. The land-water file
# (land_water.tif) is assumed to be aligned with lines of latitude and 
# longitude.  Boundary files are assumed to be named boundary.* on the gis 
# subdirectory. The land-water file is created on the input subdirectory.
#   nrows - number of rows
#   ncols - number of columns
#   xmn - minimum longitude (degrees)
#   xmx - maximum longitude (degrees)
#   ymn - minimum latitude (degrees)
#   ymx - maximum latitude (degrees)
#   no return value

def hurrecon_create_land_water(nrows, ncols, xmn, xmx, ymn, ymx):
	# get current working directory
	cwd = os.getcwd()

	# open boundaries file
	boundaries_file = cwd + "/gis/boundaries.shp"
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

	print("Cell height =", round(cell_height) , "kilometers")
	print("Cell width  =", round(cell_width), "kilometers\n")

# hurrecon_reformat_hurdat2 reformats a HURDAT2 file from the National 
# Hurricane Center for use with the HURRECON model. The input file is assumed
# to be in space-delimited text format. Two output files are created on the
# input subdirectory: hurdat2_ids.csv contains a list of hurricanes including 
# id, name, number of positions, and peak sustained wind speed (meters/second).
# hurdat2_tracks.csv contains full track information for each hurricane 
# from HURDAT2 plus columns for standard datetime and Julian day with fraction.
#   hurdat2_file - name of HURDAT2 file
#   path - optional path for HURDAT2 file
#   no return value

def hurrecon_reformat_hurdat2(hurdat2_file, path=""):
	# output files
	ids_file = "hurdat2_ids.csv"
	tracks_file = "hurdat2_tracks.csv"

	if path != "":
		if path[len(path)-1] != "/":
			path = path + "/"
		hurdat2_file = path + hurdat2_file
		ids_file = path + ids_file
		tracks_file = path + tracks_file

	# read hurdat2 file
	file_in = open(hurdat2_file, 'r') 
	hurdat = file_in.readlines() 
	nlines = len(hurdat)

	# close hurdat2 file
	file_in.close()

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
	
		# store values in ids lists
		ids_hur_id_list.append(hur_id)
		ids_name_list.append(name)
		ids_positions_list.append(positions)
		ids_wind_peak_list.append(wind_peak)

		# report progress
		x = round(line_num*100/nlines)
		if x % 10 == 0:
			print("          ", end="")
			print("\r", x, "%", end="")

	# create data frames
	ids_cols = ['hur_id', 'name', 'positions', 'wind_peak']
	ids = pd.DataFrame(data=list(zip(ids_hur_id_list, ids_name_list, ids_positions_list, ids_wind_peak_list)), 
		columns=ids_cols)

	tracks_cols = ['hur_id', 'name', 'date_time', 'jd', 'status', 'latitude', 'longitude', 'wind_max']

	tracks= pd.DataFrame(data=list(zip(tracks_hur_id_list, tracks_name_list, tracks_date_time_list, 
		tracks_jd_list, tracks_status_list, tracks_latitude_list, tracks_longitude_list, tracks_wind_max_list)), 
		columns=tracks_cols)

	# save to file
	ids.to_csv(ids_file, index=False)
	tracks.to_csv(tracks_file, index=False)

	# display number of storms
	print("\nNumber of storms = ", len(ids))
	print("Number of observations = ", len(tracks))

# hurrecon_extract_tracks extracts hurricane ids and tracks from the two
# files created by hurrecon_reformat_hurdat2 (hurdat2_ids.csv and 
# hurdat2_tracks.csv). The geographic window used to select hurricanes is 
# set by the land-water file and optionally extended by the margin parameter.
# Selection begins by identifying all positions in the window where the hurricane
# has "HU" (hurricane) status in HURDAT2.  If at least one such position exists,
# the track is extended to include one position before and one position after
# the window, if possible. If the resulting track contains at least two positions
# and the maximum sustained wind speed equals or exceeds wind_min, the track 
# is included.
#   margin - optional extension of the geographic window set by the
#     land-water file (degrees)
#   wind_min - minimum value of maximum sustained wind speed 
#     (meters/second)
#   no return value

def hurrecon_extract_tracks(margin=0, wind_min=33):
	# get current working directory
	cwd = os.getcwd()

	# output files
	ids_file = cwd + "/input/ids.csv"
	tracks_file = cwd + "/input/tracks.csv"

	# read hurdat2 ids file
	hurdat2_ids_file = cwd + "/input/hurdat2_ids.csv"
	check_file_exists(hurdat2_ids_file)
	ii = pd.read_csv(hurdat2_ids_file)
	ii_rows = len(ii)

	# read hurdat2 tracks file
	hurdat2_tracks_file = cwd + "/input/hurdat2_tracks.csv"
	check_file_exists(hurdat2_tracks_file)
	tt = pd.read_csv(hurdat2_tracks_file)
	tt_rows = len(tt)

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

	# process files
	for i in range(0, ii_rows):
		# get hurricane id & name
		hur_id = ii.hur_id[i]
		name = ii.name[i]

		index = np.where((tt.hur_id == hur_id) & (tt.latitude >= lat_min) & (tt.latitude <= lat_max) & (tt.longitude >= lon_min) & (tt.longitude <= lon_max) & (tt.status == "HU"))[0].tolist()

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

		# report progress
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

	# save to file
	ids.to_csv(ids_file, index=False)
	tracks.to_csv(tracks_file, index=False)

	# display number of storms
	print("\nNumber of storms = ", len(ids))
	print("Number of observations = ", len(tracks))


### MODELING FUNCTIONS ####################################

# hurrecon_model_site calculates wind speed (meters/second), gust speed 
# (meters/second), wind direction (degrees), and enhanced Fujita scale wind 
# damage for a given hurricane and site. If width is True, the radius of 
# maximum wind (rmw) and profile exponent (s_par) for the given hurricane are 
# used, if available.  If save is True, results are saved to a CSV file on 
# the site subdirectory; otherwise results are returned as a data frame. 
# If timing is True, the total elasped time is displayed.
#   hur_id - hurricane id
#   site_name - name of site
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   save - whether to save results to a CSV file
#   timing - whether to display total elapsed time
#   returns a data frame of results if save is False

def hurrecon_model_site(hur_id, site_name, width=False, time_step=1, save=True, 
	timing=True):

	# record total elapsed time if timing is True
	if (timing == True):
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
	track = read_hurricane_track_file(hur_id)

	# interpolate hurricane location & max wind speed
	modeled = interpolate_hurricane_location_max_wind(track, time_step)

	# interpolate hurricane speed & bearing
	modeled = interpolate_hurricane_speed_bearing(track, modeled)

	# calculate range & bearing from site to hurricane center
	modeled = calculate_site_range_bearing(modeled, site_latitude, site_longitude)

	# calculate wind speed, wind direction & enhanced Fujita scale at site
	modeled = calculate_wind_speed_direction(modeled, inflow_angle, cover_type, rmw, 
		s_par, asymmetry_factor, friction_factor, gust_factor)

	# add standard date & time
	modeled = add_standard_date_time(modeled)

	# display total elapsed time
	if timing == True:
		print(format_time_difference_ms(start_time, time.time()), " ms")

	# output
	if save == True:
		# save modeled data to CSV file
		modeled_file = cwd + "/site/" + hur_id + " " + site_name + ".csv"
		modeled.to_csv(modeled_file, index=False)
		print("Saving to", modeled_file)
	else:
		# return modeled data as data frame
		return modeled

# hurrecon_model_site_all creates a table of peak values for all hurricanes
# for a given site. If width is True, the radius of maximum wind (rmw) and 
# profile exponent (s_par) for the given hurricane are used, if available. 
# If save is True, results are saved to a CSV file on the site subdirectory; 
# otherwise results are returned as a data frame. If timing is True, the 
# total elasped time is displayed.
#   site_name - name of site
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   save - whether to save results to a CSV file
#   timing - whether to display total elapsed time
#   returns a data frame of results if save is False

def hurrecon_model_site_all(site_name, width=False, time_step=1, save=True, 
	timing=True):

	# get current working directory
	cwd = os.getcwd()

	# read ids file
	ids_file = cwd + "/input/ids.csv"
	check_file_exists(ids_file)
	ii = pd.read_csv(ids_file)
	ii_rows = len(ii)

	# initialize lists
	snam_list = []
	hnam_list = []
	dt_list = []
	wdir_list = []
	wspd_list = []
	gust_list = []
	ef_list = []
	ef0_list = []
	ef1_list = []
	ef2_list = []
	ef3_list = []
	ef4_list = []
	ef5_list = []

	# record total elasped time if timing is True
	if timing == True:
		start_time = time.time()

	# get peak values for each hurricane
	for i in range(0, ii_rows):
		# get hurricane name
		hur_id = ii.hur_id[i]

		# get modeled output
		mm = hurrecon_model_site(hur_id, site_name, width, time_step, save=False,
			timing=False)

		# get peak values
		pk = get_peak_values(hur_id, site_name, mm)
		
		snam_list.append(pk.site_name[0])
		hnam_list.append(pk.hur_id[0])
		dt_list.append(pk.date_time[0])
		wdir_list.append(pk.wind_dir[0])
		wspd_list.append(pk.wind_spd[0])
		gust_list.append(pk.gust_spd[0])
		ef_list.append(pk.ef_sca[0])
		ef0_list.append(pk.ef0[0])
		ef1_list.append(pk.ef1[0])
		ef2_list.append(pk.ef2[0])
		ef3_list.append(pk.ef3[0])
		ef4_list.append(pk.ef4[0])
		ef5_list.append(pk.ef5[0])

		# report progress
		if timing == True:
			x = round(i*100/ii_rows)
			if x % 10 == 0:
				print("          ", end="")
				print("\r", x, "%", end="")

	if timing == True:
		elapsed_time = format_time_difference_hms(start_time, time.time())
		print("\r", elapsed_time, "\n", end="")

	# create data frame for peak values
	peak_values_cols = ['site_name', 'hur_id', 'date_time', 'wind_dir', 'wind_spd',
		'gust_spd', "ef_sca", "ef0", "ef1", "ef2", "ef3", "ef4", "ef5"]

	peak_values = pd.DataFrame(data=list(zip(snam_list, hnam_list, dt_list, wdir_list, wspd_list, 
		gust_list, ef_list, ef0_list, ef1_list, ef2_list, ef3_list, ef4_list, ef5_list)), 
		columns=peak_values_cols)

	# output
	if (save == True):
		# save modeled data to CSV file
		site_peak_file = cwd + "/site/" + site_name + " Peak Values.csv"
		peak_values.to_csv(site_peak_file, index=False)
		print("Saving to", site_peak_file)
	else:
		# return modeled data as data frame
		return peak_values

# hurrecon_model_region calculates peak wind speed (meters/second), enhanced
# Fujita scale, wind direction (degrees), cardinal wind direction, gale wind
# duration (minutes), and hurricane wind duration (minutes) for a given 
# hurricane over a region. If width is True, the radius of maximum wind (rmw) 
# and profile exponent (s_par) for the given hurricane are used, if available.
# If no value is provided for time step, the time step is calculated. If water 
# is False, results are calculated for land areas only. If save is True, results
# are saved as a GeoTiff file; otherwise results are returned as a list of 
# rasters. If timing is True, the total elasped time is displayed.
#   hur_id - hurricane id
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   water - whether to caculate results over water
#   save - whether to save results to a GeoTiff file
#   timing - whether to display total elapsed time
#   returns a list of 6 rasters if save is False

def hurrecon_model_region(hur_id, width=False, time_step="", water=False, save=True, 
	timing=True):
	
	# get current working directory
	cwd = os.getcwd()

	# get time step if necessary
	if time_step == "":
		time_step = get_time_step()

	if timing:
		print("Time step =", time_step, "minutes")
 
 	# read hurricane track file
	track = read_hurricane_track_file(hur_id)

	# interpolate hurricane location & max wind speed
	modeled = interpolate_hurricane_location_max_wind(track, time_step)

	# interpolate hurricane speed & bearing
	modeled = interpolate_hurricane_speed_bearing(track, modeled)

	# get modeled values over region
	peak = get_regional_peak_wind(hur_id, modeled, width, time_step, water, timing)

	# output
	if (save == True):
		# read land-water file
		land_water_file = cwd + "/input/land_water.tif"
		check_file_exists(land_water_file)
		land_water = rio.open(land_water_file)

		# open GeoTiff file
		hur_tif_file = cwd + "/region/" + hur_id + ".tif"

		profile = land_water.profile
		profile.update(dtype='int16', nodata=-9999, count=6)
        
		hur_tif = rio.open(hur_tif_file, 'w', **profile)

		hur_tif.write(peak[0], 1)
		hur_tif.write(peak[1], 2)
		hur_tif.write(peak[2], 3)
		hur_tif.write(peak[3], 4)
		hur_tif.write(peak[4], 5)
		hur_tif.write(peak[5], 6)

		hur_tif.close()

	else:
		# return modeled values as a list of 6 arrays
		return peak

# hurrecon_model_region_all calculates peak wind speed (meters/second), 
# enhanced Fujita scale, wind direction (degrees), cardinal wind direction, 
# duration of gale winds (minutes), and duration of hurricane winds (minutes) 
# over a region for all hurricanes. If width is True, the radius of maximum 
# wind (rmw) and profile exponent (s_par) for the given hurricane are used, 
# if available. If no value is provided for time step, the time step is 
# calculated. If water is False, results are calculated for land areas only.
# Results for each hurricane are saved in a GeoTiff file on the region 
# subdirectory. Summary results for all hurricanes (summary.tif, summary.csv)
# are also calculated and saved to the region subdirectory.
#   width - whether to use width parameters for the specified hurricane
#   time_step - time step (minutes)
#   water - whether to calculate results over water
#   no return value

def hurrecon_model_region_all(width=False, time_step="", water=False):
	# get current working directory
	cwd = os.getcwd()

	# get time step if necessary
	if time_step == "":
		time_step = get_time_step()
	
	print("Time step =", time_step, "minutes")

	# read ids file
	ids_file = cwd + "/input/ids.csv"
	check_file_exists(ids_file)
	ii = pd.read_csv(ids_file)
	ii_rows = len(ii)

	# record total elasped time
	start_time = time.time()

	# get regional estimate for each hurricane
	for i in range(0, ii_rows):
		# get hurricane name
		hur_id = ii.hur_id[i]

		# report progress
		x = round(i*100/ii_rows)
		print("          ", end="")
		print("\r", x, "%", end="")

		# generate & save regional results
		hurrecon_model_region(hur_id, width, time_step, water, save=True, timing=False)

	# generate & save regional summary files
	get_regional_summary()

	# display total elapsed time
	elapsed_time = format_time_difference_hms(start_time, time.time())
	print("\r", elapsed_time, "\n", end="")


### SUMMARIZING FUNCTIONS #################################

# hurrecon_summarize_land_water displays features of the current land-water
# file (land_water.tif) in the console.
#   no return value

def hurrecon_summarize_land_water():
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

	print("Rows:", nrows, "  Columns:", ncols)
	print("Latitude:", ymn, "to", ymx, "degrees")
	print("Longitude:", xmn, "to", xmx, "degrees")
	print("Cell height:", round(cell_height), "kilometers")
	print("Cell width:", round(cell_width), "kilometers")
	print("Time Step:", time_step, "minutes")

# hurrecon_summarize_tracks displays features of the current ids file 
# (ids.csv) in the console.
#   no return value

def hurrecon_summarize_tracks():
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

	print("Number of storms =", ii_rows)
	print("Number of positions =", positions_total)
	print("Minimum peak wind =", wind_peak_min, "m/s")
	print("Maximum peak wind =", wind_peak_max, "m/s")

# hurrecon_summarize_site displays peak values for a given hurricane
# and site in the console.
#   hur_id - hurricane id
#   site_name - name of site
#   no return value

def hurrecon_summarize_site(hur_id, site_name):
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
	print(modeled_name)

	print("PEAK:", pk.date_time[0])
	print("Wind dir:", int(round(pk.wind_dir[0])), "deg")
	print("Wind spd:", int(round(pk.wind_spd[0])), "m/s")
	print("Gust spd:", int(round(pk.gust_spd[0])), "m/s")

	if pk.ef0[0] > 0:
		print("EF0:", round(pk.ef0[0], 1), "hours")
	if pk.ef1[0] > 0:
		print("EF1:", round(pk.ef1[0], 1), "hours")
	if pk.ef2[0] > 0:
		print("EF2:", round(pk.ef2[0], 1), "hours")
	if pk.ef3[0] > 0:
		print("EF3:", round(pk.ef3[0], 1), "hours")
	if pk.ef4[0] > 0:
		print("EF4:", round(pk.ef4[0], 1), "hours")
	if pk.ef5[0] > 0:
		print("EF5:", round(pk.ef5[0], 1), "hours")


### PLOTTING FUNCTIONS ####################################

# hurrecon_plot_site_ts creates a time-series plot of wind speed, gust 
# speed, or wind direction as a function of datetime for a given 
# hurricane and site. Optional start and end datetimes may be specified.
#   hur_id - hurricane id
#   site_name - name of site
#   start_datetime - optional start datetime in format YYYY-MM-DD hh:mm
#   end_datetime - optional end datetime in format YYYY-MM-DD hh:mm
#   var - variable to plot (wind_speed, gust_speed, or wind_direction)
#   no return value

def hurrecon_plot_site_ts(hur_id, site_name, start_datetime='', end_datetime='', 
	var="wind_speed"):
	
	# get current working directory
	cwd = os.getcwd()

	# get enhanced Fujita wind speeds & colors
	ef = get_fujita_wind_speeds()
	ef_col = get_fujita_colors()

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

	# get axis labels
	x_var = "dt"
	x_label = "Datetime (UTC)"

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
		print("\nvar must be wind_speed, gust_speed, or wind_direction\n")
		sys.exit()

	# subset by data range
	if (start_datetime != ""):
		sdate = start_datetime
	else:
		sdate = mm.dt[0]

	if (end_datetime != ""):
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

	plot_name = hur_id + " " + site_name

	# create plot
	plt.figure(figsize=(10,10))

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

	plt.title(plot_name, fontsize=20)
	plt.xlabel(x_label, fontsize=18)
	plt.ylabel(y_label, fontsize=18)
	plt.legend(loc='upper right')
	plt.show()
	plt.clf()

# hurrecon_plot_site_xy creates a scatter plot of wind speed or gust speed
# as a function of wind direction for a given hurricane and site. If adjust 
# is True, 360 degrees are substracted from wind direction values greater 
# than 180. Optional start and end datetimes may be specified.
#   hur_id - hurricane id
#   site_name - name of site
#   start_datetime - optional start datetime in format YYYY-MM-DD hh:mm
#   end_datetime - optional end datetime in format YYYY-MM-DD hh:mm
#   var - variable to plot (wind_speed or gust_speed)
#   adjust - whether to subtract 360 degrees from wind direction.
#   no return value

def hurrecon_plot_site_xy(hur_id, site_name, start_datetime='', end_datetime='', 
	var="wind_speed", adjust=False):

	# get current working directory
	cwd = os.getcwd()

	# get enhanced Fujita wind speeds & colors
	ef = get_fujita_wind_speeds()
	ef_col = get_fujita_colors()

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

	# get axis labels
	x_var = "wind_dir"
	x_label = "Wind Direction (deg)"

	if var == "wind_speed":
		y_var = "wind_spd"
		y_label = "Wind Speed (m/s)"
	elif var == "gust_speed":
		y_var = "gust_spd"
		y_label = "Gust Speed (m/s)"
	else:
		print("\nvar must be wind_speed or gust_speed\n")
		sys.exit()

	# adjust wind direction
	if adjust == True:
		wdir_list = []

		for i in range(0, mm_rows):
			wdir = mm.wind_dir[i]
			if wdir > 180:
				wdir = wdir - 360

			wdir_list.append(wdir)

		mm['wind_dir2'] = wdir_list
		x_var = "wind_dir2"

	if (start_datetime != ""):
		sdate = start_datetime
	else:
		sdate = mm.dt[0]

	if (end_datetime != ""):
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

	plot_name = hur_id + " " + site_name

	# create plot
	plt.figure(figsize=(10,10))

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

	plt.title(plot_name, fontsize=20)
	plt.xlabel(x_label, fontsize=18)
	plt.ylabel(y_label, fontsize=18)
	plt.legend(loc='upper right')
	plt.show()
	plt.clf()

# hurrecon_plot_site_all creates a time-series plot of peak values 
# for all hurricanes for a given site. Optional start and end years
# may be specified.
#   site_name - name of site
#   start_year - optional start year
#   end_year - optional end year
#   var - variable to plot (wind_speed, gust_speed, or wind_direction)
#   no return value

def hurrecon_plot_site_all(site_name, start_year='', end_year='', 
	var="wind_speed"):

	# get current working directory
	cwd = os.getcwd()

	# get enhanced Fujita wind speeds & colors
	ef = get_fujita_wind_speeds()
	ef_col = get_fujita_colors()

	# read data
	peak_file = cwd + "/site/" + site_name + " Peak Values.csv"
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
		print("\nvar must be wind_speed, gust_speed, or wind_direction\n")
		sys.exit()

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

	plot_name = site_name

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

	plt.title(plot_name, fontsize=20)
	plt.xlabel(x_label, fontsize=18)
	plt.ylabel(y_label, fontsize=18)
	plt.legend(loc='lower left')
	plt.show()
	plt.clf()

# hurrecon_plot_region creates regional plots of enhanced Fujita scale, 
# peak wind speed, wind direction, cardinal wind direction, gale wind 
# duration, and hurricane wind duration for a given hurricane.
#   hur_id - hurricane id
#   var - variable to plot (wind_speed, fujita_scale, wind_direction,
#     wind_compass, gale_duration, hurricane_duration)
#   no return value

def hurrecon_plot_region(hur_id, var="fujita_scale"):
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
	boundaries_file = cwd + "/gis/boundaries.shp"
	shapefile = fiona.open(boundaries_file, "r")
	features = [feature["geometry"] for feature in shapefile]
	shapefile.close()

	# get hurricane track
	track_file = cwd + "/input/tracks.csv"
	check_file_exists(track_file)
	tt = pd.read_csv(track_file)
	index = np.where(tt.hur_id == hur_id)[0].tolist()

	lat_list = [tt.latitude[index[i]] for i in range(0, len(index))]
	lon_list = [tt.longitude[index[i]] for i in range(0, len(index))]

	# get colormaps with white background
	viridis = plt.get_cmap('viridis')
	viridis.set_under('white') 	

	rainbow = plt.get_cmap('rainbow')
	rainbow.set_under('white') 	

	# create plot
	if var == "fujita_scale":
		if np.amax(hur_tif.read(2)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='grey', linewidth=0.8)
			plt.title(hur_id + ' Fujita Scale')
			img = hur_tif.read(2)
			plt.imshow(img, cmap=rainbow, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])
			cbar.set_ticklabels(['', 'None', 'EF0', 'EF1', 'EF2', 'EF3', 'EF4', 'EF5'])
			show((hur_tif, 2), cmap=rainbow, vmin=0.9)
			plt.clf()
		else:
			print("No Fujita values\n")

	elif var == "wind_speed":
		if np.amax(hur_tif.read(1)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='grey', linewidth=0.8)
			plt.title(hur_id + ' Wind Speed (m/s)')
			img = hur_tif.read(1)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			show((hur_tif, 1), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No wind speed\n")

	elif var == "wind_direction":
		if np.amax(hur_tif.read(3)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='grey', linewidth=0.8)
			plt.title(hur_id + ' Wind Direction (deg)')
			img = hur_tif.read(3)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			show((hur_tif, 3), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No wind direction\n")

	elif var == "wind_compass":
		if np.amax(hur_tif.read(4)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='grey', linewidth=0.8)
			plt.title(hur_id + ' Wind Direction')
			img = hur_tif.read(4)		
			plt.imshow(img, cmap=rainbow, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
			cbar.set_ticklabels(['', 'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
			show((hur_tif, 4), cmap=rainbow, vmin=0.9)
			plt.clf()	
		else:
			print("No wind compass\n")

	elif var == "gale_duration":
		if np.amax(hur_tif.read(5)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='grey', linewidth=0.8)
			plt.title(hur_id + ' Gale Winds (min)')
			img = hur_tif.read(5)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			show((hur_tif, 5), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No gale winds\n")
	
	elif var == "hurricane_duration":
		if np.amax(hur_tif.read(6)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.plot(lon_list, lat_list, color='grey', linewidth=0.8)
			plt.title(hur_id + ' Hurricane Winds (min)')
			img = hur_tif.read(6)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			show((hur_tif, 6), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No hurricane winds\n")

	else:
		print("\nvar must be wind_speed, fujita_scale, wind_direction, wind_compass, gale_duration, or hurricane_duration\n")
		sys.exit()

	# close GeoTiff file
	hur_tif.close()

# hurrecon_plot_region_all creates regional plots of maximum enhanced Fujita
# value and number of storms for each enhanced Fujita value for all hurricanes.
#   var - variable to plot (efmax, ef0, ef1, ef2, ef3, ef4, or ef5)
#   tracks - whether to also plot hurricane tracks
#   no return value

def hurrecon_plot_region_all(var="efmax", tracks=False):
	# get current working directory
	cwd = os.getcwd()

	# read GeoTiff file
	sum_tif_file = cwd + "/region/" + "summary.tif"
	check_file_exists(sum_tif_file)
	sum_tif = rio.open(sum_tif_file)

	# get extent
	lat_min = sum_tif.bounds.bottom
	lat_max = sum_tif.bounds.top
	lon_min = sum_tif.bounds.left
	lon_max = sum_tif.bounds.right

	# read boundaries file
	boundaries_file = cwd + "/gis/boundaries.shp"
	shapefile = fiona.open(boundaries_file, "r")
	features = [feature["geometry"] for feature in shapefile]
	shapefile.close()

	# get colormaps with white background
	viridis = plt.get_cmap('viridis')
	viridis.set_under('white') 	

	rainbow = plt.get_cmap('rainbow')
	rainbow.set_under('white') 	

	# get hurricane tracks
	if tracks == True:
		ids_file = cwd + "/input/ids.csv"
		check_file_exists(ids_file)
		ii = pd.read_csv(ids_file)

		track_file = cwd + "/input/tracks.csv"
		check_file_exists(track_file)
		tt = pd.read_csv(track_file)

		summary_file = cwd + "/region/summary.csv"
		check_file_exists(summary_file)
		kk = pd.read_csv(summary_file)

	# create plot
	if var == "efmax":
		if np.amax(sum_tif.read(1)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title('Maximum Fujita Scale')
			img = sum_tif.read(1)		
			plt.imshow(img, cmap=rainbow, vmin=0.9)
			cbar = plt.colorbar(shrink=0.3)
			cbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])
			cbar.set_ticklabels(['', 'None', 'EF0', 'EF1', 'EF2', 'EF3', 'EF4', 'EF5'])
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 0
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						plt.plot(xx.longitude, xx.latitude, color='grey', linewidth=0.8)
			show((sum_tif, 1), cmap=rainbow, vmin=0.9)
			plt.clf()	
		else:
			print("No Fujita values\n")

	elif var == "ef0":
		if np.amax(sum_tif.read(2)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title('Fujita Scale 0 Storms')
			img = sum_tif.read(2)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 0
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						plt.plot(xx.longitude, xx.latitude, color='grey', linewidth=0.8)
			show((sum_tif, 2), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No F0 values\n")

	elif var == "ef1":
		if np.amax(sum_tif.read(3)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title('Fujita Scale 1 Storms')
			img = sum_tif.read(3)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 1
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						plt.plot(xx.longitude, xx.latitude, color='grey', linewidth=0.8)
			show((sum_tif, 3), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No F1 values\n")

	elif var == "ef2":
		if np.amax(sum_tif.read(4)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title('Fujita Scale 2 Storms')
			img = sum_tif.read(4)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 2
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						plt.plot(xx.longitude, xx.latitude, color='grey', linewidth=0.8)
			show((sum_tif, 4), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No F2 values\n")

	elif var == "ef3":
		if np.amax(sum_tif.read(5)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title('Fujita Scale 3 Storms')
			img = sum_tif.read(5)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 3
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						plt.plot(xx.longitude, xx.latitude, color='grey', linewidth=0.8)
			show((sum_tif, 5), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No F3 values\n")
	
	elif var == "ef4":
		if np.amax(sum_tif.read(6)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title('Fujita Scale 4 Storms')
			img = sum_tif.read(6)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 4
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						plt.plot(xx.longitude, xx.latitude, color='grey', linewidth=0.8)
			show((sum_tif, 6), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No F4 values\n")

	elif var == "ef5":
		if np.amax(sum_tif.read(7)) > 0:
			fig, ax = plt.subplots(figsize=(15, 15))
			plt.axis([lon_min, lon_max, lat_min, lat_max])
			plt.xlabel('Longitude (degrees)')
			plt.ylabel('Latitude (degrees)')
			patches = [PolygonPatch(feature, edgecolor="black", facecolor="none", linewidth=1) for feature in features]
			ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))
			plt.title('Fujita Scale 5 Storms')
			img = sum_tif.read(7)		
			plt.imshow(img, cmap=viridis, vmin=0.9)
			plt.colorbar(shrink=0.3)
			if tracks == True:
				for i in range(0, len(ii)):
					hur_id = ii.loc[i, "hur_id"]
					fuj_min = 5
					xx = get_track_lat_lon(hur_id, fuj_min, tt, kk)
					if len(xx) > 0:
						plt.plot(xx.longitude, xx.latitude, color='grey', linewidth=0.8)
			show((sum_tif, 7), cmap=viridis, vmin=0.9)
			plt.clf()	
		else:
			print("No F5 values\n")

	else:
		print("\nvar must be efmax, ef0, ef1, ef2, ef3, ef4, or ef5\n")
		sys.exit()

	# close GeoTiff file
	sum_tif.close()
