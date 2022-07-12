## HURRECON Model

The HURRECON model estimates wind speed, wind direction, enhanced Fujita scale
wind damage, and duration of EF0 to EF5 winds as a function of hurricane location 
and maximum sustained wind speed. Results may be generated for a single site or 
an entire region. Hurricane track and intensity data may be imported directly 
from the US National Hurricane Center's HURDAT2 database.

HURRECON is available in both R (HurreconR) and Python (HurreconPython) versions. 
The model is an updated version of the original HURRECON model written in Borland 
Pascal for use with Idrisi (see below for details).

Please note: both versions are currently under development and subject to change.

## Getting Started

Here are the basic steps for using the model. Please see below for more details.

1. Download the R or Python version of the model from GitHub.
2. Create a directory for a set of related model runs with subdirectories for input 
and output files (as described below).
3. Create a site file (<i>sites.csv</i>) with geographic coordinates for one or more 
study sites.
4. Download or create geographic and political boundary shapefiles for the desired 
geographic region. The coordinate system should be latitude/longitude (degrees). 
Rename these files so the first name of each file is "boundaries".
5. Create a land-water file (<i>land_water.tif</i>) for the region.
6. Create a parameter file (<i>parameters.csv</i>) with parameters for all hurricanes 
and (optionally) for particular hurricanes.
7. Create an input hurricane track file (<i>input_tracks.csv</i>) for the geographic 
region. If desired, this file can be created directly from HURDAT2.
8. Run the model to create site and regional estimates. Use the plot functions to 
view model results.

## Details

All user functions begin with "hurrecon". Input CSV files have a single header
line that contains variable names (see details below). Model results are 
generated for two cover types: water and land. Datetimes are assumed to be 
UTC (GMT). Raster files are in GeoTiff format. Geographic coordinates are assumed 
to be latitude/longitude (degrees). The following measurement units are used throughout:

```
bearing, direction - degrees
distance - kilometers
speed - meters/second
time - minutes or hours
```

The user specifies a path (<i>hur_path</i>) for a given set of model runs. Input
and output files are stored on the following subdirectories of this path:

```
hur_path/input
hur_path/region
hur_path/region-all
hur_path/region-dt
hur_path/site
hur_path/site-all
hur_path/vector
```

The <i>input</i> subdirectory contains input files. The <i>site</i> and <i>region</i> 
subdirectories contain site and regional output files, respectively. Shapefiles that 
contain geographic and political boundaries for viewing regional results are stored 
on the <i>vector</i> subdirectory.

The following input files are required:

```
sites.csv
boundaries.*
land_water.tif
parameters.csv
input_tracks.csv
```

All input files (except boundary files) are located on the <i>input</i> subdirectory.

The sites file contains the name, location, and cover type (water = 1, land = 2)
of each study site. Variables: <i>site_name</i>, <i>latitude</i>, <i>longitude</i>, 
<i>cover_type</i>.

The boundary files are vector shapefiles that are used for creating maps of regional
results. These files are located on the <i>vector</i> subdirectory.

The land-water file is a raster GeoTiff file that specifies the cover type 
(water = 1, land = 2) for each cell across a region. The geographic coordinates 
and the number of rows and columns of the land-water file are used to set the 
geographic window and spatial resolution for regional modeling.

The parameters file contains model parameters (radius of maximum winds and scaling
parameter) for all hurricanes and (optionally) for individual hurricanes. Variables:
<i>hur_id</i>, <i>rmw</i>, <i>s_par</i>. This file must contain at least one record 
(<i>hur_id</i> = ALL) that specifies the default values of <i>rmw</i> and <i>s_par</i>. 
Values typically range from 20 to 100 km for <i>rmw</i> and from 1.2 to 1.5 for <i>s_par</i>, 
depending on the region.

The input tracks file contains location and maximum wind speed for each position of
each hurricane for a given set of model runs. Variables: <i>hur_id</i>, <i>name</i>, 
<i>date_time</i>, <i>jd</i>, <i>status</i>, <i>latitude</i>, <i>longitude</i>, 
<i>wind_max</i>.

The input tracks file may be created directly from HURDAT2. Use the 
<i>hurrecon_reformat_hurdat2</i> function to reformat a HUTDAT2 file as 
<i>hurdat2_tracks.csv</i>, rename this file to <i>input_tracks.csv</i>, and copy this 
file to the <i>input</i> subdirectory.

The <i>hurrecon_extract_tracks</i> function is used to extract the data needed for a 
particular set of model runs. This function uses the input tracks file and the 
land-water file to create input files (<i>ids.csv</i>, <i>tracks.csv</i>, 
<i>tracks-all.csv</i>) required for the <i>hurrecon_model</i> functions.

Examples of input files may be found on the <i>inst/input</i> subdirectory (R) or 
<i>data</i> subdirectory (Python).

To run the model, create the above directories, copy the input files to their
respective subdirectories, and run <i>hurrecon.R</i> (R) or <i>hurrecon.py</i> (Python). 

The R version may also be installed as an R package using the devtools package:

```
devtools::install_github("hurrecon-model/HurreconR")
```

This has the advantage of providing readily accessible help messages for each 
function.

## Model Functions

```
hurrecon_reformat_hurdat2

hurrecon_set_path
hurrecon_get_path

hurrecon_create_land_water
hurrecon_extract_tracks

hurrecon_model_site
hurrecon_model_site_all
hurrecon_model_region
hurrecon_model_region_dt
hurrecon_model_region_all

hurrecon_summarize_land_water
hurrecon_summarize_tracks
hurrecon_summarize_site

hurrecon_plot_site
hurrecon_plot_site_all
hurrecon_plot_tracks
hurrecon_plot_region
hurrecon_plot_region_dt
hurrecon_plot_region_all
```

The <i>hurrecon_reformat_hurdat2</i> function reformats data from HURDAT2 
for use with HURRECON. This is normally a one-time operation for a
given version of HURDAT2. Hurricane IDs in HURDAT2 are reformatted 
to facilitate sorting by year (e.g. AL031935 becomes AL1935-03).

The <i>hurrecon_set_path</i> function sets the path for the current set of 
model runs. The <i>hurrecon_get_path</i> function returns the current path. 
Use <i>hurrecon_set_path</i> before using other functions.

The <i>hurrecon_create_land_water</i> function creates a land-water raster file
in GeoTiff format using the specified minimum & maximum latitude & longitude,
the number of rows & columns, and vector boundary files in shapefile 
format used to set the cover type of each cell. The land-water file is used
by other functions to get the spatial parameters required for regional 
modeling.  The <i>hurrecon_extract_tracks</i> function extracts data from an input 
tracks file (which may be derived from HURDAT2) for use with a particular 
land-water file. Optional parameters may be used to broaden the geographic area
or adjust the minimum hurricane intensity when selecting hurricane tracks.

The <i>hurrecon_model</i> functions generate output for a single hurricane and a 
single site (all datetimes), all hurricanes for a single site (peak values), 
a single hurricane for a specified geographic region (peak values or specified
datetime), and all hurricanes for a specified geographic region (peak values).
If save is TRUE (default), results are written to the <i>site</i>, <i>site-all</i>, 
<i>region</i>, <i>region-all</i>, or <i>region-dt</i> subdirectory as CSV or GeoTiff
files. The default time step for site results is 1 minute. The default time step 
for regional results is calculated as the time required to traverse one cell in
the vertical direction at 20 meters per second, rounded to one of these values:
1, 2, 3, 5, 10, 15, 30, or 60 minutes.

The <i>hurrecon_summarize_land_water</i> function displays information about the current
land-water file. The <i>hurrecon_summarize_tracks</i> function displays information about
the current track files. The <i>hurrecon_summarize_site</i> function displays peak values 
for a single hurricane and a single site.

The <i>hurrecon_plot_site</i> functions create time-series and scatter plots for a single 
hurricane and time-series plots for all hurricanes for a given site. The 
<i>hurrecon_plot_tracks</i> function creates a map of the land-water file with selected 
hurricane tracks. The <i>hurrecon_plot_region</i> functions create maps of regional results 
for a single hurricane or for all hurricanes.

## Examples

Sample commands for the 1935 Florida Keys hurricane and Miami FL:

```
hurrecon_reformat_hurdat2(hurdat2_file="hurdat2-1851-2021-041922.txt")
[copy hurdat2_tracks.csv to input_tracks.csv on input subdirectory]

hurrecon_set_path("c:/hurrecon/east_20km")
hurrecon_get_path()

hurrecon_create_land_water(nrows=150, ncols=180, xmn=-100, xmx=-59, ymn=23, ymx=50)
hurrecon_extract_tracks(wind_min=70)

hurrecon_model_site("AL1935-03", "Miami FL")
hurrecon_model_site_all("Miami FL")
hurrecon_model_region("AL1935-03")
hurrecon_model_region_dt("AL1935-03", "1935-09-03T12:00")
hurrecon_model_region_all()

hurrecon_summarize_land_water()
hurrecon_summarize_tracks()
hurrecon_summarize_site("AL1935-03", "Miami FL")

hurrecon_plot_site("AL1935-03", "Miami FL")
hurrecon_plot_site_all("Miami FL")
hurrecon_plot_tracks()
hurrecon_plot_region("AL1935-03")
hurrecon_plot_region_dt("AL1935-03", "1935-09-03T12:00")
hurrecon_plot_region_all()
```

## Model Equations

The sustained wind speed (Vs) at any point P in the northern hemisphere is estimated as:

```
[1] Vs = F[Vm - S(1 - sin T)Vh/2] * Sqrt[(Rm/R)^B * exp(1 - (Rm/R)^B)]

where:

F = scaling parameter for the effects of friction (water = 1.0, land = 0.8)
Vm = maximum sustained wind speed over water anywhere in hurricane
S = scaling parameter for asymmetry due to forward motion of hurricane (1.0)
T = clockwise angle between forward motion of hurricane and radial line to P
Vh = forward speed of hurricane
Rm = radius of maximum winds
R = radial distance from hurricane center to point P
B = scaling parameter controlling shape of wind profile curve
```

The peak wind gust speed (Vg) at point P is estimated as:

```
[2] Vg = G * Vs

where:

G = gust factor (water = 1.2, land = 1.5)
```

The wind direction (D) in degrees at point P is estimated as:

```
[3] D = Az - 90 - I

where:

Az = azimuth from point P to hurricane center
I = cross-isobar inflow angle (water = 20 degrees, land = 40 degrees)
```

In the southern hemisphere:

```
T = counterclockwise angle between forward motion of hurricane and radial line to P
D = Az + 90 + I
```

For more details, see publications below.

## History

The original HURRECON model was written in Borland Pascal and depended on Idrisi 
for spatial visualization. The model was used in published studies of the ecological 
impacts of historical hurricanes in New England and Puerto Rico:

* Boose, E. R., Chamberlin, K. E., Foster, D. R. 2001. Landscape and regional impacts 
of hurricanes in New England. Ecological Monographs 71: 27-48.
doi:10.1890/0012-9615(2001)071[0027:LARIOH]2.0.CO;2.

* Boose, E. R., Serrano, M. I., Foster, D. R. 2004. Landscape and regional impacts of 
hurricanes in Puerto Rico. Ecological Monographs 74: 335-352. doi:10.1890/02-4057.

New features in the updated version of HURRECON include support for: (1) estimating 
wind damage on the enhanced Fujita scale, (2) importing hurricane track and intensity 
data directly from HURDAT2, (3) creating a land-water file with user-selected 
geographic coordinates and spatial resolution, and (4) creating plots of site and 
regional results.

The model equations for estimating wind speed and direction, including parameter values
for inflow angle, friction factor, and wind gust factor (over land and water), are 
unchanged from the original HURRECON model.

