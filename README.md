## HURRECON Model

The HURRECON model estimates wind speed, wind direction, enhanced Fujita scale
wind damage, duration of gale winds, and duration of hurricane winds as a function 
of hurricane location and maximum sustained wind speed. Results may be generated for 
a single site or an entire region. Hurricane track and intensity data may be imported
directly from the US National Hurricane Center's HURDAT2 database.

HURRECON is available in both R (HurreconR) and Python (HurreconPython) versions. 
The model is an updated version of the original HURRECON model written in Borland 
Pascal for use with Idrisi (see below for details).

Please note: both versions are currently under development and subject to change.

## Getting Started

All user functions begin with "hurrecon". Input CSV files have a single header
line that contains variable names (see input file examples). Model results are 
generated for two cover types: water and land. All datetimes are assumed to be 
UTC (GMT). The following measurement units are used throughout:

```{r}
latitude, longitude, bearing, direction - degrees
distance - kilometers
speed - meters/second
time - minutes or hours
```

The user specifies a directory (hur_dir) for a given set of model runs. Input
and output files are stored on the following subdirectories of this directory:

```{r}
hur_dir/input
hur_dir/region
hur_dir/region-dt
hur_dir/region-all
hur_dir/site
hur_dir/site-all
hur_dir/vector
```

The input subdirectory contains input files. The site and region subdirectories
contain site and regional output files, respectively. Shapefiles that contain
geographic and political boundaries for viewing regional results are stored on the 
vector subdirectory.

The following input files are required:

```{r}
sites.csv
ids.csv
tracks.csv
parameters.csv
land-water.tif
boundary.*
```

All input files (except boundary files) are located on the input subdirectory.

The sites file contains the name, location, and cover type (water = 1, land = 2)
of each study site. Variables: site_name, latitude, longitude, cover_type.

The ids file contains the hurricane id, name, number of positions, and peak 
sustained wind speed for each hurricane. Variables: hur_id, name, positions,
wind_peak.

The tracks file contains location and maximum wind speed for each location of
each hurricane. Variables: hur_id, name, date_time, jd, status, latitude,
longitude, wind_max.

The ids and tracks files may be created directly from HURDAT2. Use the 
hurrecon_reformat_hurdat2 function to reformat the HURDAT2 file as hurdat2_ids.csv
and hurdat2_tracks.csv, copy these files to the input directory, and use the
hurrecon_extract_tracks function to create ids.csv and tracks.csv.

The parameters file contains model parameters (radius of maximum winds and profile
exponent) for all hurricanes and (optionally) for individual hurricanes. Variables:
hur_id, rmw, s_par. This file must contain at least one record (hur_id = ALL) that 
specifies the default values of rmw and s_par.

The land-water file is a raster GeoTiff file that specifies the cover type 
(water = 1, land = 2) for each cell across a region. The geographic coordinates 
and number of rows and columns of the land-water file are used to set the geographic
window and spatial resolution for regional modeling.

The boundary files are vector shapefiles that are used for creating maps of regional
results. Files must be renamed so the first name of each file is "boundary".
These files are located on the vector subdirectory.

Examples of input files may be found on the inst/extdata subdirectory (R) or data
subdirectory (Python).

To run the model, create the above directories, copy the input files to their
respective subdirectories, and source hurrecon.R (R) or hurrecon.py (Python). 

The R version may also be installed as an R package using the devtools package:

```{r}	
devtools::install_github("hurrecon-model/HurreconR")
```

This has the advantage of providing readily accessible help messages for each 
function.

## Model Functions

```{r}	
hurrecon_reformat_hurdat2

hurrecon_set_path

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
hurrecon_plot_region
hurrecon_plot_region_dt
hurrecon_plot_region_all
```

The hurrecon_reformat_hurdat2 function reformats data from HURDAT2 
for use with HURRECON. This is normally a one-time operation for a
given version of HURDAT2.

The hurrecon_set_path function sets the current working directory to 
the desired directory for the current set of model runs.

The hurrecon_create_land_water function creates a land-water raster file
in GeoTiff format using the specified minimum & maximum latitude & longitude,
the number of rows & columns, and vector boundary files in shapefile 
format used to set the cover type of each cell. The land-water file is used
by other functions to get the spatial parameters required for regional 
modeling.  The hurrecon_extract_tracks function extracts data from 
the reformatted HURDAT2 file for use with a particular land-water file.
Optional parameters may be used to broaden the geographic area or adjust 
the minimum hurricane intensity when selecting hurricane tracks.

The hurrecon_model functions generate output for a single hurricane and a 
single site (site), all hurricanes for a single site (site_all), one hurricane
for the current geographic region (region), one hurricane for the current 
geographic region at a specified datetime (region-dt), or all hurricanes for
the current geographic region (region_all). If save is TRUE (default), results
are written to the site, site-all, region, region-dt, or region-all subdirectory 
as CSV or GeoTiff files, respectively. The default time step for site results is 
1 minute. The default time step for regional results is calculated as the time 
required to traverse one cell in the vertical direction at 20 meters per second, 
rounded to one of these values: 1, 2, 3, 5, 10, 15, 30, or 60 minutes.

The hurrecon_summarize_land_water function displays features of the current
land-water file. The hurrecon_summarize_tracks function displays features of
the current ids and tracks files. The hurrecon_summarize_site function displays
summary values for a single hurricane and a single site.

The hurrecon_plot_site functions create a time series plot of wind speed, gust
speed, or wind direction or a scatter plot of wind or gust speed as a function of
wind direction for a single hurricane, or a time series plot of peak
wind speed, gust speed, or wind direction for all hurricanes for a given site.
The hurrecon_plot_region functions create maps of regional results for a 
single hurricane, a single hurricane at a specified datetime, or for all hurricanes.

## Examples

Sample commands for the 1935 Florida Keys hurricane and Miami FL:

```{r}
hurrecon_reformat_hurdat2(hurdat2_file="hurdat2-1851-2019-042820.txt")

hurrecon_set_path("c:/hurrecon/r/east")

hurrecon_create_land_water(nrows=100, ncols=120, xmn=-100, xmx=-59, ymn=23, ymx=50)
hurrecon_extract_tracks(wind_min=70)

hurrecon_model_site("AL031935", "Miami FL")
hurrecon_model_site_all("Miami FL")
hurrecon_model_region("AL031935")
hurrecon_model_region_dt("AL031935", "1935-09-03T05:27")
hurrecon_model_region_all()

hurrecon_summarize_land_water()
hurrecon_summarize_tracks()
hurrecon_summarize_site("AL031935", "Miami FL")

hurrecon_plot_site("AL031935", "Miami FL")
hurrecon_plot_site_all("Miami FL")
hurrecon_plot_region("AL031935")
hurrecon_plot_region_dt("AL031935", "1935-09-03T05:27")
hurrecon_plot_region_all()
```

## History

The original HURRECON model was written in Borland Pascal and depended on Idrisi 
for spatial visualization. The model was used in published studies of the ecological 
impacts of historical hurricanes in New England and Puerto Rico:

* Boose, E. R., Chamberlin, K. E., Foster, D. R. 2001. Landscape and regional impacts 
of hurricanes in New England. Ecological Monographs 71: 27-48.

* Boose, E. R., Serrano, M. I., Foster, D. R. 2004. Landscape and regional impacts of 
hurricanes in Puerto Rico. Ecological Monographs 74: 335-352.

New features in the updated version of HURRECON include support for: (1) estimating 
wind damage on the enhanced Fujita scale, (2) importing hurricane track and intensity 
data directly from HURDAT2, (3) creating a land-water file with user-selected 
geographic coordinates and spatial resolution, and (4) creating plots of site and 
regional results.

The model equations for estimating wind speed and direction are unchanged from the 
original HURRECON model and are explained in the publications above.

