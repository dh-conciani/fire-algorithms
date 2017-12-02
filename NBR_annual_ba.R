# Algorithm to extract burned area from NBR index and temporal/spatial fire regime
# Dhemerson Conciani (dh.conciani@gmail.com)

# Input:
# NBR index (single band raster) - (calculated by you or provided by ESPA / USGS)

# Read libraries
library (data.table) # Make data frame with rasters
library (lubridate)  # Extract temporal variables
library (rgdal)      # Raster Tools
library (raster)     # Raster Tools
library (stringr)    # Parse
library (tools)      # Data manipulation tools
#library (xlsx)      # Write final table

###  LOADING DATA AND SETTING PARAMETERS ###
# Read NBR index rasters
nbr        <- list.files(path = 'H:/nbr/index/itirapina', pattern = 'sr_nbr.tif$', full.names = T)
rs_list    <- lapply(nbr,raster)                    # Transform rasters into list
nbr_names  <- file_path_sans_ext(basename(nbr))     # Extract filenames
nbr.lenght <- length(nbr)                           # Extract lenght
base_lists = 1:nbr.lenght

# Parse dates from NBR's and transform into factor variable
parse  <- sapply(strsplit(nbr_names, split='_', fixed=TRUE), function(x) (x[3]))  # Parse brute date string 
date   <- ymd(parse)                                                              # Transform string into date
date   <- as.character(date)                                                      # Transform date in strings separated by '-'
years  <- sapply(strsplit(date, split='-', fixed=TRUE), function(x) (x[1]))       # Extract year string, the first sep. 
months <- sapply(strsplit(date, split='-', fixed=TRUE), function(x) (x[2]))       # Extract month string, the second sep. 
years   = as.factor (years)
months  = as.factor (months)
rm (parse, date)

# Import .shp's from interest areas (covering the NBR index raster) and write strings
eecl      = readOGR (dsn='H:/shps/limites', layer='dissolve')      # Ecological Station of Itirapina
eei       = readOGR (dsn='H:/shps/limites', layer='Experimental')  # Experimental Station of Itirapina
buffer    = readOGR (dsn='H:/shps/buffer7km', layer='mask')        # 7 km buffer zone

# Create reclassification matrix (1= burned area; 0= unburned area)
m        <- c(-9999.9999, -1000, 1,  -1000, 20000, 0)   # Define thresholds and classes to reclassify (this is the principal parameter to reclass NBR)
rclmat   <- matrix(m, ncol=3, byrow=TRUE)               # Make reclassification table

# Create annual burned area reclassificatin matrix
m.annual      <- c(1, 99999, 1)                        # equal or > than 1 = (1)
rclmat.annual <- matrix(m.annual, ncol=3, byrow=TRUE)  # make annual reclassification table 

# Define functions
reproject_raster  = function (x) projectRaster (x, crs = ref_proj)        # Reproject raster list using a reference cordinate system
crop_to_area      = function (x) crop (x, buffer)
reclass_nbr       = function (x) reclassify (x, rclmat)                   # Reclassify rasters and make binary images
#reclass_annual_ba = function (x) reclassify (x, rclmat.annual)           # Reclassify classiified NBR to binary annual Burned Area
eecl_mask         = function (x) mask (x, eecl)                           # Extract data only to EEcl
eei_mask          = function (x) mask (x, eei)                            # Extract data only to EEI
buffer_mask       = function (x) mask (x, buffer)                         # Extract data only to Buffer
area_fun          = function (x) sum(values(x)==1, na.rm=TRUE) * 30 * 30  # Calc burned area based on binary images
### END DEFINITIONS ###

### START OF CALCS STEP ### 
# 1. Reproject rasters to reference projection
ref_proj  <- crs(buffer)                                # Extract projection reference
reproj_list <- lapply (rs_list, reproject_raster)       # Reproject rasters to reference projection
rm (rs_list, ref_proj)

# Crop raster to selected area
croped_list = lapply (reproj_list, crop_to_area)        #Crop to Buffer area
rm(reproj_list)

# 2. Reclassify rasters to binary images (1= burned area; 0= unburned area)
rcl_list <- lapply (croped_list, reclass_nbr)             # Perform reclass and make binary raster's
rm (m, rclmat) 

# 3. Stack classified rasters and reclassify to annual burned area
stacked_nbr          <- stack (rcl_list)                         # stack classified rasters
annual_ba            <- stackApply(stacked_nbr,years, fun=sum)   # sum stack, parameters = input_stack, indices (in this case, year factors), function
plot (annual_ba)                                                 # check somatory
rcl_annual_ba        <- reclassify (annual_ba, rclmat.annual)    # reclassify to annual burned area
names(rcl_annual_ba) <- levels (years)                           # name rasters with year information
plot(rcl_annual_ba)                                              # check annual burned area
annual_ba            <- unstack (rcl_annual_ba)                  # unstack 
rm (stacked_nbr, rcl_annual_ba)

# 3. Extract only interest areas 
eecl_list    <- lapply (annual_ba, eecl_mask)    # Extract to EEcl
eei_list     <- lapply (annual_ba, eei_mask)     # Extract to EEI
buffer_list  <- lapply (annual_ba, buffer_mask)  # Extract to Buffer

#### STATISTICS #####
# 4. Calc burned area (in hectares) and percentage of burned area
# a. Ecological Station of Itirapina
local      = rep("EEcl",length(eecl_list))                                            # Local name              
area       = rep(2678, length(eecl_list))                                             # Area of local
ba_eecl   <- lapply (eecl_list, area_fun)                                              # List burned area in EEcl
             burned_area= as.matrix (ba_eecl)                                          # Make a matrix
             burned_area= as.numeric (burned_area)                                     # Convert to numeric variable
             df_eecl    = data.frame (levels(years), burned_area, local, area)         # Structure data frame
             df_eecl$burned_area = df_eecl$burned_area * 0.0001                        # Convert m2 to ha
             df_eecl[, "percent_ba"] <- df_eecl[, "burned_area"] / df_eecl[, "area"]   # Normalize burned area (burned area / local area)
             df_eecl$percent_ba= df_eecl$percent_ba * 100                              # Multiply *100 (percentage)

# b. Experimental Station of Itirapina
local      = rep("EEI",length(eei_list))                                               # Local name              
area       = rep(3131, length(eei_list))                                               # Area of local
ba_eei   <- lapply (eei_list, area_fun)                                                # List burned area in EEcl
             burned_area= as.matrix (ba_eei)                                           # Make a matrix
             burned_area= as.numeric (burned_area)                                     # Convert to numeric variable
             df_eei    = data.frame (levels(years), burned_area, local, area)          # Structure data frame
             df_eei$burned_area = df_eei$burned_area * 0.0001                          # Convert m2 to ha
             df_eei[, "percent_ba"] <- df_eei[, "burned_area"] / df_eei[, "area"]      # Normalize burned area (burned area / local area)
             df_eei$percent_ba= df_eei$percent_ba * 100                                # Multiply *100 (percentage)
             
# c. Buffer of Itirapina
local      = rep("Buffer",length(buffer_list))                                              # Local name              
area       = rep(42340, length(buffer_list))                                                # Area of local
ba_buffer   <- lapply (buffer_list, area_fun)                                               # List burned area in EEcl
             burned_area= as.matrix (ba_buffer)                                             # Make a matrix
             burned_area= as.numeric (burned_area)                                          # Convert to numeric variable
             df_buffer    = data.frame (levels(years), burned_area, local, area)            # Structure data frame
             df_buffer$burned_area = df_buffer$burned_area * 0.0001                         # Convert m2 to ha
             df_buffer[, "percent_ba"] <- df_buffer[, "burned_area"] / df_buffer[, "area"]  # Normalize burned area (burned area / local area)
             df_buffer$percent_ba= df_buffer$percent_ba * 100                               # Multiply *100 (percentage)
             
# 5. Make a single data frame with all data
df <- rbind (df_eecl, df_eei, df_buffer)  
setwd ("H:/nbr/results")
write.xlsx(df, file="itirapina.xlsx", sheetName="sheet1")


plot (df$levels.years., df$percent_ba)
