##### AUTHOR:     JILL THORR, M.S.
##### SCRIPT:     MODEL VALIDATION: HYCOM-TSIS (GOMb0.04) MAX DEPTH ESTIMATES
##### VERSION:    1
##### DATE:       5 MARCH 2025
##### OBJECTIVE:  RUN EOF ANALYSIS ON THE GRIDDED HYCOM PRODUCTS (NCODA + TSIS)
##### 			      TO ASSESS MODES OF SUBSURFACE TEMPERATURE VARIABILITY
##### 			      WITHIN THE LATEX CONTINENTAL SHELF (1/25 RESOLUTION) 
##### NOTES:* EOF USED FOR CHAPTER 3
#####			  * SEE HYCOM_TSISreanalysis_DownloadMaxDepthEstimates.R
#####			   	FOR HOW DATA WERE ACCESSED AND DOWNLOADED

# Set your working directory to the repo root before running.
# All paths below are relative to that root.

# SET THE VARIABLE
VAR <- 'watertemp' 
#VAR <- 'salinity'
#VAR <- 'wvelocity'


# SET UP ENVIRONMENT
options(scipen = 999)
options(digits = 8)    # required for precise joining

# OPEN LIBRARIES
library(ggplot2,     quietly = TRUE)         # for plots
library(dplyr,       quietly = TRUE)         # for data manipulation
library(lubridate,   quietly = TRUE)	       # for formatting month
library(foreach,     quietly = TRUE)         # for running in parallels
library(doParallel,  quietly = TRUE)         # for running in parallels
library(parallel,    quietly = TRUE)         # for running in parallels

# OPEN SOURCE FILES
source('./src/function_unregistercluster.R') # unregister cluster after parallel processing
source('./src/function_loadRData.R')         # specify name when loading RData 
source('./src/function_SVD.R')               # singular value decomposition of high-dimensional data
source('./src/function_calcEOFandPC.R')      # calculate Empirical Orthogonal functions of decomposed data
fontfam <- 'DejaVu LGC Serif Condensed'      # for plot fonts


# LOAD DATA --------------------------------------------------------------------

# Load spatial reference data
llx       <- read.csv('./data/HYCOM_MaxDepth/HYCOMTSIS_spRefwithMaxLayerZ_pointsInLATEXcontShelf.csv')
llid      <- read.csv('./data/HYCOM_MaxDepth/HYCOMTSIS_spRefwithMaxLayerZ.csv')
llid$keep <- ifelse(llid$TSISLLID %in% llx$TSISLLID, 'Y', 'N')   # T == points in study domain
length(unique(llid[which(llid$keep == 'Y'),]$TSISLLID))          # [1] 9081
dim(llid); dim(llx)
rm(llx)


# Load the information about the layer with max depth (of the 40 vertical layers)
ref    <- read.csv('./data/HYCOM_MaxDepth/maxLayer_HYCOMTSIS_data.csv'); dim(ref)        # [1] 219  109
ref.rn <- read.csv('./data/HYCOM_MaxDepth/maxLayer_HYCOMTSIS_rownames.csv'); dim(ref.rn) # [1] 219 <-- row names
ref.cn <- read.csv('./data/HYCOM_MaxDepth/maxLayer_HYCOMTSIS_colnames.csv'); dim(ref.cn) # [1] 109 <-- column names
rownames(ref) <- ref.rn[,1]
colnames(ref) <- ref.cn[,1]


# Set variable specific labels and colors for EOF and PC plots
if(VAR == 'salinity')   {varName <- 'Bottom Salinity';    varUnit <- ''}
if(VAR == 'watertemp')  {varName <- 'Bottom Temperature'; varUnit <- '(\u00B0C)'}
if(VAR == 'wvelocity')  {varName <- 'Vertical Velocity';  varUnit <- '(m/s)'}


# Load the time and space references for the HYCOM model data
dtref     <- read.csv('./data/HYCOM_MaxDepth/HYCOMNCODAandTSIS_1996to2022_timeSteps.csv') %>% mutate(DATE = as.Date(DATE)) # time-steps
llid      <- read.csv('./data/HYCOM_MaxDepth/HYCOMTSIS_LLIDref_pointsInLATEXcontShelf.csv')[,1:4]                          # spatial references -- drop keep column
llidref   <- read.csv('./data/HYCOM_MaxDepth/HYCOMNCODAandTSIS_1996to2022_LLIDref.csv')                                    # spatial reference between TSIS and NCODA models   
ncodallid <- read.csv('./data/HYCOM_MaxDepth/HYCOMexp50p1_LLIDref.csv') %>% dplyr::rename(lon=LON,lat=LAT,NCODALLID = LLID)# spatial reference to the standard LAT/LON used for plotting
llid      <- llid[which(llid$LLID %in% llidref$TSISLLID),]                                                                 # filter to just the points in continental shelf using tsisllid
llid      <- llid %>% dplyr::rename(TSISLLID = LLID)  %>%                                                                  # spatial reference that combines the two models
             left_join(llidref, by = 'TSISLLID')      %>% 
             dplyr::select(-c(lon,lat,Z_layer))       %>%
             left_join(ncodallid, by = 'NCODALLID')


# Load detrended residuals and merge lists of each cell (LLID) to matrix where:
# columns = LLID (cell)
# rows    = dtref (day)
f     <- './data/HYCOM_MaxDepth/DetrendedResidualsOutput/detrendedresiduals_HYCOMNCODAandTSIS_1996to2022_nbc_BottomTemp.RData'
resid <- do.call(cbind, loadRData(f[1]))


# ANALYSES ---------------------------------------------------------------------

# Decompose the data then calculate EOF, PCs, and variance explained by each mode
if(!(tolower(eoffn <- paste0('EOFresults_1996-2022_',gsub(' ','',VAR),'.RData')) %in% tolower(list.files('./data/HYCOM_MaxDepth/')))) {
  dcmp   <- f.SVD(A = resid)                                  # decompose the detrended residuals with singular value decomposition
  eofRes <- f.calcEOFandPC(svdresults = dcmp, A = resid)      # calculate the EOF from the SVD results
  eofRes$llid <- llid                                         # spatially assign EOF coordinates with NCODA lat/lon values
  eofRes$dtref<- dtref                                        # add dates associated with PC scores
  save(eofRes, file = paste0('./data/HYCOM_MaxDepth/',eoffn)) # save locally
} else {
  eofRes <- loadRData(paste0('./data/HYCOM_MaxDepth/',eoffn)) # open EOF results if already saved locally
}

# PLOTS ------------------------------------------------------------------------

# Specify colors for EOF maps
hex       <- c('#FF0000', '#FFA500', '#FFFF00', '#008000', '#9999FF', '#000066') 

# Specify color for PC line time-series plots (light = raw values, dark = gaussian smoothed)
if(VAR == 'salinity') {linecolor <- c('#bde4e4','#357363')}  # PC line color = green [light HEX, dark HEX]
if(VAR == 'watertemp'){linecolor <- c('#caa5eb','#181d62')}  # PC line color = purple
if(VAR == 'wvelocity'){linecolor <- c('#feb3a4','#fc4e2a')}  # PC line color = orange

# Function to plot EOFs
f.plotEOF <- function(X) {
  # X = EOF Mode (saved as columns)
  ggplot()                                                                                                     +
    geom_tile(aes(x = llid$lon, y = llid$lat, fill = eofRes$EOF[,X]))                                          +
    scale_y_continuous(limits = c(26, 29.8), breaks = seq(26, 30, 1), expand = c(0,0))                         +
    scale_x_continuous(limits = c(-97.49, -89), breaks = seq(-97, -89, 1), expand = c(0,0))                    +
    scale_fill_gradientn(colors = rev(hex), limits = c(-1, 1))                                                 +
    geom_vline(xintercept = seq(-97.5,-89, 0.25), lty = 'dashed', alpha = 0.2)                                 +
    geom_hline(yintercept = seq(26,29.8,0.25), lty = 'dashed', alpha = 0.2)                                    +
    theme_classic(base_size = 14)                                                                              + 
    theme(text = element_text(family = fontfam), axis.title = element_blank(), legend.title = element_blank()) +
    labs(title = paste0('[',varName,']', ' Mode ', X, ' (Variance Explained = ', round(eofRes$VarExp[X], 1), '%)'))
}

# Function to plot PCs with raw and smoothed values
f.plotPC <- function(X) {
  # X = PC scores associated with the Mode (saved as rows)
  ggplot()                       															              +
    geom_line(aes(x = dtref$DATE, y = eofRes$PC[X,]), color = linecolor[1], lwd = 0.2)                   +
    geom_line(aes(x = dtref$DATE, y = smoother::smth.gaussian(eofRes$PC[X,], window = 200)), color = linecolor[2])+
    scale_x_date(date_labels = "%Y", date_breaks = '1 year', expand = expansion(add = c(100,7)))      +
    theme_classic(base_size = 14) 																	  + 
    theme(text         = element_text(family = fontfam),
          axis.text.x  = element_text(angle  = 90,hjust = 1,vjust = 0.5))+
    labs(title = '', y = 'Principal Component', x = 'Year')
}

# Plot EOF and PC Modes that explain at least 5% of the variation in the data 
outplot   <- list(); i = 1
for(i in 1:length(which(eofRes$VarExp > 5))) {outplot[[i]] <- cowplot::plot_grid(f.plotEOF(X = i),f.plotPC(X = i), nrow = 1)}

# Save plots locally in a single png file
phigh <- length(outplot)*4; pwide <- 17  # specify file height and width; each plot should be 4 in high 
pfile <- paste0('./data/HYCOM_MaxDepth/Plots/HYCOMTSIS_1996-2022_',VAR,'_EOFandPCplots_varexp5percent_GausSmoother_200daywindow_wGrid.png')
png(pfile, width = pwide, height = phigh, units = 'in',res = 600)
if(length(outplot) == 3) {print(cowplot::plot_grid(outplot[[1]],outplot[[2]],outplot[[3]], nrow = 3))} 
if(length(outplot) == 4) {print(cowplot::plot_grid(outplot[[1]],outplot[[2]],outplot[[3]],outplot[[4]], nrow = 4))} 
dev.off()
