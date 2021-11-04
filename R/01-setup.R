### Set up 

library(here)
library(reticulate)
library(tidyverse)
library(UpSetR)

# load reticulate and python client
library(reticulate)
synapse <- reticulate::import("synapseclient")
syn <- synapse$Synapse()

# login to Synapse
syn$login()

### Get de-identified Rosmap data
rosmap_obj <- syn$get("syn26436148")
rosmap <- read_csv(rosmap_obj$path)


