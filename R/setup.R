#set up
library(here)
library(reticulate)
library(tidyverse)

# load reticulate and python client
synapse <- reticulate::import("synapseclient")
syn <- synapse$Synapse()

# login to Synapse
syn$login()
