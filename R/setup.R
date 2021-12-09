#set up
library(here)
library(reticulate)
library(tidyverse)
library(fastDummies)
library(ComplexUpset)

# load reticulate and python client
synapse <- reticulate::import("synapseclient")
syn <- synapse$Synapse()

# login to Synapse
syn$login()
