#This script combines the two files output from the Event calling script, also sorts the events and adds some more information.
#The Spo11 hits from Pan et al are used to calculate Spo11 hits per event.

#load required package
library(plyr)

#source the function
source("event_sorter.R")

# Read and process file for calculating Spo11 hits per event 
# Only run these lines once:

Pan <- read.delim("FullMap.Cer3_Pan_HA_1_4h_c.txt")
TotalHits=sum(Pan$Watson+Pan$Crick)
Pan$TotalHpM=(Pan$Watson+Pan$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)

##read in the Event table to be sorted##



##read in the Event table to be sorted##
##Change the filenames as required##
b<-read.delim("Events_TW1_1500") # the name of the file describing the events from NEW event caller
u<-read.delim("UEvents_TW1_1500") #the name of the file describing the events from OLD event caller - for Unclassified events

event_sorter(b, u, Pan, "NEvents_TW1_1500")


##read in the Event table to be sorted##
##Change the filenames as required##
b<-read.delim("Events_TN3_1500") # the name of the file describing the events from NEW event caller
u<-read.delim("UEvents_TN3_1500") #the name of the file describing the events from OLD event caller - for Unclassified events

event_sorter(b, u, Pan, "NEvents_TN3_1500")


##read in the Event table to be sorted##
##Change the filenames as required##
b<-read.delim("Events_TN4_1500") # the name of the file describing the events from NEW event caller
u<-read.delim("UEvents_TN4_1500") #the name of the file describing the events from OLD event caller - for Unclassified events

event_sorter(b, u, Pan, "NEvents_TN4_1500")


##read in the Event table to be sorted##
##Change the filenames as required##
b<-read.delim("Events_TN6_1500") # the name of the file describing the events from NEW event caller
u<-read.delim("UEvents_TN6_1500") #the name of the file describing the events from OLD event caller - for Unclassified events

event_sorter(b, u, Pan, "NEvents_TN6_1500")


##read in the Event table to be sorted##
##Change the filenames as required##
b<-read.delim("Events_TN8_1500") # the name of the file describing the events from NEW event caller
u<-read.delim("UEvents_TN8_1500") #the name of the file describing the events from OLD event caller - for Unclassified events

event_sorter(b, u, Pan, "NEvents_TN8_1500")


##read in the Event table to be sorted##
##Change the filenames as required##
b<-read.delim("Events_OMN2_1500") # the name of the file describing the events from NEW event caller
u<-read.delim("UEvents_OMN2_1500") #the name of the file describing the events from OLD event caller - for Unclassified events

event_sorter(b, u, Pan, "NEvents_OMN2_1500")

##read in the Event table to be sorted##
##Change the filenames as required##
b<-read.delim("Events_OMN6_1500") # the name of the file describing the events from NEW event caller
u<-read.delim("UEvents_OMN6_1500") #the name of the file describing the events from OLD event caller - for Unclassified events

event_sorter(b, u, Pan, "NEvents_OMN6_1500")