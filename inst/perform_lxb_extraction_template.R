#!/usr/bin/Rscript

# Wrapper script to perform the extraction of lxb data, simply calls functions defined in the 'extract_lxb_data.R' script in the correct order

# Get the plate variable which describes the experimental setup
# Specific to the dataset at hand, process the annotations to a plate format and provide the correspondence between the bead IDs ('region') and the proteins they identify ('analyte')
region_analytes = matrix(c(14, "SMAD2",
                          18, "GSK3b",
                          27, "MEK",
                          34, "JNK",
                          35, "p90RSK",
                          36, "p38",
                          38, "ERK",
                          44, "EGFR",
                          46, "mTOR",
                          53, "p53",
                          54, "PI3Kp85",
                          56, "cJUN",
                          67, "IkBa",
                          74, "RPS6",
                          75, "AKT"
                           ), nrow=2)
region = region_analytes[1,]
analyte = region_analytes[2,]

library("gdata") # To read xls files
plate = matrix( "", nrow=8, ncol=12, dimnames=list(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) ) # Empty plate
# Load the annotations of the wells
plate[c("A", "B", "C"),1] = "BSA"
plate[c("A", "B", "C"),2] = "NGF"
plate[c("A", "B"),3] = "EGF"
plate[c("A", "B", "C"),4] = "BDNF"
plate[c("A", "B"),5] = "IGF1"
plate[c("A", "B"),6] = "TGFb"
plate[c("C", "D"),6] = "TNFa"
plate[c("C", "D"),c(3, 5, 7)] = "blank"
plate["D", 1] = "BioRad_TNF"
plate["D", 2] = "BioRad_phosphatase"
plate["D", 4] = "BioRad_EGF"

# Define controls
plate[grepl("BSA", plate)] = "control"

library(extractlxb)
library(lxb)
# Retrieve the folder to analyse
print("Extracting treatments")
esetup = extractExperimentalSetup(plate)

# Load the data
print("Loading the lxb files")
if (!exists("cargs")) {
    cargs = commandArgs(trailingOnly=T)
} else if (is.character(cargs)) {
    cargs = unlist(strsplit(cargs, " "))
}
# Retrieve the name of the folder from the command line
do_plots = FALSE
for (arg in cargs) {
    if (grepl("_lxb$", arg)) {
        lxb_folder = paste0(arg, "/")
    } else if (grepl("_lxb/$", arg)) {
        lxb_folder = arg
    } else if (grepl("--plot", arg)) {
        do_plots = TRUE
    }
}
if (!exists("lxb_folder")) {
    if (length(cargs) > 0) {
        lxb_folder = cargs[1]
    } else {
        lxb_folder = "170629_widr_lxb" # For lazy use of the script
    }
}
if (!grepl("/$", lxb_folder)) {
    lxb_folder = paste0(lxb_folder, "/")
}
full_data = readLxb(paste0(lxb_folder, "*.lxb"), text=FALSE)
lxb_name = gsub("[0-9]{6}_", "", lxb_folder)
lxb_name = gsub("_lxb", "", gsub("/$", "", lxb_name))

print("Processing bead data")
beads = readBeads(full_data, region, analyte, esetup$controls, esetup$blanks, esetup$wptb, lxb_folder, esetup$externals, do_plots)
print("Bead data processed")
if (file.exists("post_processing.R")) {
    print("Post processing the data")
    source("post_processing.R")
}
writeMIDAS(beads, lxb_name, esetup$inhibitors, esetup$stimulators, lxb_name)
print(paste0("Data writen in MIDAS format in folder ", lxb_name))

pdf(paste0(lxb_name, "/", "beads_distribution_", lxb_name, ".pdf"), width=10)
analyseBeadDistributions(full_data, region, analyte, esetup$treatments, bn_limit=10)
dev.off()

