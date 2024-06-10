# step 1: create dir
dir.create("./6035873")
dir.create("./GSE141259")
dir.create("./GSE152501")
dir.create("./GSE163465")
dir.create("./GSE180420")
dir.create("./GSE186986")
dir.create("./GSE200843")
dir.create("./GSE205037")
dir.create("./GSE205690")

dir.create("./6035873/rawData")
dir.create("./GSE141259/rawData")
dir.create("./GSE152501/rawData")
dir.create("./GSE163465/rawData")
dir.create("./GSE180420/rawData")
dir.create("./GSE186986/rawData")
dir.create("./GSE200843/rawData")
dir.create("./GSE205037/rawData")
dir.create("./GSE205690/rawData")

dir.create("./6035873/RData")
dir.create("./GSE141259/RData")
dir.create("./GSE152501/RData")
dir.create("./GSE163465/RData")
dir.create("./GSE180420/RData")
dir.create("./GSE186986/RData")
dir.create("./GSE200843/RData")
dir.create("./GSE205037/RData")
dir.create("./GSE205690/RData")

# step 2: download raw data
# raw data was downloaded from GEO and Zenodo (Table S2) to each "rawData" dir

# step 3: prepare scripts
file.copy(from = "./sc_script_for_6035873_liver.R", to = "./6035873/sc_script_for_6035873_liver.R")
file.copy(from = "./sc_script_for_GSE141259_lung.R", to = "./GSE141259/sc_script_for_GSE141259_lung.R")
file.copy(from = "./sc_script_for_GSE152501_airway.R", to = "./GSE152501/sc_script_for_GSE152501_airway.R")
file.copy(from = "./sc_script_for_GSE163465_heart.R", to = "./GSE163465/sc_script_for_GSE163465_heart.R")
file.copy(from = "./sc_script_for_GSE180420_kidney.R", to = "./GSE180420/sc_script_for_GSE180420_kidney.R")
file.copy(from = "./sc_script_for_GSE186986_skin.R", to = "./GSE186986/sc_script_for_GSE186986_skin.R")
file.copy(from = "./sc_script_for_GSE200843_joints.R", to = "./GSE200843/sc_script_for_GSE200843_joints.R")
file.copy(from = "./sc_script_for_GSE205037_spinal_cord.R", to = "./GSE205037/sc_script_for_GSE205037_spinal_cord.R")
file.copy(from = "./sc_script_for_GSE205690_muscle.R", to = "./GSE205690/sc_script_for_GSE205690_muscle.R")

file.remove("./sc_script_for_6035873_liver.R",
            "./sc_script_for_GSE141259_lung.R",
            "./sc_script_for_GSE152501_airway.R",
            "./sc_script_for_GSE163465_heart.R",
            "./sc_script_for_GSE180420_kidney.R",
            "./sc_script_for_GSE186986_skin.R",
            "./sc_script_for_GSE200843_joints.R",
            "./sc_script_for_GSE205037_spinal_cord.R",
            "./sc_script_for_GSE205690_muscle.R")
# step 4: run scripts in each dir
# step 5: run "sc_script_for_plot.R"

