library(optimx)        
library(FD)       
library(snow) 
library(parallel)
library(BioGeoBEARS)

setwd("~/Cirsium/biogeo")
getwd()

calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

trfn = "new.cirsium.tre"
moref(trfn)
tr = read.tree(trfn)
tr

geogfn = "geo.biome.data"
moref(geogfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
max_range_size = 5

# Run DEC
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object$master_table
BioGeoBEARS_run_object$speedup = TRUE         
BioGeoBEARS_run_object$use_optimx = "GenSA"  
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    
#BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "cirsium_biomes_DEC.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  load(resfn)
  resDEC = res
}



#######################################################
# PDF plots
#######################################################
pdffn = "Cirsium_biomes_DEC.pdf"
pdf(pdffn, width=12, height=16)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Cirsium (biomes)"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), 
                                plotwhat="text", label.offset=0.05, tipcex=0.7, statecex=0.7, splitcex=0.6, 
                                titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                                tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
                         label.offset=0.05, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8,
                         plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, 
                         tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

# Run DIVALIKE
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object$master_table
BioGeoBEARS_run_object$on_NaN_error = -1e50 
BioGeoBEARS_run_object$speedup = TRUE         
BioGeoBEARS_run_object$use_optimx = "GenSA"  
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Cirsium_biomes_DIVALIKE.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}



pdffn = "Cirsium_biomes_DIVALIKE.pdf"
pdf(pdffn, width=12, height=16)

#######################################################
# Plot ancestral states - DIVALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE on Cirsium (biomes)"

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.05, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, 
                                plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, 
                                tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
                         label.offset=0.05, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, 
                         plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, 
                         tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

# Run Bayarelike
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object$master_table
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE         
BioGeoBEARS_run_object$use_optimx = "GenSA"  
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Cirsium_biomes_BAYAREALIKE.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}


pdffn = "Cirsium_biomes_BAYAREALIKE.pdf"
pdf(pdffn, width=12, height=16)

#######################################################
# Plot ancestral states - BAYAREALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Cirsium (biomes)"

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.05, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, 
                                plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, 
                                tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
                         label.offset=0.05, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

##################

