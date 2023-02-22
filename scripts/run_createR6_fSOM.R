# example script to render processing.Rmd directly
library(rmarkdown)

## SPADE clustering
### - either from Cytobank or done after pre-processing script

input.dir <- file.path("/projects","cytofProcess")

### set report parameters
report_params <- list(cytob.id = "17932", 
                      spade.id = "FlowSOM_16_40",
                      live.gate = "Stringent singlets",
                      fcs.dir = file.path(input.dir,"17932","FlowSOM_16_40"),
                      annotation = file.path(input.dir,"17932","curated","experiment_17932_annotations_clinical.rds"),
                      markers = file.path(input.dir,"17932","FlowSOM_16_40","markers.txt"),
                      gating = file.path(input.dir,"17932","FlowSOM_16_40","gating_v31.xml"), # copy the new gating there for consistency
                      latest.gating = file.path(input.dir,"17932","CytExp_17932_Gates_v31.xml"),
                      mst = file.path(input.dir,"17932","FlowSOM_16_40","mst.gml"),
                      metas = file.path(input.dir,"17932","FlowSOM_16_40","bubbles.rds"),
                      bubbles = file.path(input.dir,"17932","FlowSOM_16_40","bubbles_named.rds"),
                      info = file.path(input.dir,"17932","FlowSOM_16_40","info.yaml"),
                      overwrite = TRUE,
                      author = "Nevena Zivanovic")

render("reports/createR6.Rmd", 
       params = report_params, 
       output_dir = file.path("/projects/cytofProcess",
                              report_params$cytob.id),
       output_file = "cytof_fSOM.html",
       knit_root_dir = ".")
