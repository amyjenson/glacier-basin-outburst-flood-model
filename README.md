# glacier-basin-outburst flood model

This is an glacier outburst flood model combined with a glacier flow model to describe long term changes in outburst floods as the main trunk glacier retreats and glacier catchment geometry changes. 

To run the full model, first the glacier flow model must be run from 'run_glacier_flow_model.m' where the output for each year will be saved. After the glacier model is run and the all years of glacial retreat have been saved, the outburst flood model can be run from 'Running_NF_with_GC.m'. 

To change parameters for the glacier flow model, see 'glacier_flow_model.m' and 'run_glacier_flow_model.m'. For a detailed description of parameters for glacier flow model see Enderlin et al. (2013) and Carnahan et al. (2019). To change parameters for the glacier outburst flood model, see 'FullNyeFowlerForEnviro_updated.m' and 'Running_NF_with_GC.m'. For a detailed description of parameters in the outburst flood model see Kingslake (2013). To change basin parameters such as shape, storage capacity, remnant ice, etc., see 'LakeWithIce_with_GC.m'. For a detailed description of parameters for in the outburst flood model for basin geometry and for all remaining parameters, see the preprint of the manuscript, 'Long-period variability in ice-dammed glacier outburst floods due to evolving catchment geometry' found at https://doi.org/10.5194/tc-2021-141. 


These files have been created or adapted by Amy J. Jenson and Jason M. Amundson from Kingslake (2013), Enderlin et al. (2013), and Carnahan et al. (2019). 
