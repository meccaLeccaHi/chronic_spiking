# chronic_spiking

## Face-learning pipeline  
### Analysis pipeline for single-unit recordings of face-cells in macaque brain.  

These scripts perform unsupervised spike detection and sorting using using wavelets and super-paramagnetic clustering performed by *Wave_clus* (http://www.vis.caltech.edu/~rodri/Wave_clus/Wave_clus_home.htm **/** https://github.com/csn-le/wave_clus)

#### waveClus_shell.m
_Magnum opus_ script for running entire processing/analysis pipeline (mostly for reproducibility).

#### setFilePath.m
_Usage: [paths, monkeys] = setFilePath(proj)_
setting up appropriate paths, depending on the project specified and the machine being used.

## Processing
##

#### buildDatabase.m
_Usage: buildDatabase(header,paths,varargin)_
Script for selecting parsing method for particular project, and applying it to each file in the header file.

#### parseData.m
_*Relies on proprietary TDT drivers (http://www.tdt.com/nightly-software-updates.html)*_
_Usage: Builds data structure for all projects
dat = parseData(header,paths,varargin)_
spike-response data adhering to format.

#### buildSpikeDatabase.m
_Usage: snames = buildSpikeDatabase(header,paths)_
Script for aggregating datasets for single neurons over multiple days.

#### initFaceLearningHeader
_Usage: header = initFaceLearningHeader(expt)_
Function for initializing a blank header structure for faceLearning project.

#### decodeFaceLearning.m
_Usage: [type face step] = decodeFaceLearning(stim,h)_
Decoder function for faceLearning project look-up-table of stimuli.

#### initRadTanHeader
_Usage: header = initRadTanHeader(expt)_
Function for initializing a blank header structure for radTan project.

#### decodeRadTan.m
_Usage: [type face step] = decodeRadTan(stim,h)_
Decoder function for radTan project look-up-table of stimuli.


## Visualization
##

#### plot_spheres.m
Uses sphere function to represent neural responses projected in multi-dimensonal space.

####insp_dat_fine.m
_Usage: insp_dat_fine(paths,monkey,suffix)_
Cursory analysis script [faceLearning project]

####insp_dat_identTraj.m
_Usage: insp_dat_fine(paths,monkey,suffix)_
Performs analysis w.r.t. the identity levels of stimuli [faceLearning project]

####insp_dat_rtTraj.m
_Usage: insp_dat_rtTraj(paths,monkey,suffix)_
Performs analysis w.r.t. the identity levels of stimuli [radTan project]
Includes comparisons across both radial and tangential axes

### Population-level

####insp_dat_identTrajPop.m
_Usage: insp_dat_identTrajPop(paths,monkey,suffix)_
Itiratively reads files containing all neurons for each monkey and makes population plots [faceLearning project]

####insp_dat_rtPop.m
_Usage: Usage: insp_dat_identTraj(paths,monkey,suffix)_
Itiratively reads files containing all neurons for each monkey and makes population plots [faceLearning project]




