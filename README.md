# chronic_spiking
#### Matlab scripts for analysis of single-unit recordings of face-cells in the macaque brain.  

These scripts perform unsupervised spike detection and sorting using using wavelets and super-paramagnetic clustering performed using:  
**_Wave_clus_** *(https://github.com/csn-le/wave_clus)*

## Processing scripts

#### waveClus_shell.m
_Magnum opus_ script for running entire processing/analysis pipeline (mostly for reproducibility).

#### setFilePath.m
_Usage: [paths, monkeys] = setFilePath(proj)_  
Sets up appropriate paths, depending on the project specified and the local machine.

#### buildDatabase.m
_Usage: buildDatabase(header,paths,varargin)_  
Script for selecting parsing method for particular project, and applying it to each file in the header file.

#### parseData.m
_Usage: dat = parseData(header,paths,varargin)_
Parses data spiking-response for each stimulus.  
_*Relies on proprietary TDT drivers (http://www.tdt.com/nightly-software-updates.html)*_

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

## Visualization scripts

### _Face-learning project_

#### plot_spheres.m
Uses sphere function to represent neural responses projected in multi-dimensonal space.

#### insp_dat_fine.m
_Usage: insp_dat_fine(paths,monkey,suffix)_  
Cursory analysis script (plots rasters and spike-density functions).

#### insp_dat_identTraj.m
_Usage: insp_dat_fine(paths,monkey,suffix)_  
Performs analysis with respect to the identity levels of stimuli.

#### insp_dat_identTrajPop.m  

_Usage: insp_dat_identTrajPop(paths,monkey,suffix)_  
Iteratively reads files containing all neurons for each monkey and makes population plots.

### _Radial-tangential project_

#### insp_dat_rtTraj.m
_Usage: insp_dat_rtTraj(paths,monkey,suffix)_  
Performs analysis w.r.t. the identity levels of stimuli
Includes comparisons across both radial and tangential axes.

#### insp_dat_rtPop.m  
_Usage: Usage: insp_dat_identTraj(paths,monkey,suffix)_  
Iteratively reads files containing all neurons for each monkey and makes population plots.  

**Created by Dr Adam Jones  
Laboratory of Neurophysiology,  
National Institutes of Health,  
Bethesda, MD, USA** 


