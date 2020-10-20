# TEPs-MEPs
## The correspondence between EMG and EEG measures of changes in cortical excitability following transcranial magnetic stimulation

Please read our paper, and if you use this code, please cite our paper:

M.Biabani, A.Fornito, J.Coxon, B.Fulcher, N.Rogasch. [The correspondence between EMG and EEG measures of changes in cortical excitability following transcranial magnetic stimulation].

This repository provides the code for reproducing a range of analyses on EEG and EMG signals recorded in response to transcranial magnetic stimulation (TMS) over motor cortex. The pipeline runs a set of analyses exploring how cortical responses to TMS measured by EMG (MEPs evoked by paired-pulse TMS) and EEG (TEPs elicited by single-pulse TMS) relate in shape and amplitude. 

### Dependencies
 [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php), [TESA](https://nigelrogasch.github.io/TESA/), and [FieldTrip](http://www.fieldtriptoolbox.org/) toolboxes.
 
### Data files
In order to reproduce data used in this project please download the data from [this figshare repository](https://doi.org/10.26180/5cff2a0fc38e9). 

Please create a directory in the root directory of this repository and call it **TEPs-MEPs**. Within this directory create two more directories: **Inputs** and **Outputs** and in each one again create two more directories: **Biphasic** and **Monophasic**.

Each file has an experiment ID (_Biphasic; _Monophasic) at the end of its name. Save the files into their the corresponding folders in Inputs directory and remove the ID of the file. (e.g., save TEPs_Biphasic.mat into Inputs > Biphasic folder and rename it to TEPs.mat). There is one file called **neighbour_template.mat** that shoulde be saved into both Biphasic and Monophasic folders. 

The data provided here are the results of pre-processing pipeline. The code for pre-processing is available at [this github repository](https://github.com/BMHLab/TEPs-PEPs). If you require the raw data please contact [Mana Biabani](mailto:mana.biabani@gmail.com).

### Data processing
After retrieving data, all of the processing analyses done for this project can be carried out by running processingPipeline_MB.m. This script is located in the Code folder. Before starting the analysis move both processingPipeline_MB.m and AddPaths.mat to the root repository (TEPs-MEPs folder). The script takes the clean data from a selected experiment (Biphasic or Monophasic) in Inputs, runs each analysis step and saves the outputs of each step to the corresponding folder in Outputs.

In case of any questions please contact [Mana Biabani](mailto:mana.biabani@gmail.com).
