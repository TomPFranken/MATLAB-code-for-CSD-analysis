This is MATLAB code to perform CSD (Current Source Density) analysis on data recorded with multielectrode probes such as Neuropixels in visual cortex. Current source density can be used to estimate the position of superficial, input and deep laminar compartments in the cortex (see Franken and Reynolds, eLife 2021; an example is shown in the included powerpoint). 
The code can be tested on the example file provided here: 10.6084/m9.figshare.30133237
The script plotCSD.m results in the CSD pattern.
The same data can be used to perform additional analysis to identify layers, i.e. the spike-phase-approach of Davis et al. (https://elifesciences.org/articles/84512)
For that analysis one needs to use the spike times. The last paragraph in plotCSD.m loads in the spike times from the data.



