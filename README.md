This repository contains the code used to generate the processed data and figures in "Increasing central serotonin with 5-HTP disrupts the inhibition of social gaze in non-human primates" (Weinberg-Wolf et al., 2021).

Reaction time, performance, and saccade-metric plots are generated in `hwwa_store_run_plot_approach_avoid_behav.m`, in the `script/` subdirectory. Gaze heatmaps are generated in `hwwa_target_onset_heatmap`, also in the `script/` subdirectory.

Before running these scripts, grab these dependencies and the source data:
* [shared_utils](https://github.com/nfagan/shared_utils)
* [categorical](https://github.com/nfagan/categorical)
* [dsp3](https://github.com/nfagan/dsp3)
* [data](...)

Once downloaded, unzip the source data and provide the path to the folder like so: 
```
conf = hwwa.config.load();
conf.PATHS.data_root = 'path/to/unzipped/folder';
hwwa.config.save(conf);
```