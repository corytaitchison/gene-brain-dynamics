# SSP Project 2019: How the brain's activity dynamics are shaped by gene expression.

## Files

`final.m`: Main file used to correlate the gene data (neuroreceptor subsets) with the HCTSA time series features.

`final_cluster.m`: Plots a t-SNE cluster plot of gene expression for each subset.

### Analysis Files

- `create_acronyms.m`: Creates a table of acronyms for each gene subset
- `filter_region.m`: Only retain rows that match the region that given by {'Region'} (e.g. 'Isocortex')
- `gd_tsd_merge.m`: Merges the gene and time series data sets, ensuring rows (areas) are consistent between both
- `LoadGeneExpressionData.m`: Loads gene data from the file (not included)
- `remove_nans.m`: Removes nans from rows and columns above the given threshold

### Final Data

Contains the gene ontology reference lists for creating the gene subsets.

## Example Usage

```Matlab
addpath('data');
addpath('other_files');
addpath('analysis_files');

% Load time series data (not included)
load('100subj_avgHCTSA.mat', 'TS_DataMat', 'Operations');
load('100subj_joinedStructInfo.mat');

% Run program
[gf, gfp, features, operationsF] = final(TS_DataMat, joinedStructInfo, [2, 4, 6, 8, 9, 12, 14], 0.15, {'Isocortex'}, Operations, false);
```
