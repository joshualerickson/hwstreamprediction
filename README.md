# hwstreamprediction
Predicting headwater streams in northwest Montana.

Workflow Order below. Please take any needed data (e.g., .shp, .csv, .tif files) from [google drive](https://drive.google.com/drive/folders/100in8JlxDCbRLTPCzQ-p9D54SQxIP37d?usp=sharing).

1. **PreProc** - Pre-processing the data. Correlation, distributions, outliers, etc.

2. **ffs_script_clustering** - performing numerous spatial holdouts with ffs.

3. **tuning** - Perform tuning on all model types (topoclimatic and topo-only) as well as 10 cluster dependance structures.

4. **nhd_add_resamps** - Bring it all together. Calculate resample statistics including NHDPlus.

5. **bootstrapping** - Create bootstrap distribution confidence intervals (standard error and percentile).