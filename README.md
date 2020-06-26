# hwstreamprediction
Predicting headwater streams in northwest Montana.

Workflow Order below. Please take any needed data (e.g., .shp, .csv, .tif files) from [google drive](https://drive.google.com/drive/folders/100in8JlxDCbRLTPCzQ-p9D54SQxIP37d?usp=sharing).

1. **PreProc** - Pre-processing the data. Correlation, distributions, outliers, etc.

2. **dataSplit** - Splitting the data into feature selection and training/evaluation. Also creating spatial dependance folds.

3. **customRF** - Adapt a custom random forest model.

4. **ffs** - Perform feature selection with the three dependance structures. Total of 3 models.

5. **tuning** - Perform tuning on all model types (topoclimatic and topo-only) as well as all dependance structures. Total of 6 models.

6. **nhd_add_resamps** - Bring it all together. Calculate resample statistics including NHDPlus.

7. **bootstrapping** - Create bootstrap distribution confidence intervals (standard error and percentile).
