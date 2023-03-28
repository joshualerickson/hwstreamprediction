# hwstreamprediction  

Modeling the distribution of headwater streams using topoclimatic indices, remote sensing and machine learning.  

**Abstract**  

Headwater streams are ecologically important components of mountain ecosystems. However, they are difficult to map and may not be accurately represented in existing spatial datasets. We used topographically resolved climatic water balance data and satellite indices retrieved from Google Earth Engine to model the occurrence (presence or absence) of headwater streams across Northwest Montana. We use a multi-scale feature selection (MSFS) procedure and boosted regression tree models/machine learning algorithms to identify variables associated with headwater stream occurrence.  Models that included climatic water balance deficit were more accurate than using only terrain indices and improved upon estimates of stream extent represented by the National Hydrography Dataset (NHDPlus). Including topoclimate captured the varying effect of upslope accumulated area across a strong moisture gradient. Multi-scale cross-validation (MSCV), coupled with a MSFS algorithm allowed us to find a parsimonious model that was not immediately evident using standard cross-validation procedures. NHDPlus flow lines had higher specificity but were less accurate overall, indicating these data may underpredict where streams occur. More accurate spatial model predictions of stream occurrence have potential for immediate application in land and water resource management, where significant field time can be spent identifying potential stream impacts prior to contracting and planning.    

## Methods  

Workflow Order below. Please access .tif files from [google drive](https://drive.google.com/drive/folders/100in8JlxDCbRLTPCzQ-p9D54SQxIP37d?usp=sharing) if wanting to explore.  

Some of the data needed will also be in `data` folder.

1. **PreProc.Rmd** - Pre-processing the data. Correlation, distributions, outliers, etc.

2. **ffs_script_clustering.R** - performing numerous spatial holdouts with ffs (MSFS).

3. **tuning_final.R** - Perform tuning on all model types (topoclimatic and topo-only) as well as 10 cluster scales (MSCV).

4. **final_results.R** - Statistic metrics used for final tables and graphing as well as figures generated in paper.
