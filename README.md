#### ReCIDE: Robust estimation of cell type proportions by integrating single-reference-based deconvolutions 

Author: Minghan Li

Environment: R

Correspondence:  If you have any questions, please ask on GitHub or contact mhli23@m.fudan.edu.cn. 

# Requirement:

Before using ReCIDE, please make sure the following packages are installed in your R environment. Sometimes, different versions of R packages may not cause significant issues, but it is beneficial to keep R package versions consistent.

| R package | versions | R package | versions |
| :-------: | :------: | :-------: | :------: |
|   DWLS    |  0.1.0   | ggbiplot  |   0.55   |
|  Seurat   |  4.4.0   | pbmcapply |  1.5.1   |
| fastSave  |  0.1.0   |  mclust   |  6.0.0   |
|   COSG    |  0.9.0   |  gmodels  | 2.18.1.1 |
|  stringr  |  1.5.0   | parallel  |  4.3.1   |
|   dplyr   |  1.1.3   |           |          |
|  ggplot2  |  3.4.3   |           |          |
| tidyverse |  2.0.0   |           |          |
|           |          |           |          |

# Usage

### 1. Load package

 Load the required dependencies: 

```R
library(DWLS)
library(Seurat)
library(fastSave)
library(COSG)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(pbmcapply)
library(mclust)
library(gmodels)
```

Download "ReCIDE_all_function.R" and load it: 

```
source('/your_path/ReCIDE_all_function.R')
```

### 2. Load DEMO data:

The DEMO data is from the cross datasets test of PFC, and all data can be obtained from https://drive.google.com/drive/folders/1TRz7Xj4Ddvt3NMTOHbpiTJBkTbDuaXzf?usp=drive_link. In addition to the reference single-cell RNA-Seq matrix, other data can also be obtained from  [ReCIDE/DEMO_data at main · limingham/ReCIDE (github.com)](https://github.com/limingham/ReCIDE/tree/main/DEMO_data) .

```
##Load single-cell reference and pseudo-bulk datasets
SC_ref = readRDS.lbzip2('/your_path/DEMO_data/REF_PFC.rdsFS',n.cores=50)
EXP_df = readRDS('/your_path/DEMO_data/EXP_PFC_df.rds')

##Obtaining reference cell type labels and subject labels
##Please note that the cell type labels should not contain any symbols other than "_".

celltype = SC_ref@meta.data[['cell type label']] ##In the demo data, it is major.celltype.
subject = SC_ref@meta.data[['subject labels']] ##In the demo data, it is subject.
```

(Optional) If computational resources are insufficient (<256GB), subjects in SC_ref can be sampled.

```
set.seed(123)
sample_subject=sample(unique(SC_ref@meta.data[["subject"]]),10)
SC_ref=subset(SC_ref,subject %in% sample_subject)
```

### 3. Construction of reference MGMs

```
ref_MGM.list = MGM_build(SC_ref = SC_ref,
                        EXP_df = EXP_df,
                        label_celltype = celltype,
                        label_subject = subject,
                        n_cores = 50)
```

SC_ref: single-cell expression matrix (row: gene, col: cell) ，Seurat or data.frame object

EXP_df: bulk expression matrix (row: gene, col: subject) ，data.frame object

label_celltype: Celltype labels corresponding to each cell in SC_ref 

label_subject: Subject labels corresponding to each cell in SC_ref 

n_cores: Number of parallel cores, default 5

### 4.  Deconvolution using each MGM  

```
ReCIDE_results = ReCIDE_deconvolution(Sig_list = ref_MGM.list,
                                      EXP_df = EXP_df,
                                      Method = 'DWLS',
                                      n_cores = 50,
                                      n_celltype = 0)
```

Sig_list:  list of reference MGMs output in step 3

EXP_df: bulk expression matrix (row: gene, col: subject) ，data.frame，Same as in step 3

Method: Deconvolution kernel, can be DWLS, CIBERSORT or FARDEEP, default is DWLS

n_cores: Number of parallel cores, default 5

n_celltype: Subjects in the reference set with a cell type number greater than n_celltype are used for deconvolution to prevent a situation where a cell type number is too low for a particular reference subject. The default is 0.

### 5. Results and Visualization

The final deconvolution results are stored in `ReCIDE_results[["results_final_df"]`, where Rowname is the cell type and Colname is the subject name.

```
ReCIDE_results[["results_final_df"]]
```





