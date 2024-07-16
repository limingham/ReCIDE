library(Seurat)
library(fastSave)


source("/home/lmh/ReCIDE/benchmark_syq/main_run_all.R")
source("your_function_path1/ReCIDE_all_function.R")
source("your_function_path2/ReCIDE_all_function.R")

results_path <- "your_results_path"
if (dir.exists(results_path)) {
  cat("Folder exists:", results_path, "\n")
} else {
  cat("Folder does not exist\n")
  dir.create(results_path)
  cat("Folder created:", results_path, "\n")
}






seurat_all = readRDS.lbzip2('your_data_path/Seurat_data.rdsFS',n.cores = 50)

bulk_data <- readRDS("your_data_path/bulk_data.rds")
bulk_data[] <- lapply(bulk_data, function(x) as.integer(as.character(x)))


celltype1 = "cellType"
subject1 = "subject_id"

dir_ref_recide = paste0(results_path,'/ref_ReCIDE.rds')
dir_results_recide = paste0(results_path,'/results_ReCIDE.rds')
dir_results_bisque = paste0(results_path,'/results_bisque.rds')
dir_results_music = paste0(results_path,'/results_music.rds')
dir_results_bayes = paste0(results_path,'/results_bayes.rds')
dir_ref_CIBERSORT = paste0(results_path,'/ref_CIBERSORT.rds')
dir_results_CIBERSORT = paste0(results_path,'/results_CIBERSORT.txt')
dir_ref_DWLS = paste0(results_path,'/ref_DWLS.rds')
dir_results_DWLS = paste0(results_path,'/results_DWLS.rds')

######run_ReCIDE
func_run_ReCIDE(SC_ref = seurat_all,
                EXP_df = bulk_data,
                celltype = celltype1,
                subject = subject1,
                dir_ref = dir_ref_recide,
                dir_results = dir_results_recide)


######run_bisque
func_run_bisque(ref_seurat = seurat_all,
                bulkdata = bulk_data,
                celltype = celltype1,
                subject = subject1,
                dir_results = dir_results_bisque)
gc()

######run_music
func_run_music(ref_seurat = seurat_all,
                EXP_df = bulk_data,
                celltype = celltype1,
                subject = subject1,
                dir_results = dir_results_music)
gc()




######run_CIBERSORT
func_run_CIBERSORT(combined.data = seurat_all,
                  bulk.mtx = bulk_data,
                  celltype = celltype1,
                  subject = subject1,
                  dir_ref = dir_ref_CIBERSORT,
                  dir_results = dir_results_CIBERSORT)
gc()


######run_DWLS
func_run_DWLS(scdata_test = seurat_all,
              bulkdata = bulk_data,
              celltype = celltype1,
              dir_DWLS_ref = dir_ref_DWLS,
              dir_DWLS_results = dir_results_DWLS)
gc()


######run_bayes
func_run_bayes(ref_seurat = seurat_all,
               bulk.mtx = bulk_data,
               celltype = celltype1,
               dir_results = dir_results_bayes,
               key1 = NULL)

# rm(list = ls())
gc()
