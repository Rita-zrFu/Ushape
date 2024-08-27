rm(list = ls())
library(openxlsx)
# summarize rda files from simulations into data.frame

# function to write results into excel
write_to_excel <- function(file_path, start_row, data, col_names) {
  # Load existing workbook or create a new one
  if (file.exists(file_path)) {
    wb <- loadWorkbook(file_path)
  } else {
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet1")
  }
  if(start_row ==1){
    writeData(wb, sheet = 1, x=as.data.frame(matrix(col_names, nrow=1)), startRow = start_row, startCol = 1, colNames = FALSE, rowNames = FALSE)
    start_row = start_row+1
    }
   writeData(wb, sheet = 1, x = data, startRow = start_row, startCol = 1, colNames = FALSE, rowNames = FALSE)
  saveWorkbook(wb, file_path, overwrite = TRUE)
}

mydir = ".../simulation_fin"
setwd(mydir)
  # column names of the sheet
mycolnames = c("sample size", "scenario", " ", "b0", "logb1", "a1", "a2", 
               "critical point", "critical region lower bound", "critical region upper bound",
               "C index", "survival probability")
file_path <- "example.xlsx"  # Path to your Excel file

### point estimate and inference
allfiles = list.files()
start_row = 1
for (i in 1:length(allfiles)) { ## for each scenario
  senario  = allfiles[i]
  subdir = paste0(mydir, "/", senario)
  setwd(subdir)
  print(subdir)
  point_est = NULL; 
  cindex_est = NULL; 
  cpoint_est = NULL; 
  cregin_est = NULL; 
  S_est = NULL; 
  S_true_est = NULL
  sd_par = NULL
  ci_par_mat = NULL; 
  ci_cpt_mat = NULL; 
  ci_crup_mat = NULL; 
  ci_crlo_mat = NULL
  conv = NULL
  CI.index = NULL
  CI.cp = NULL
  CI.crup = NULL
  CI.crlo = NULL
  bias_ls = NULL; 
  bias_cp_ls = NULL; 
  bias_cr_ls = NULL
  
  for (j in 1:length(list.files())) { ## for each replication
    load(list.files()[j])
    point_est = rbind(point_est, est_par)
    cindex_est  = append(cindex_est, est_cindex) 
    cpoint_est = append(cpoint_est, est_cpoint)
    cregin_est = rbind(cregin_est, est_cregin)
    sd_par = rbind(sd_par, sd_est)
    bias_ls = rbind(bias_ls,bias)
    bias_cp_ls = append(bias_cp_ls, bias_cp)
    bias_cr_ls = rbind(bias_cr_ls, bias_cr)
    ci_cpt_mat = rbind(ci_cpt_mat, ci_cpoint)
    ci_par_mat = rbind(ci_par_mat, ci_par)
    ci_crup_mat = rbind(ci_crup_mat, ci_cr_up)
    ci_crlo_mat = rbind(ci_crlo_mat, ci_cr_lo)
    conv = append(conv, converg)
    CI.index = rbind(CI.index, CI_index_par)
    CI.cp = append(CI.cp, CI_index_cp)
    CI.crup = append(CI.crup, CI_index_cr_up)
    CI.crlo = append(CI.crlo, CI_index_cr_lo)
    S_est = append(S_est, S)
  }

    setwd(mydir)
############### add in how you calculated the true value of C and S ######################
  ### cindex_true = 
  ### S_true = 
    vector_data <- matrix(c(apply(point_est, 2, mean) -theta_true, mean(bias_cp_ls), apply(bias_cr_ls, 2, mean), 
                            mean(cindex_est)-cindex_true, mean(S_est,na.rm = TRUE) - S_true, 
                     apply(point_est, 2, sd), sd(cpoint_est), apply(cregin_est, 2, sd), sd(cindex_est), sd(S_est, na.rm=TRUE),
                     apply(sd_par, 2, mean)[,1:4], apply(sd_par, 2, mean)[,7:9], "-", "-",  
                     apply(CI.index, 2, mean), mean(CI.cp), mean(CI.crlo), mean(CI.crup), "-", "-"), nrow=4, byrow=TRUE)  # Define the data to write

  write_to_excel(file_path, start_row, vector_data, mycolnames)  # Write data to the i-th row
  start_row = start_row + length(all.files)
}

