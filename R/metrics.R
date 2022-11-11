#' @title Performs the modeling of the observed 
#' data and returns the fit metrics of the 
#' studied model
#' @name metrics
#' @description A function that, based on the 
#' observed data, the independent variable 
#' (e.g. time in h) and the dependent
#'  variable (e.g. CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  production in g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^-1}}),
#'   performs the modeling of the fermentation 
#'   curve based on the chosen model (**5PL**, **Gompertz**, or **4PL**)
#'   and returns the model fit **metrics**. 
#'  
#'  As a result, the fit metrics for the 
#'  chosen model are returned in the form of 
#'  data.frame: **Correlation**, **R\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}**,**Residual sum of squares (RSS\ifelse{html}{\out{<sub>min</sub>}}{\eqn{_min}})** and **Residual standard error**.
#' 
#' @param data Data frame to be analyzed. The data frame must be in the following order:
#' \itemize{
#' \item **First**: All columns containing the independent 
#' variable (e.g. \emph{time in hours})
#' \item **Second**: All columns containing dependent variables
#'  (e.g. \emph{CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^-1}} 
#'  production})
#' \item **Header**: Columns must contain a 
#' header. If the treatment **ID** is in 
#' the header, this **ID** will be used 
#' to **identify** the metrics for each 
#' analyzed curve.} 
#' @param model Model to be adjusted. Argument for model:
#' \itemize{
#' \item **Model = 1**. 5PL Model (five-parameter logistic (5PL) model)
#' \item **Model = 2**. Gompertz Model
#' \item **Model = 3**. 4PL Model (four-parameter logistic (4PL) model)
#' } 
#' @param save.xls If TRUE, an xlsx file containing 
#' the metrics will be saved in the working directory. If FALSE, 
#' the xlsx file will not be saved.
#' @param dir.save Directory path where 
#' the xlsx file is to be saved.
#' @param xls.name File name. Must contain the 
#' format. For example, "Metrics.xlsx".
#' @param startA Starting estimate of the value of A for 5PL model.
#' @param startB Starting estimate of the value of B for 5PL model.
#' @param startC Starting estimate of the value of C for 5PL model.
#' @param startD Starting estimate of the value of D for 5PL model.
#' @param startG Starting estimate of the value of G for 5PL model.
#' @details 
#' Curve fitting from the observed data is 
#' performed by the nlsLM() function in 
#' the 'minpack.lm' package.
#'     
#' @return The **metrics** from the analyzed model 
#' are returned in a **data.frame**. In addition, 
#' a **"Metrics.xlsx" file** can be generated, 
#' containing the **model fit metrics** for each fermentation 
#' curve studied: **Correlation**; 
#' **R\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}**;
#' **Residual standard error**;
#' **Residual sum of squares (RSS\ifelse{html}{\out{<sub>min</sub>}}{\eqn{_min}})**.
#' 
#' @author Angelo Gava
#' @examples
#' 
#' #Creating a data.frame. 
#' #First, columns containing independent variable.
#' #Second, columns containing dependent variable.
#' #The data frame created presents two 
#' #fermentation curves for two yeasts with 
#' #different times and carbon dioxide production.
#' 
#' df <- data.frame('Time_Yeast_A' = seq(0,280, by=6.23),
#'                  'Time_Yeast_B' = seq(0,170, by=3.7777778),
#'                  'CO2_Production_Yeast_A' = c(0,0.97,4.04,9.62,13.44,17.50,
#'                                               24.03,27.46,33.75,36.40,40.80,
#'                                              44.24,48.01,50.85,54.85,57.51,
#'                                              61.73,65.43,66.50,72.41,75.47,
#'                                              77.22,78.49,79.26,80.31,81.04,
#'                                              81.89,82.28,82.56,83.13,83.62,
#'                                              84.11,84.47,85.02,85.31,85.61,
#'                                              86.05,86.27,85.29,86.81,86.94,
#'                                              87.13,87.33,87.45,87.85),
#'                  'CO2_Production_Yeast_B' = c(0,0.41,0.70,3.05,15.61,18.41,
#'                                               21.37,23.23,28.28,41.28,43.98,
#'                                               49.54,54.43,60.40,63.75,69.29,
#'                                               76.54,78.38,80.91,83.72,84.66,
#'                                               85.39,85.81,86.92,87.38,87.61,
#'                                               88.38,88.57,88.72,88.82,89.22,
#'                                               89.32,89.52,89.71,89.92,90.11,
#'                                               90.31,90.50,90.70,90.90,91.09,
#'                                               91.29,91.49,91.68,91.88))
#' 
#' #Using the metrics() function to find the 
#' #model fit metrics
#' 
#' metrics(data = df,
#' model = 1, 
#' startA = 0,
#' startB = 1.5,
#' startC = 500,
#' startD = 92, 
#' startG = 1500,
#' save.xls = FALSE) #5PL Model adopted
#' 
#' metrics(data = df,
#' model = 2,
#' startA = 92,
#' startB = 1.5,
#' startC = 0,
#' startD = NA, 
#' startG = NA, 
#' save.xls = FALSE) #Gompertz Model adopted
#' 
#' metrics(data = df,
#' model = 3,
#' startA = 0,
#' startB = 2.5,
#' startC = 10,
#' startD = 92, 
#' startG = NA, 
#' save.xls = FALSE) #4PL Model adopted 
#' 
#' #Saving an xlsx file. In this example, 
#' #we will use saving a temporary file in 
#' #the temporary file directories.
#' 
#'    
#' 
#' @import minpack.lm
#' @import openxlsx
#' @export 
#' 
metrics <- function(data, 
                    model,
                    save.xls = FALSE,
                    dir.save,
                    xls.name,
                    startA,
                    startB,
                    startC,
                    startD,
                    startG) {
  fit_5PL <- function(data, var_dep, var_indep){
    minpack.lm::nlsLM(var_dep ~ d + ((a-d)/((1+((var_indep/c)^b))^g)), 
                      data = data, start = list(a = startA, b = startB, c=startC,d=startD, g = startG),
                      control = minpack.lm::nls.lm.control(maxiter = 500))
  }
  fit_gompertz <- function(data, var_dep, var_indep) {
    minpack.lm::nlsLM(var_dep ~ a*exp(-exp(-c*var_indep+b)),
                      data = data,
                      start = list(a = startA, b= startB,c = startC),
                      control = minpack.lm::nls.lm.control(maxiter = 200))
  }
  fit_4PL <- function(data, var_dep, var_indep){
    minpack.lm::nlsLM(var_dep ~ d+(a-d)/(1+(var_indep/c)^b), 
                      data = data, start = list(a = startA, b = startB, c=startC,d=startD),
                      control = minpack.lm::nls.lm.control(maxiter = 200))
  }
  if (model == 1){
    metrics_5PL <- data.frame("R2" = rep(NA, ncol(data)/2), "Correlation" = rep(NA, ncol(data)/2), "Residual sum of squares (RSS)" = rep(NA, ncol(data)/2),"Residual standard error" = rep(NA, ncol(data)/2), check.names = FALSE)
    for (i in 1:(ncol(data)/2)) {
      model_5PL <- fit_5PL(data = data, 
                           var_indep = data[,i], 
                           var_dep = data[,i + ((ncol(data))/2)])
      r2_5PL <-stats::cor(data[i + (ncol(data)/2)],stats::predict(model_5PL)) * stats::cor(data[i + (ncol(data)/2)],stats::predict(model_5PL))
      Correlation_5PL <- stats::cor(data[i + (ncol(data)/2)],stats::predict(model_5PL))
      Residual_standard_error_5PL <- summary(model_5PL)[3]
      metrics_5PL[i, "R2"] <- r2_5PL
      metrics_5PL[i, "Correlation"] <- Correlation_5PL
      metrics_5PL[i,"Residual sum of squares (RSS)"]<- sum(stats::resid(model_5PL)^2)
      metrics_5PL[i, "Residual standard error"] <- Residual_standard_error_5PL
      row.names(metrics_5PL) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
    }
    print(metrics_5PL)
    if (save.xls == TRUE){
    openxlsx::write.xlsx(metrics_5PL, file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                         col.names = TRUE, row.names = TRUE, append = FALSE)}
  }
  else if (model == 2){
    metrics_gompertz <- data.frame("R2" = rep(NA, ncol(data)/2), "Correlation" = rep(NA, ncol(data)/2), "Residual sum of squares (RSS)" = rep(NA, ncol(data)/2),"Residual standard error" = rep(NA, ncol(data)/2), check.names = FALSE) 
    for (i in 1:(ncol(data)/2)) {
      model_Gompertz <- fit_gompertz(data = data, 
                                     var_indep = data[,i], 
                                     var_dep = data[,i + ((ncol(data))/2)])
      r2_gompertz <-stats::cor(data[i + (ncol(data)/2)],stats::predict(model_Gompertz)) * stats::cor(data[i + (ncol(data)/2)],stats::predict(model_Gompertz))
      Correlation_gompertz <- stats::cor(data[i + (ncol(data)/2)],stats::predict(model_Gompertz))
      Residual_standard_error_gompertz <- summary(model_Gompertz)[3]
      metrics_gompertz[i, "R2"] <- r2_gompertz
      metrics_gompertz[i, "Correlation"] <- Correlation_gompertz
      metrics_gompertz[i,"Residual sum of squares (RSS)"]<- sum(stats::resid(model_Gompertz)^2)
      metrics_gompertz[i, "Residual standard error"] <- Residual_standard_error_gompertz
      row.names(metrics_gompertz) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
    }
    print(metrics_gompertz)
    if (save.xls == TRUE){
    openxlsx::write.xlsx(metrics_gompertz, file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                         col.names = TRUE, row.names = TRUE, append = FALSE)}
    
  }
  else if (model == 3){
    metrics_4PL <- data.frame("R2" = rep(NA, ncol(data)/2), "Correlation" = rep(NA, ncol(data)/2), "Residual sum of squares (RSS)" = rep(NA, ncol(data)/2),"Residual standard error" = rep(NA, ncol(data)/2), check.names =FALSE)
    for (i in 1:(ncol(data)/2)) {
      model_4PL <- fit_4PL(data = data, 
                           var_indep = data[,i], 
                           var_dep = data[,i + ((ncol(data))/2)])
      r2_4PL <-stats::cor(data[i + (ncol(data)/2)],stats::predict(model_4PL)) * stats::cor(data[i + (ncol(data)/2)],stats::predict(model_4PL))
      Correlation_4PL <- stats::cor(data[i + (ncol(data)/2)],stats::predict(model_4PL))
      Residual_standard_error_4PL <- summary(model_4PL)[3]
      metrics_4PL[i, "R2"] <- r2_4PL
      metrics_4PL[i, "Correlation"] <- Correlation_4PL
      metrics_4PL[i,"Residual sum of squares (RSS)"]<- sum(stats::resid(model_4PL)^2)
      metrics_4PL[i, "Residual standard error"] <- Residual_standard_error_4PL
      row.names(metrics_4PL) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
    }
    print(metrics_4PL)
    if (save.xls == TRUE){
    openxlsx::write.xlsx(metrics_4PL, file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                         col.names = TRUE, row.names = TRUE, append = FALSE)}
  }
  else {print ("Model not found!!")}
}
