#' @title Get the model's predicted values
#' @name pred
#' 
#' @description A function that, based on the 
#' observed data, the independent variable 
#' (e.g. time in h) and the dependent
#'  variable (e.g. CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  production in g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^{-1}}}),
#'   performs the modeling of the fermentation 
#'   curve based on the chosen model(s) (\strong{5PL}, \strong{Gompertz}, or/and \strong{4PL}). 
#'   
#'  From the analyzed data, this function will 
#'  provide the predicted data for each evaluated 
#'  fermentation curve.
#' @param data 
#' Data frame to be analyzed. 
#' The data frame must be in the 
#' following order:
#'\itemize{
#'\item \strong{First}: All columns containing the independent 
#' variable (e.g. \emph{time in hours})
#'\item \strong{Second}: All columns containing dependent variables
#'  (e.g. \emph{CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^{-1}}} 
#'  production})
#'\item \strong{Header}: Columns must contain a 
#' header. If the treatment \strong{ID} is in 
#' the header, this \strong{ID} will 
#' be used to name the graphics 
#' PDF files for each 
#' analyzed curve.
#' }
#' @param model Model or models to be adjusted:
#'\itemize{
#'\item \strong{Model = 1}. 5PL Model.
#'\item \strong{Model = 2}. Gompertz Model.
#'\item \strong{Model = 3}. 4PL Model.
#'}
#' @param startA Starting estimate of the value of A for model.
#' @param startB Starting estimate of the value of B for model.
#' @param startC Starting estimate of the value of C for model.
#' @param startD Starting estimate of the value of D for model.
#' @param startG Starting estimate of the value of G for model.
#' @param save.xls If TRUE, an xlsx file containing 
#' the predicted values of each curve will be saved 
#' in the working directory. If it is FALSE, the 
#' xlsx file will not be saved.
#' @param dir.save Directory path where 
#' the xlsx file is to be saved.
#' @param xls.name File name. Must contain the 
#' format. For example, "Predicted Values.xlsx".
#'@details 
#'Curve fitting from the observed data is 
#' performed by the nlsLM() function in 
#' the 'minpack.lm' package.
#' 
#' 
#' @return The predicted values of each analyzed 
#' curve will be returned in a data.frame. 
#' In addition, a file \strong{"Predicted Values.xlsx"} 
#' can be generated, containing the predicted values of 
#' each fermentation curve studied.
#'   
#'@author Angelo Gava
#'   
#'@examples
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
#' #Using the pred() function to find the 
#' #predicted valuesaccording to the adopted model.
#' 
#' pred(data = df,
#' model = 1, 
#' startA = 0,
#' startB = 1.5,
#' startC = 500,
#' startD = 92, 
#' startG = 1500,
#' save.xls = FALSE) #5PL Model adopted
#' 
#' pred(data = df,
#' model = 2,
#' startA = 92,
#' startB = 1.5,
#' startC = 0,
#' startD = NA, 
#' startG = NA, 
#' save.xls = FALSE) #Gompertz Model adopted
#' 
#' pred(data = df,
#' startA = 0,
#' startB = 2.5,
#' startC = 10,
#' startD = 92, 
#' startG = NA,
#' model = 3, 
#' save.xls = FALSE) #4PL Model adopted
#' 
#' #Saving an xlsx file. In this example, 
#' #we will use saving a temporary file in 
#' #the temporary file directories.
#' @import minpack.lm
#' @import openxlsx
#' @export
pred <- function(data, 
                 model,
                 startA,
                 startB,
                 startC,
                 startD,
                 startG,
                 save.xls = FALSE,
                 dir.save,
                 xls.name
) {
  
  fit_5PL <- function(data, var_dep, var_indep, treatment_name){
    predict_5PL <- stats::predict(minpack.lm::nlsLM(var_dep ~ d + ((a-d)/((1+((var_indep/c)^b))^g)), 
                                                    data = data, start = list(a = startA, b = startB, c=startC,d=startD, g = startG),
                                                    control = minpack.lm::nls.lm.control(maxiter = 500)))
  }
  fit_gompertz <- function(data, var_dep, var_indep, treatment_name){
    predict_gompertz <- stats::predict(minpack.lm::nlsLM(var_dep ~ a*exp(-exp(-c*var_indep+b)),
                                                         data = data,
                                                         start = list(a = startA, c = startC, b= startB),
                                                         control = minpack.lm::nls.lm.control(maxiter = 500)))
  }
  fit_4PL <- function(data, var_dep, var_indep, treatment_name){
    predict_4PL <- stats::predict(minpack.lm::nlsLM(var_dep ~ d+(a-d)/(1+(var_indep/c)^b), 
                                                    data = data, start = list(a = startA, b = startB, c=startC,d=startD),
                                                    control = minpack.lm::nls.lm.control(maxiter = 500)))
  }
  if(model == 1){
    pred_5PL <- data[,1:(ncol(data)/2)]
    for (i in 1:(ncol(data)/2)) {
      values_5PL <- fit_5PL(data = data, 
                            var_indep = data[,i], 
                            var_dep = data[,i + ((ncol(data))/2)])
      pred_5PL[(ncol(data)/2) + i] <- values_5PL
      colnames(pred_5PL)[(ncol(data)/2)+i] <- colnames(pred_5PL)[i]
    }
    if (save.xls == TRUE){
      openxlsx::write.xlsx(pred_5PL ,file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                           col.names = TRUE, row.names = TRUE, append = FALSE)}
    pred_5PL
    print(pred_5PL)}
  else if (model == 2){
    pred_gompertz <- data[,1:(ncol(data)/2)]
    for (i in 1:(ncol(data)/2)) {
      values_gompertz <- fit_gompertz(data = data, 
                                      var_indep = data[,i], 
                                      var_dep = data[,i + (((ncol(data))/2))])
      pred_gompertz[(ncol(data)/2) + i] <- values_gompertz
      colnames(pred_gompertz)[(ncol(data)/2)+i] <- colnames(pred_gompertz)[i]
    }
    if (save.xls == TRUE){
      openxlsx::write.xlsx(pred_gompertz ,file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                           col.names = TRUE, row.names = TRUE, append = FALSE)}
    pred_gompertz
    print(pred_gompertz)}
  else if (model == 3){
    pred_4PL <- data[,1:(ncol(data)/2)]
    for (i in 1:(ncol(data)/2)) {
      values_4PL <- fit_4PL(data = data, 
                            var_indep = data[,i], 
                            var_dep = data[,i + (((ncol(data))/2))])
      pred_4PL[(ncol(data)/2) + i] <- values_4PL
      colnames(pred_4PL)[(ncol(data)/2)+i] <- colnames(pred_4PL)[i]
    }
    if (save.xls == TRUE){
      openxlsx::write.xlsx(pred_4PL,file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                           col.names = TRUE, row.names = TRUE, append = FALSE)}
    pred_4PL
    print(pred_4PL)}
  else {print ("Model not found!!")}
}