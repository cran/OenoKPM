#' @title Performs the modeling of the observed 
#' data and returns the coefficients 
#' @name coeff
#' @description A function that, based on the 
#' observed data, the independent variable 
#' (e.g. time in h) and the dependent
#'  variable (e.g. CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  production in g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^-1}}),
#'   performs the modeling of the fermentation 
#'   curve based on the chosen model (**5PL**, **Gompertz**, or **4PL**). 
#'   
#'   As a result, the curve-fitted coefficients of 
#'   the chosen equation (model) are returned 
#'   in the form of data frame.
#'     
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
#' to **identify** the coefficients for each 
#' analyzed curve.} 
#' @param model Model to be adjusted. Argument for model:
#' \itemize{
#' \item **Model = 1**. 5PL Model (five-parameter logistic (5PL) model)
#' \item **Model = 2**. Gompertz Model
#' \item **Model = 3**. 4PL Model (four-parameter logistic (4PL) model)
#' } 
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
#'  
#' @return The coefficients of the analyzed model 
#'are returned in a data.frame.
#'
#' According to the model adopted, 
#' the following coefficients are presented:
#' 
#' 5PL: **a**, **b**, **c**, **d** and **g**
#' 
#' Gompertz: **a**, **b** and **c**
#' 
#' 4PL: **a**, **b**, **c** and **d**
#' 
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
#' #Using the coeff() function to return the fitted model coefficients
#' 
#' coeff(data = df,
#' model = 1,
#' startA = 0,
#' startB = 1.5,
#' startC = 500,
#' startD = 92, 
#' startG = 1500) #5PL Model adopted
#' 
#' coeff(data = df,
#' model = 2, 
#' startA = 92,
#' startB = 1.5,
#' startC = 0,
#' startD = NA, 
#' startG = NA) #Gompertz Model adopted
#' 
#' coeff(data = df,
#' model = 3, 
#' startA = 0,
#' startB = 2.5,
#' startC = 10,
#' startD = 92, 
#' startG = NA) #4PL Model adopted
#'
#'
#'@import minpack.lm
#'
#'
#'@export 
coeff <- function(data, model,startA,startB,startC,startD,startG){
  fit_5PL <- function(data, var_dep, var_indep){
    minpack.lm::nlsLM(var_dep ~ d + ((a-d)/((1+((var_indep/c)^b))^g)), 
                      data = data, start = list(a = startA, b = startB, c=startC,d=startD, g = startG),
                      control = minpack.lm::nls.lm.control(maxiter = 200))
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
  
  if(model == 1) { Coeff_5PL <- data.frame(a = rep(NA, ncol(data)/2),b = rep(NA, ncol(data)/2),c = rep(NA, ncol(data)/2),d = rep(NA, ncol(data)/2), g = rep(NA, ncol(data)/2), check.names = FALSE) 
  
  for (i in 1:(ncol(data)/2)) {
    model_5PL <- fit_5PL(data = data, 
                         var_indep = data[,i], 
                         var_dep = data[,i + ((ncol(data))/2)])
    Coeff_5PL[i,] <- as.data.frame(t(stats::coefficients(model_5PL)))
    row.names(Coeff_5PL) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
  }
  Coeff_5PL
  print(Coeff_5PL)} 
  else if (model == 2){
    Coeff_gompertz <- data.frame(a = rep(NA, ncol(data)/2),b = rep(NA, ncol(data)/2),c = rep(NA, ncol(data)/2), check.names = FALSE) 
    
    for (i in 1:(ncol(data)/2)) {
      model_gompertz <- fit_gompertz(data = data, 
                                     var_indep = data[,i], 
                                     var_dep = data[,i + ((ncol(data))/2)])
      Coeff_gompertz[i,] <- as.data.frame(t(stats::coefficients(model_gompertz)))
      row.names(Coeff_gompertz) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
    }
    Coeff_gompertz
    print(Coeff_gompertz)
  }
  else if (model == 3){
    Coeff_4PL <- data.frame(a = rep(NA, ncol(data)/2),b = rep(NA, ncol(data)/2),c = rep(NA, ncol(data)/2), d=rep(NA, ncol(data)/2), check.names = FALSE) 
    
    for (i in 1:(ncol(data)/2)) {
      model_4PL <- fit_4PL(data = data, 
                           var_indep = data[,i], 
                           var_dep = data[,i + ((ncol(data))/2)])
      Coeff_4PL[i,] <- as.data.frame(t(stats::coefficients(model_4PL)))
      row.names(Coeff_4PL) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
    }
    Coeff_4PL
    print(Coeff_4PL)
  }
  else {print ("Model not found!!")}
}
