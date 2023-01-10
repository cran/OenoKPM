#' @title Calculates kinetic parameters as a 
#' function of model fit for \ifelse{html}{\out{CO<sub>2</sub>}}{\eqn{CO2}} production 
#' as a function of time
#' @name kp
#'@description A function that, based on the 
#' observed data, the independent variable 
#' (e.g. time in h) and the dependent
#'  variable (e.g. CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  production in g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^{-1}}}),
#'   performs the modeling of the fermentation 
#'   curve based on the chosen model (\strong{5PL}, \strong{Gompertz}, or \strong{4PL}). 
#'
#'   Next, the coefficients are used in 
#'   mathematical formulas to obtain the 
#'   following kinetic parameters:
#'   
#' \strong{t\ifelse{html}{\out{<sub>Lag</sub>}}{\eqn{_{Lag}}}} - Duration of the latency phase for CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} production;
#'
#' \strong{V\ifelse{html}{\out{<sub>max</sub>}}{\eqn{_{max}}}} - Maximum rate of production of CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}};
#'
#' \strong{t\ifelse{html}{\out{<sub>Vmax</sub>}}{\eqn{_{V_{max}}}}} - Moment in which maximum fermentation rate occurs;
#'
#' \strong{CO\ifelse{html}{\out{<sub>2Vmax</sub>}}{\eqn{_{2_{Vmax}}}}} - CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} Produced until Maximum fermentation rate occurs;
#'
#' \strong{Y\ifelse{html}{\out{<sub>max</sub>}}{\eqn{_{max}}}} - Maximum production of carbon dioxide (CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}});
#'
#' @param data Data frame to be analyzed. 
#' The data frame must be in the 
#' following order:
#' \itemize{
#' \item \strong{First}: All columns containing the independent 
#' variable (e.g. \emph{time in hours})
#' \item \strong{Second}: All columns containing dependent variables
#'  (e.g. \emph{CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^{-1}}} 
#'  production})
#' \item \strong{Header}: Columns must contain a 
#' header. If the treatment \strong{ID} is in 
#' the header, this \strong{ID} will be used 
#' to \strong{identify} the coefficients 
#' and kinetic parameters for each 
#' analyzed curve.}
#' 
#' @param model Model to be adjusted. Argument for model:
#' \itemize{
#' \item \strong{Model = 1}. 5PL Model (five-parameter logistic (5PL) model).
#' \item \strong{Model = 2}. Gompertz Model.
#' \item \strong{Model = 3}. 4PL Model (four-parameter logistic (4PL) model).
#' } 
#' @param save.xls If TRUE, an xlsx file containing 
#' the coefficients and kinetic parameters will 
#' be saved in the working directory. If FALSE, 
#' the xlsx file will not be saved.
#' @param dir.save Directory path where 
#' the xlsx file is to be saved.
#' @param xls.name File name. Must contain the 
#' format. For example, "Parameters.xlsx".
#' @param startA Starting estimate of the value of A for model.
#' @param startB Starting estimate of the value of B for model.
#' @param startC Starting estimate of the value of C for model.
#' @param startD Starting estimate of the value of D for model.
#' @param startG Starting estimate of the value of G for model.
#' @details 
#' 
#' Curve fitting from the observed data is 
#' performed by the nlsLM() function in 
#' the 'minpack.lm' package.
#' 
#' You can see our article for more details on 
#'the mathematical formulas used to obtain each 
#'kinetic parameter (Gava \emph{et al}., 2020). In addition, feel free to 
#'use it as a reference in your works.
#'
#' @return The analyzed model \strong{coefficients} 
#' and the calculated \strong{kinetic parameters} 
#' are returned in a data.frame. In addition, 
#' a \strong{"Parameters.xlsx" file} can be generated, 
#' containing the coefficients and kinetic 
#' parameters of each studied fermentation curve.
#'
#' @references Gava, A., Borsato, D., & Ficagna, E. 
#' (2020). Effect of mixture of fining agents on 
#' the fermentation kinetics of base wine for 
#' sparkling wine production: Use of methodology 
#' for modeling. \emph{LWT}, \emph{131}, 109660.
#' \doi{https://doi.org/10.1016/j.lwt.2020.109660}
#'
#'Zwietering, M. H., Jongenburger, I., Rombouts, F. M., 
#'& Van't Riet, K. J. A. E. M. (1990). 
#'Modeling of the bacterial growth curve. 
#'\emph{Applied and environmental microbiology}, \emph{56}(6),
#' 1875-1881. \doi{https://doi.org/10.1128/aem.56.6.1875-1881.1990}
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
#' #Using the kp() function to find the 
#' #coefficients and kinetic parameters 
#' #according to the adopted model.
#' 
#' kp(data = df,
#' model = 1, 
#' startA = 0,
#' startB = 1.5,
#' startC = 500,
#' startD = 92, 
#' startG = 1500,
#' save.xls = FALSE) #5PL Model adopted
#' 
#' kp(data = df,model = 2,
#' startA = 92,
#' startB = 1.5,
#' startC = 0,
#' startD = NA, 
#' startG = NA, 
#' save.xls = FALSE) #Gompertz Model adopted
#' 
#' kp(data = df,
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
#' 
#' 
#' @import minpack.lm
#' @import openxlsx
#' 
#' @export 
kp <- function(data, 
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
                      control = minpack.lm::nls.lm.control(maxiter = 500))
  }
  fit_4PL <- function(data, var_dep, var_indep){
    minpack.lm::nlsLM(var_dep ~ d+(a-d)/(1+(var_indep/c)^b), 
                      data = data, start = list(a = startA, b = startB, c=startC,d=startD),
                      control = minpack.lm::nls.lm.control(maxiter = 500))
  }
  if(model == 1){
    Coeff_5PL <- data.frame(a = rep(NA, ncol(data)/2),b = rep(NA, ncol(data)/2),c = rep(NA, ncol(data)/2),d = rep(NA, ncol(data)/2), g = rep(NA, ncol(data)/2)) 
    for (i in 1:(ncol(data)/2)) {
      model_5PL <- fit_5PL(data = data, 
                           var_indep = data[,i], 
                           var_dep = data[,i + ((ncol(data))/2)])
      Coeff_5PL[i,] <- as.data.frame(t(stats::coefficients(model_5PL)))
      row.names(Coeff_5PL) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
    }
    Parameters_5PL <- Coeff_5PL
    Parameters_5PL[] <- lapply(Parameters_5PL, function(x) as.numeric(as.character(x))) 
    Parameters_5PL$'tVmax (h)' <- with(Parameters_5PL, exp(log(((-1+b))/(b*g+1))/b)*c) 
    Parameters_5PL$'tLag (h)' <- with(Parameters_5PL,(d+(a-d)/((1+exp(log((-1+b)/(b*g+1))/b)^b)^g))/(a-d)/g/(exp(log((-1+b)/(b*g+1))/b)^b)/b*(1+exp(log((-1+b)/(b*g+1))/b)^b)^g*exp(log((-1+b)/(b*g+1))/b)*c*(1+exp(log((-1+b)/(b*g+1))/b)^b)+exp(log((-1+b)/(b*g+1))/b)*c) 
    Parameters_5PL$'Vmax (g/L/h)'<- with(Parameters_5PL,-(a-d)*g*exp(log((-1+b)/(b*g+1))/b)^b*b/((1+exp(log((-1+b)/(b*g+1))/b)^b)^g)/exp(log((-1+b)/(b*g+1))/b)/c/(1+exp(log((-1+b)/(b*g+1))/b)^b))
    Parameters_5PL$'CO2Vmax (g/L)'<-with(Parameters_5PL, d+(a-d)/((1+exp(log((-1+b)/(b*g+1))/b)^b)^g))
    Parameters_5PL$'Ymax (g/L)'<- Parameters_5PL$d 
    if (save.xls == TRUE){
    openxlsx::write.xlsx(Parameters_5PL ,file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                         col.names = TRUE, row.names = TRUE, append = FALSE)}
    Parameters_5PL
    print(Parameters_5PL)}
  else if(model == 2){
    Coeff_gompertz <- data.frame(a = rep(NA, ncol(data)/2),b = rep(NA, ncol(data)/2),c = rep(NA, ncol(data)/2), check.names = FALSE) 
    for (i in 1:(ncol(data)/2)) {
      model_gompertz <- fit_gompertz(data = data, 
                                     var_indep = data[,i], 
                                     var_dep = data[,i + ((ncol(data))/2)])
      Coeff_gompertz[i,] <- as.data.frame(t(stats::coefficients(model_gompertz)))
      row.names(Coeff_gompertz) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
    }
    Parameters_gompertz <- Coeff_gompertz
    Parameters_gompertz[] <- lapply(Parameters_gompertz, function(x) as.numeric(as.character(x)))
    Parameters_gompertz$'tVmax (h)' <- with(Parameters_gompertz, b/c)
    Parameters_gompertz$'tLag (h)' <- with(Parameters_gompertz, (b-1)/c)
    Parameters_gompertz$'Vmax (g/L/h)'<- with(Parameters_gompertz,a*c*exp(-1))
    Parameters_gompertz$'CO2Vmax (g/L)'<- with(Parameters_gompertz, a*exp(-1))
    Parameters_gompertz$'Ymax (g/L)'<- Parameters_gompertz$a
    if (save.xls == TRUE){
    openxlsx::write.xlsx(Parameters_gompertz, file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                         col.names = TRUE, row.names = TRUE, append = FALSE)}
    Parameters_gompertz
    print(Parameters_gompertz)
  } 
  else if(model == 3){
    Coeff_4PL <- data.frame(a = rep(NA, ncol(data)/2),b = rep(NA, ncol(data)/2),c = rep(NA, ncol(data)/2), d=rep(NA, ncol(data)/2), check.names = FALSE) 
    for (i in 1:(ncol(data)/2)) {
      model_4PL <- fit_4PL(data = data, 
                           var_indep = data[,i], 
                           var_dep = data[,i + ((ncol(data))/2)])
      Coeff_4PL[i,] <- as.data.frame(t(stats::coefficients(model_4PL)))
      row.names(Coeff_4PL) <- colnames(data[(((ncol(data))/2))+ 1:(ncol(data)/2)])
    }
    Parameters_4PL <- Coeff_4PL
    Parameters_4PL[] <- lapply(Parameters_4PL, function(x) as.numeric(as.character(x)))
    Parameters_4PL$'tVmax (h)' <- with(Parameters_4PL, exp(log((-1+b)/(b+1))/b)*c)
    Parameters_4PL$'tLag (h)' <- with(Parameters_4PL, c*(d*(exp(log((-1+b)/(b+1))/b)^b)^2+((b+1)*a-d*(-1+b))*exp(log((-1+b)/(b+1))/b)^b+a)*exp(log((-1+b)/(b+1))/b)/(a-d)/(exp(log((-1+b)/(b+1))/b)^b)/b)
    Parameters_4PL$'Vmax (g/L/h)'<- with(Parameters_4PL, -(a-d)/(1+exp(log((-1+b)/(b+1))/b)^b)^2*exp(log((-1+b)/(b+1))/b)^b*b/exp(log((-1+b)/(b+1))/b)/c)
    Parameters_4PL$'CO2Vmax (g/L)'<- with(Parameters_4PL, d+(a-d)/(1+exp(log((-1+b)/(b+1))/b)^b))
    Parameters_4PL$'Ymax (g/L)'<- Parameters_4PL$d
    if (save.xls == TRUE){
    openxlsx::write.xlsx(Parameters_4PL, file = paste(dir.save,xls.name, sep = ''), sheetName = "Sheet1", 
                         col.names = TRUE, row.names = TRUE, append = FALSE)}
    Parameters_4PL
    print(Parameters_4PL)
  } else {print ("Model not found!!")}
}
