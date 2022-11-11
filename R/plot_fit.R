#' @title Plot graphs with observed data 
#' and predicted data from models
#' @name plot_fit
#' 
#' @description A function that, based on the 
#' observed data, the independent variable 
#' (e.g. time in h) and the dependent
#'  variable (e.g. CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  production in g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^-1}}),
#'   performs the modeling of the fermentation 
#'   curve based on the chosen model(s) (**5PL**, **Gompertz**, or/and **4PL**). 
#'   
#' From the observed data and predicted data, 
#' whether from one or all of the available 
#' models, this function will plot a graph 
#' for each fermentation curve evaluated. 
#' The chart will have the following basic 
#' structure:
#'
#' **X axis**: \emph{fermentation time}
#'
#' **Y axis**: \emph{CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} production}
#'
#' **Observed data**: \emph{Scatterplot with dots}. Plot with geom_point function from ggplot2 package.
#'
#' **Predicted data**: \emph{Smoothed line}. Plot with the stat_smooth function from the ggplot2 package.  
#'  
#' @param data 
#' Data frame to be analyzed. 
#' The data frame must be in the 
#' following order:
#'\itemize{
#'\item **First**: All columns containing the independent 
#' variable (e.g. \emph{time in hours})
#'\item **Second**: All columns containing dependent variables
#'  (e.g. \emph{CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} 
#'  g L\ifelse{html}{\out{<sup>-1</sup>}}{\eqn{^-1}} 
#'  production})
#'\item **Header**: Columns must contain a 
#' header. If the treatment **ID** is in 
#' the header, this **ID** will 
#' be used to name the graphics 
#' PDF files for each 
#' analyzed curve.}
#' 
#' @param models Model or models to be adjusted:
#'\itemize{
#'\item **Models = 1**. Only the 5PL Model.
#'\item **Models = 2**. Only the Gompertz Model.
#'\item **Models = 3**. Only the 4PL Model.
#'\item **Models = 4**. 5PL and Gompertz Models.
#'\item **Models = 5**. 5PL and 4PL Models.
#'\item **Models = 6**. Gompertz and 4PL Models.
#'\item **Models = 7**. 5PL, Gompertz and 4PL Models.}
#' @param startA Starting estimate of the value of A for 5PL model.
#' @param startB Starting estimate of the value of B for 5PL model.
#' @param startC Starting estimate of the value of C for 5PL model.
#' @param startD Starting estimate of the value of D for 5PL model.
#' @param startG Starting estimate of the value of G for 5PL model.
#' @param col Plot color of observed data in points. For example, "black".
#' 
#' @param col1 Plot color of predicted data from **model 1** (\emph{5PL Model}). For example, "red".
#'  
#' @param col2 Plot color of predicted data from **model 2** (\emph{Gompertz Model}). For example, "blue".
#' 
#' @param col3 Plot color of predicted data from **model 3** (\emph{4PL Model}). For example, "green".
#'
#' @param axisX X Axis Title. Character vector (or expression).
#'
#' @param axisY Y Axis Title. Character vector (or expression).
#' 
#' @param breaksX 
#' One of:
#' 
#'\itemize{
#'\item Use ggplot2::waiver() for the default X-axis breaks calculated by the transform object.
#'\item Numerical vector of X-Axis scale positions. For example, seq(0,200,20).}
#'@param limitsX
#' 
#'   One of:
#'\itemize{
#'\item NULL to use the default X-Axis scaling range.
#'\item A numeric vector of length two, giving the limits of the X-axis scale. For example, c(0,200).}
#' @param breaksY
#' One of:
#'\itemize{
#'\item Use ggplot2::waiver() for the default Y-axis breaks calculated by the transform object.
#'\item A numerical vector of Y-Axis scale positions. For example, seq(0,100,10).}
#' 
#'@param limitsY
#'  
#'  One of:
#'\itemize{
#'\item NULL to use the default Y-Axis scaling range.
#'\item A numeric vector of length two, giving the limits of the Y-axis scale. For example, c(0,100).}
#' @param font Base font family
#' 
#' @param font.size Base font size, given in pts.
#' 
#' @param legend.position The position of the caption ("none", "left", "right", "bottom", "top", or a numeric vector of two elements (X,Y).
#' 
#' @param show.R2 If TRUE, plots the R\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}} of the plotted predicted models on the graph. If FALSE, do not plot the R\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}} of the plotted predicted models.
#'
#' @param save.PDF If TRUE, create a folder (directory)
#'  and save each graphic in PDF format. If FALSE, 
#'  it does not create a directory or save the 
#'  graphics in PDF.
#'  
#'@param dir.save Path of the directory in 
#'  which a new folder (directory) will be
#'   created for saving graphics in PDF format.
#'  
#' @param dir.name Folder name (directory name)
#'  to be created within the working directory 
#'  for saving PDF graphics. Character vector.
#'
#' @param width.PDF Width, in cm, of the graphic to be saved in a PDF file.
#' 
#' @param height.PDF Height, in cm, of the graphic to be saved in a PDF file.
#'@details 
#'Curve fitting from the observed data is 
#' performed by the nlsLM() function in 
#' the 'minpack.lm' package.
#' 
#' Graphs are 
#' plotted using the various functions in 
#' the 'ggplot2' package.
#' 
#' 
#' @return Elegant graphics plotted according 
#' to observed and predicted data. In addition, 
#' a folder (directory) can be created, 
#' in which the PDF graphics will be saved,
#'  if desired.
#' In this folder, 
#' the graph of each analyzed 
#' fermentation curve is saved in PDF format,
#'  with the dimensions stipulated 
#'  in the \emph{width.PDF} and \emph{height.PDF}
#'   arguments. The name of each PDF 
#'   file will be extracted from the 
#'   header of the dependent variable used 
#'   for the graph.See more in the examples.
#'   
#'@author Angelo Gava
#'   
#'@examples
#'   
#' #################Example 1#################
#' #Using only required arguments
#' 
#' #Creating a data.frame. 
#' #First, columns containing independent variable.
#' #Second, columns containing dependent variable.
#' #The data frame created presents two 
#' #fermentation curves for two yeasts with 
#' #different times and carbon dioxide production.
#'  
#' df <- data.frame('Time_Yeast_A' = seq(0,280, by=6.23),
#'  'Time_Yeast_B' = seq(0,170, by=3.7777778),
#'  'CO2_Production_Yeast_A' = c(0,0.97,4.04,9.62,13.44,17.50,
#'                               24.03,27.46,33.75,36.40,40.80,
#'                               44.24,48.01,50.85,54.85,57.51,
#'                               61.73,65.43,66.50,72.41,75.47,
#'                               77.22,78.49,79.26,80.31,81.04,
#'                               81.89,82.28,82.56,83.13,83.62,
#'                               84.11,84.47,85.02,85.31,85.61,
#'                               86.05,86.27,85.29,86.81,86.94,
#'                               87.13,87.33,87.45,87.85),
#'  'CO2_Production_Yeast_B' = c(0,0.41,0.70,3.05,15.61,18.41,
#'                               21.37,23.23,28.28,41.28,43.98,
#'                               49.54,54.43,60.40,63.75,69.29,
#'                               76.54,78.38,80.91,83.72,84.66,
#'                               85.39,85.81,86.92,87.38,87.61,
#'                               88.38,88.57,88.72,88.82,89.22,
#'                               89.32,89.52,89.71,89.92,90.11,
#'                               90.31,90.50,90.70,90.90,91.09,
#'                               91.29,91.49,91.68,91.88))
#' 
#' #Using the plot_fit function to 
#' #generate elegants graphs PDF files 
#' #containing both observed data and 
#' #predicted data.
#' 
#' 
#' #Graph plotted only with Model 5PL 
#' #fit (models = 1)
#' 
#' plot_fit(data = df,
#'          models = 1,
#'          startA = 0,
#'          startB = 1.5,
#'          startC = 500,
#'          startD = 92, 
#'          startG = 1500)
#'          
#' 
#' #Graph plotted with 5PL and Gompertz 
#' #model fits (models = 4)
#' 
#' plot_fit(data = df,
#'          models = 4,
#'          startA = 0,
#'          startB = 1.5,
#'          startC = 500,
#'          startD = 92, 
#'          startG = 1500)
#'          
#' 
#' 
#' #################Example 2#################
#' #Using the various function arguments to
#' #customize the graph.
#' 
#' #Creating a data.frame. 
#' #First, columns containing independent variable.
#' #Second, columns containing dependent variable.
#' #The data frame created presents two 
#' #fermentation curves for two yeasts with 
#' #different times and carbon dioxide production.
#' 
#' df <- data.frame('Time_Treatment_A' = seq(0,200, by=6.45),
#'                  'Time_Treatment_B' = seq(0,200, by=6.45),
#'                  'CO2_Production_Treatment_A' = c(0,0.47,0.78,3.23,19.15,22.86,
#'                                                   26.81,29.36,36.14,52.61,55.58,
#'                                                   61.38,66.25,71.83,74.8,78.88,
#'                                                   83.47,84.48,85.94,87.45,87.98,
#'                                                   88.42,88.68,89.40,89.72,89.87,
#'                                                   90.41,90.51,90.62,90.70,91.05,
#'                                                   91.185),
#'                  'CO2_Production_Treatment_B' = c(0,0.19,0.39,1.36,9.23,11.29,
#'                                                   13.58,15.06,19.34,30.92,33.28,
#'                                                   37.98,42.14,47.17,50.00,54.28,
#'                                                   60.92,62.80,65.54,69.74,71.52,
#'                                                   73.07,73.98,76.75,77.79,78.70,
#'                                                   80.65,81.48,82.07,82.47,84.04,
#'                                                   84.60))
#' 
#' #Using the plot_fit function to 
#' #generate elegants graphs PDF files 
#' #containing both observed data and 
#' #predicted data.
#' 
#' #Graph plotted only with Model 5PL 
#' #fit (models = 1)
#' #Do not show R^2
#' 
#' plot_fit(data = df,
#'          startA = 0,
#'          startB = 1.5,
#'          startC = 500,
#'          startD = 92, 
#'          startG = 1500, 
#'          models = 1, 
#'          col = "red", #Color of observed data (points)
#'          col1 = "blue", #Predicted data color from model 1 (line). Model = 1 <- 5PL Model
#'          axisX =  "Fermentation time (h)", #Title X-Axis
#'          axisY = "Carbon dioxide production (g/L)", #Title Y-Axis
#'          breaksX = seq(0,200,20), #X-Axis scale (positions). 0,20,40,60,80,...
#'          limitsX = c(0,200), #X-Axis Limits
#'          breaksY = seq(0,90,5),#Y-Axis scale (positions). 0,5,10,15,20,...
#'          limitsY = c(0,95), #Y-Axis Limits
#'          font = "serif",
#'          font.size = 12,
#'          legend.position = "right",
#'          show.R2 = FALSE) #Do not show R^2
#'        
#' 
#' 
#' 
#' #Graph plotted with 5PL and 4PL 
#' #model fits (models = 5)
#' #Show R^2
#' \dontrun{
#' plot_fit(data = df, 
#'          models = 5, 
#'          startA = 0,
#'          startB = 1.5,
#'          startC = 500,
#'          startD = 92, 
#'          startG = 1500,
#'          col = "#000000", #Color of observed data (points)
#'          col1 = "#FF0000", #Predicted data color from model 1 (line). Model = 1 <- 5PL Model
#'          col3 = "#0B6121",#Predicted data color from model 3 (line). Model = 3 <- 4PL Model
#'          axisX =  "Time (h)", #Title X-Axis
#'          axisY = "CO2 production (g/L)", #Title Y-Axis
#'          breaksX = seq(0,200,20), #X-Axis scale (positions). 0,20,40,60,80,...
#'          limitsX = c(0,200), #X-Axis Limits
#'          breaksY = seq(0,90,10),#Y-Axis scale (positions). 0,10,20,30,40,...
#'          limitsY = c(0,95), #Y-Axis Limits
#'          font = "serif",
#'          font.size = 14,
#'          legend.position = "bottom",
#'          show.R2 = TRUE) #Show R^2
#'  }        
#' 
#' #Graph plotted with 5PL, Gompertz and 4PL 
#' #model fits (models = 7)
#' #Do not show R^2
#' \dontrun{
#' plot_fit(data = df, 
#'          models = 7,
#'          startA = 0,
#'          startB = 1.5,
#'          startC = 500,
#'          startD = 92, 
#'          startG = 1500, 
#'          col = "#FF0000", #Color of observed data (points)
#'          col1 = "#FF00FF", #Predicted data color from model 1 (line). Model = 1 <- 5PL Model
#'          col2 = "#0101DF",#Predicted data color from model 2 (line). Model = 2 <- Gompertz Model
#'          col3 = "#088A08",#Predicted data color from model 3 (line). Model = 3 <- 4PL Model
#'          axisX =  "Time (h)", #Title X-Axis
#'          axisY = "Carbon dioxide production (g/L)", #Title Y-Axis
#'          breaksX = seq(0,200,20), #X-Axis scale (positions). 0,20,40,60,80,...
#'          limitsX = c(0,200), #X-Axis Limits
#'          breaksY = seq(0,90,10),#Y-Axis scale (positions). 0,10,20,30,40,...
#'          limitsY = c(0,95), #Y-Axis Limits
#'          font = "serif",
#'          font.size = 14,
#'          legend.position = "top",
#'          show.R2 = FALSE) #Do not show R^2
#'  }        
#' 
#' @import ggplot2
#' @import grDevices
#' @import minpack.lm
#' @export
plot_fit <- function(data, 
                      models,
                     startA,
                     startB,
                     startC,
                     startD,
                     startG,
                      col = "black",
                      col1 = "red",
                      col2 = "cornflowerblue",
                      col3 = "forestgreen",
                      axisX =  "Time (hours)",
                      axisY = expression(paste("CO"["2"]*" Production (g L"^{"-1"}*")")),
                      breaksX = seq(0,196,24),
                      limitsX = c(0,196),
                      breaksY = seq(0,100,10),
                      limitsY = c(0,100),
                      font = "serif",
                      font.size = 14,
                      legend.position = "top",
                      show.R2 = FALSE,
                      save.PDF = FALSE,
                      dir.save,
                      dir.name = "Graphics",
                      width.PDF = 15,
                      height.PDF = 12
) {
  graph_5PL <- function(data, var_dep, var_indep, treatment_name){
    model_5PL <- minpack.lm::nlsLM(var_dep ~ d + ((a-d)/((1+((var_indep/c)^b))^g)), 
                                   data = data, start = list(a = startA, b = startB, c=startC,d=startD, g = startG),
                                   control = minpack.lm::nls.lm.control(maxiter = 500))
    r2_5PL <-stats::cor(var_dep,stats::predict(model_5PL)) * stats::cor(var_dep,stats::predict(model_5PL))
    predicted_5PL <- stats::predict(model_5PL)
    models <- c("5PL Model" = col1)
    graph_5PL <- ggplot2::ggplot(data, ggplot2::aes_string(var_indep, var_dep))+
      ggplot2::geom_point(col= col, size=2.5, pch=21, bg = col)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_5PL, color = "5PL Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::theme_bw(base_family = font)+
      ggplot2::theme(legend.position = legend.position, legend.title = ggplot2::element_blank())+
      ggplot2::labs(color = "Legend")+
      ggplot2::scale_color_manual(values = models)+
      ggplot2::xlab(axisX) + 
      ggplot2::ylab(axisY)+
      ggplot2::scale_x_continuous(breaks = breaksX, limits = limitsX) +
      ggplot2::scale_y_continuous(breaks = breaksY, limits = limitsY)+
      ggplot2::theme(text=ggplot2::element_text(size = font.size))
    
    if(show.R2 == TRUE){graph_5PL <- graph_5PL + ggplot2::geom_text (x = 25, y = 90, parse = TRUE ,ggplot2::aes(label = "R2"))+
      ggplot2::geom_text (x=33, y=90, ggplot2::aes(label = " = "))+
      ggplot2::geom_text(x = 51, y=90, ggplot2::aes(label = round(r2_5PL, 7)))}
    if(save.PDF == TRUE){
    dir.create(paste(dir.save,dir.name, sep = ''))   
    ggplot2::ggsave(graph_5PL, filename = paste("Treatment", treatment_name, ".pdf"), width = width.PDF, height = height.PDF, units = 'cm', device = grDevices::cairo_pdf, path = paste(dir.save,dir.name, sep = ''))}
    print(graph_5PL)
  }
  
  graph_gompertz <- function(data, var_dep, var_indep, treatment_name){
    model_gompertz <- minpack.lm::nlsLM(var_dep ~ a*exp(-exp(-c*var_indep+b)),
                                        data = data,
                                        start = list(a = max(var_dep), c = min(var_dep), b= 1.5),
                                        control = minpack.lm::nls.lm.control(maxiter = 200))
    r2_gompertz <-stats::cor(var_dep,stats::predict(model_gompertz)) * stats::cor(var_dep,stats::predict(model_gompertz))
    predicted_gompertz <- stats::predict(model_gompertz)
    models <- c("Gompertz Model" = col2)
    graph_gompertz <- ggplot2::ggplot(data, ggplot2::aes_string(var_indep, var_dep))+
      ggplot2::geom_point(col= col, size=2.5, pch=21, bg = col)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_gompertz, color="Gompertz Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::theme_bw(base_family = font)+
      ggplot2::theme(legend.position = legend.position, legend.title = ggplot2::element_blank())+
      ggplot2::labs(color = "Legend")+
      ggplot2::scale_color_manual(values = models)+
      ggplot2::xlab(axisX) + 
      ggplot2::ylab(axisY)+
      ggplot2::scale_x_continuous(breaks = breaksX, limits = limitsX) +
      ggplot2::scale_y_continuous(breaks = breaksY, limits = limitsY)+
      ggplot2::theme(text=ggplot2::element_text(size = font.size))
    
    if(show.R2 == TRUE){graph_gompertz <- graph_gompertz + ggplot2::geom_text (x = 25, y = 90, parse = TRUE ,ggplot2::aes(label = paste("R^2")))+
      ggplot2::geom_text (x=33, y=90, ggplot2::aes(label = " = "))+
      ggplot2::geom_text(x = 51, y=90, ggplot2::aes(label = round(r2_gompertz, 7)))}
    if(save.PDF == TRUE){
    dir.create(paste(dir.save,dir.name, sep = ''))
    ggplot2::ggsave(graph_gompertz, filename = paste("Treatment", treatment_name, ".pdf"), width = width.PDF, height = height.PDF, units = 'cm', device = grDevices::cairo_pdf, path = paste(dir.save,dir.name, sep = ''))}
    print(graph_gompertz)
  }
  
  graph_4PL <- function(data, var_dep, var_indep, treatment_name){
    model_4PL <- minpack.lm::nlsLM(var_dep ~ d+(a-d)/(1+(var_indep/c)^b), 
                                   data = data, start = list(a = min(var_dep), b = 2.5, c=10,d=max(var_dep)),
                                   control = minpack.lm::nls.lm.control(maxiter = 200))
    r2_4PL <-stats::cor(var_dep,stats::predict(model_4PL)) * stats::cor(var_dep,stats::predict(model_4PL))
    predicted_4PL <- stats::predict(model_4PL)
    models <- c("4PL Model" = col3)
    graph_4PL <- ggplot2::ggplot(data, ggplot2::aes_string(var_indep, var_dep))+
      ggplot2::geom_point(col= col, size=2.5, pch=21, bg=col)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_4PL, color = "4PL Model"), se=FALSE, linetype = "longdash")+
      ggplot2::theme_bw(base_family = font)+
      ggplot2::theme(legend.position = legend.position, legend.title = ggplot2::element_blank())+
      ggplot2::labs(color = "Legend")+
      ggplot2::scale_color_manual(values = models)+
      ggplot2::xlab(axisX) + 
      ggplot2::ylab(axisY)+
      ggplot2::scale_x_continuous(breaks = breaksX, limits = limitsX) +
      ggplot2::scale_y_continuous(breaks = breaksY, limits = limitsY)+
      ggplot2::theme(text=ggplot2::element_text(size = font.size))
    
    if(show.R2 == TRUE){graph_4PL <- graph_4PL + ggplot2::geom_text (x = 25, y = 90, parse = TRUE ,ggplot2::aes(label = paste("R^2")))+
      ggplot2::geom_text (x=33, y=90, ggplot2::aes(label = " = "))+
      ggplot2::geom_text(x = 51, y=90, ggplot2::aes(label = round(r2_4PL, 7)))}
    if(save.PDF == TRUE){
    dir.create(paste(dir.save,dir.name, sep = ''))
    ggplot2::ggsave(graph_4PL, filename = paste("Treatment", treatment_name, ".pdf"), width = width.PDF, height = height.PDF, units = 'cm', device = grDevices::cairo_pdf, path = paste(dir.save,dir.name, sep = ''))}
    print(graph_4PL)
  }
  
  graph_5PLgompertz <- function(data, var_dep, var_indep, treatment_name){
    model_5PL <- minpack.lm::nlsLM(var_dep ~ d + ((a-d)/((1+((var_indep/c)^b))^g)), 
                                   data = data, start = list(a = startA, b = startB, c=startC,d=startD, g = startG),
                                   control = minpack.lm::nls.lm.control(maxiter = 200))
    model_gompertz <- minpack.lm::nlsLM(var_dep ~ a*exp(-exp(-c*var_indep+b)),
                                        data = data,
                                        start = list(a = max(var_dep), c = min(var_dep), b= 1.5),
                                        control = minpack.lm::nls.lm.control(maxiter = 200))
    r2_5PL <-stats::cor(var_dep,stats::predict(model_5PL)) * stats::cor(var_dep,stats::predict(model_5PL))
    predicted_5PL <- stats::predict(model_5PL)
    r2_gompertz <-stats::cor(var_dep,stats::predict(model_gompertz)) * stats::cor(var_dep,stats::predict(model_gompertz))
    predicted_gompertz <- stats::predict(model_gompertz)
    models <- c("5PL Model" = col1, "Gompertz Model" = col2)
    graph_5PLgompertz <- ggplot2::ggplot(data, ggplot2::aes_string(var_indep, var_dep))+
      ggplot2::geom_point(col= col, size=2.5, pch=21, bg=col)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_gompertz, color="Gompertz Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_5PL, color="5PL Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::theme_bw(base_family = font)+
      ggplot2::theme(legend.position = legend.position, legend.title = ggplot2::element_blank())+
      ggplot2::xlab(axisX) + 
      ggplot2::ylab(axisY)+
      ggplot2::labs(color = "Legend")+
      ggplot2::scale_color_manual(values = models)+
      ggplot2::scale_x_continuous(breaks = breaksX, limits = limitsX)+
      ggplot2::scale_y_continuous(breaks = breaksY, limits = limitsY)+
      ggplot2::theme(text=ggplot2::element_text(size = font.size))
    
    
    if(show.R2 == TRUE){graph_5PLgompertz <- graph_5PLgompertz + ggplot2::geom_text (x=6, y=90, ggplot2::aes(label = "Gompertz"), size = 4)+
      ggplot2::geom_text (x = 25, y = 90, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=90, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text(x = 51, y=90, ggplot2::aes(label = round(r2_gompertz, 7)), size = 4)+
      ggplot2::geom_text (x = 25, y = 85, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=85, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text (x=10, y=85, ggplot2::aes(label = "5PL"), size = 4)+
      ggplot2::geom_text(x = 51, y=85, ggplot2::aes(label = round(r2_5PL, 7)), size = 4)}
    if(save.PDF == TRUE){
    dir.create(paste(dir.save,dir.name, sep = ''))
    ggplot2::ggsave(graph_5PLgompertz, filename = paste("Treatment", treatment_name, ".pdf"), width = width.PDF, height = height.PDF, units = 'cm', device = grDevices::cairo_pdf, path = paste(dir.save,dir.name, sep = ''))}
    print(graph_5PLgompertz)
  }
  
  graph_5PL4PL <- function(data, var_dep, var_indep, treatment_name){
    model_5PL <- minpack.lm::nlsLM(var_dep ~ d + ((a-d)/((1+((var_indep/c)^b))^g)), 
                                   data = data, start = list(a = startA, b = startB, c=startC,d=startD, g = startG),
                                   control = minpack.lm::nls.lm.control(maxiter = 500))
    model_4PL <- minpack.lm::nlsLM(var_dep ~ d+(a-d)/(1+(var_indep/c)^b), 
                                   data = data, start = list(a = min(var_dep), b = 2.5, c=10,d=max(var_dep)),
                                   control = minpack.lm::nls.lm.control(maxiter = 200))
    r2_5PL <-stats::cor(var_dep,stats::predict(model_5PL)) * stats::cor(var_dep,stats::predict(model_5PL))
    predicted_5PL <- stats::predict(model_5PL)
    r2_4PL <-stats::cor(var_dep,stats::predict(model_4PL)) * stats::cor(var_dep,stats::predict(model_4PL))
    predicted_4PL <- stats::predict(model_4PL)
    models <- c("5PL Model" = col1, "4PL Model" = col3)
    graph_5PL4PL <- ggplot2::ggplot(data, ggplot2::aes_string(var_indep, var_dep))+
      ggplot2::geom_point(col= col, size=2.5, pch=21, bg=col)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_4PL, color="4PL Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_5PL, color="5PL Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::theme_bw(base_family = font)+
      ggplot2::theme(legend.position = legend.position, legend.title = ggplot2::element_blank())+
      ggplot2::xlab(axisX) + 
      ggplot2::ylab(axisY)+
      ggplot2::labs(color = "Legend")+
      ggplot2::scale_color_manual(values = models)+
      ggplot2::scale_x_continuous(breaks = breaksX, limits = limitsX) +
      ggplot2::scale_y_continuous(breaks = breaksY, limits = limitsY)+
      ggplot2::theme(text=ggplot2::element_text(size = font.size))
    
    if(show.R2 == TRUE){graph_5PL4PL <- graph_5PL4PL + ggplot2::geom_text (x=15, y=90, ggplot2::aes(label = "4PL"), size = 4)+
      ggplot2::geom_text (x = 25, y = 90, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=90, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text(x = 51, y=90, ggplot2::aes(label = round(r2_4PL, 7)), size = 4)+
      ggplot2::geom_text (x = 25, y = 85, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=85, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text (x=15, y=85, ggplot2::aes(label = "5PL"), size = 4)+
      ggplot2::geom_text(x = 51, y=85, ggplot2::aes(label = round(r2_5PL, 7)), size = 4)
    }
    if(save.PDF == TRUE){
    dir.create(paste(dir.save,dir.name, sep = ''))
    ggplot2::ggsave(graph_5PL4PL, filename = paste("Treatment", treatment_name, ".pdf"), width = width.PDF, height = height.PDF, units = 'cm', device = grDevices::cairo_pdf, path = paste(dir.save,dir.name, sep = ''))}
    print(graph_5PL4PL)
  }
  
  graph_Gompertz4PL <- function(data, var_dep, var_indep, treatment_name){
    model_gompertz <- minpack.lm::nlsLM(var_dep ~ a*exp(-exp(-c*var_indep+b)),
                                        data = data,
                                        start = list(a = max(var_dep), c = min(var_dep), b= 1.5),
                                        control = minpack.lm::nls.lm.control(maxiter = 200))
    model_4PL <- minpack.lm::nlsLM(var_dep ~ d+(a-d)/(1+(var_indep/c)^b), 
                                   data = data, start = list(a = min(var_dep), b = 2.5, c=10,d=max(var_dep)),
                                   control = minpack.lm::nls.lm.control(maxiter = 200))
    r2_gompertz <-stats::cor(var_dep,stats::predict(model_gompertz)) * stats::cor(var_dep,stats::predict(model_gompertz))
    predicted_gompertz <- stats::predict(model_gompertz)
    r2_4PL <-stats::cor(var_dep,stats::predict(model_4PL)) * stats::cor(var_dep,stats::predict(model_4PL))
    predicted_4PL <- stats::predict(model_4PL)
    models <- c("Gompertz Model" = col2, "4PL Model" = col3)
    graph_Gompertz4PL <- ggplot2::ggplot(data, ggplot2::aes_string(var_indep, var_dep))+
      ggplot2::geom_point(col= col, size=2.5, pch=21, bg=col)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_4PL, color="4PL Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_gompertz, color="Gompertz Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::theme_bw(base_family = font)+
      ggplot2::theme(legend.position = legend.position, legend.title = ggplot2::element_blank())+
      ggplot2::xlab(axisX) + 
      ggplot2::ylab(axisY)+
      ggplot2::labs(color = "Legend")+
      ggplot2::scale_color_manual(values = models)+
      ggplot2::scale_x_continuous(breaks = breaksX, limits = limitsX) +
      ggplot2::scale_y_continuous(breaks = breaksY, limits = limitsY)+
      ggplot2::theme(text=ggplot2::element_text(size = font.size))
    
    
    if(show.R2 == TRUE){graph_Gompertz4PL <- graph_Gompertz4PL + ggplot2::geom_text (x=15, y=90, ggplot2::aes(label = "4PL"), size = 4)+
      ggplot2::geom_text (x = 25, y = 90, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=90, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text(x = 51, y=90, ggplot2::aes(label = round(r2_4PL, 7)), size = 4)+
      ggplot2::geom_text (x = 25, y = 85, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=85, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text (x=5, y=85, ggplot2::aes(label = "Gompertz"), size = 4)+
      ggplot2::geom_text(x = 51, y=85, ggplot2::aes(label = round(r2_gompertz, 7)), size = 4)}
    if(save.PDF == TRUE){
    dir.create(paste(dir.save,dir.name, sep = ''))
    ggplot2::ggsave(graph_Gompertz4PL, filename = paste("Treatment", treatment_name, ".pdf"), width = width.PDF, height = height.PDF, units = 'cm', device = grDevices::cairo_pdf, path = paste(dir.save,dir.name, sep = ''))}
    print(graph_Gompertz4PL)
  }
  
  graph_5PLGompertz4PL <- function(data, var_dep, var_indep, treatment_name){
    model_5PL <- minpack.lm::nlsLM(var_dep ~ d + ((a-d)/((1+((var_indep/c)^b))^g)), 
                                   data = data, start = list(a = startA, b = startB, c=startC,d=startD, g = startG),
                                   control = minpack.lm::nls.lm.control(maxiter = 500))
    model_gompertz <- minpack.lm::nlsLM(var_dep ~ a*exp(-exp(-c*var_indep+b)),
                                        data = data,
                                        start = list(a = max(var_dep), c = min(var_dep), b= 1.5),
                                        control = minpack.lm::nls.lm.control(maxiter = 200))
    model_4PL <- minpack.lm::nlsLM(var_dep ~ d+(a-d)/(1+(var_indep/c)^b), 
                                   data = data, start = list(a = min(var_dep), b = 2.5, c=10,d=max(var_dep)),
                                   control = minpack.lm::nls.lm.control(maxiter = 200))
    r2_5PL <-stats::cor(var_dep,stats::predict(model_5PL)) * stats::cor(var_dep,stats::predict(model_5PL))
    predicted_5PL <- stats::predict(model_5PL)
    r2_gompertz <-stats::cor(var_dep,stats::predict(model_gompertz)) * stats::cor(var_dep,stats::predict(model_gompertz))
    predicted_gompertz <- stats::predict(model_gompertz)
    r2_4PL <-stats::cor(var_dep,stats::predict(model_4PL)) * stats::cor(var_dep,stats::predict(model_4PL))
    predicted_4PL <- stats::predict(model_4PL)
    models <- c("5PL Model" = col1, "Gompertz Model" = col2, "4PL Model" = col3)
    graph_5PLGompertz4PL <- ggplot2::ggplot(data, ggplot2::aes_string(var_indep, var_dep))+
      ggplot2::geom_point(col= col, size=2.5, pch=21, bg=col)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_5PL, color="5PL Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_gompertz, color="Gompertz Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::stat_smooth(ggplot2::aes(var_indep, predicted_4PL, color="4PL Model"), se=FALSE, linetype = "longdash", show.legend = TRUE)+
      ggplot2::theme_bw(base_family = font)+
      ggplot2::theme(legend.position = legend.position, legend.title = ggplot2::element_blank())+
      ggplot2::xlab(axisX) + 
      ggplot2::ylab(axisY)+
      ggplot2::labs(color = "Legend")+
      ggplot2::scale_color_manual(values = models)+
      ggplot2::scale_x_continuous(breaks = breaksX, limits = limitsX) +
      ggplot2::scale_y_continuous(breaks = breaksY, limits = limitsY)+
      ggplot2::theme(text=ggplot2::element_text(size = font.size))
    
    
    if(show.R2 == TRUE){graph_5PLGompertz4PL <- graph_5PLGompertz4PL + ggplot2::geom_text (x=15, y=90, ggplot2::aes(label = "5PL"), size = 4)+
      ggplot2::geom_text (x = 25, y = 90, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=90, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text(x = 51, y=90, ggplot2::aes(label = round(r2_5PL, 7)), size = 4)+
      ggplot2::geom_text (x = 25, y = 85, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=85, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text (x=5, y=85, ggplot2::aes(label = "Gompertz"), size = 4)+
      ggplot2::geom_text(x = 51, y=85, ggplot2::aes(label = round(r2_gompertz, 7)), size = 4)+
      ggplot2::geom_text (x=15, y=80, ggplot2::aes(label = "4PL"), size = 4)+
      ggplot2::geom_text (x = 25, y = 80, parse = TRUE ,ggplot2::aes(label = paste("R^2")), size = 4)+
      ggplot2::geom_text (x=33, y=80, ggplot2::aes(label = " = "), size = 4)+
      ggplot2::geom_text(x = 51, y=80, ggplot2::aes(label = round(r2_4PL, 7)), size = 4)}
    if(save.PDF == TRUE){
    dir.create(paste(dir.save,dir.name, sep = ''))
    ggplot2::ggsave(graph_5PLGompertz4PL, filename = paste("Treatment", treatment_name, ".pdf"), width = width.PDF, height = height.PDF, units = 'cm', device = grDevices::cairo_pdf, path = paste(dir.save,dir.name, sep = ''))}
    print(graph_5PLGompertz4PL)
    
  }
  
  if(models == 1){
    for (i in 1:(ncol(data)/2)) {graph_5PL(data = data, 
                                           var_dep = data[,i + (ncol(data)/2)],
                                           var_indep = data[,i],
                                           treatment_name = colnames(data)[i + (ncol(data)/2)])}
  }
  else if (models == 2){for (i in 1:(ncol(data)/2)) {graph_gompertz(data = data, 
                                                                   var_dep = data[,i + (ncol(data)/2)],
                                                                   var_indep = data[,i],
                                                                   treatment_name = colnames(data)[i + (ncol(data)/2)])}
    
  }
  else if (models == 3){for (i in 1:(ncol(data)/2)) {graph_4PL(data = data, 
                                                              var_dep = data[,i + (ncol(data)/2)],
                                                              var_indep = data[,i],
                                                              treatment_name = colnames(data)[i + (ncol(data)/2)])}}
  else if (models == 4){for (i in 1:(ncol(data)/2)) {graph_5PLgompertz(data = data, 
                                                                      var_dep = data[,i + (ncol(data)/2)],
                                                                      var_indep = data[,i],
                                                                      treatment_name = colnames(data)[i + (ncol(data)/2)])}}
  else if (models == 5){
    for (i in 1:(ncol(data)/2)) {graph_5PL4PL(data = data, 
                                              var_dep = data[,i + (ncol(data)/2)],
                                              var_indep = data[,i],
                                              treatment_name = colnames(data)[i + (ncol(data)/2)])}}
  else if (models == 6){
    for (i in 1:(ncol(data)/2)) {graph_Gompertz4PL(data = data, 
                                                   var_dep = data[,i + (ncol(data)/2)],
                                                   var_indep = data[,i],
                                                   treatment_name = colnames(data)[i + (ncol(data)/2)])}}
  else if (models == 7){
    for (i in 1:(ncol(data)/2)) {graph_5PLGompertz4PL(data = data, 
                                                      var_dep = data[,i + (ncol(data)/2)],
                                                      var_indep = data[,i],
                                                      treatment_name = colnames(data)[i + (ncol(data)/2)])}}
  else {print("Models not found!!")}}