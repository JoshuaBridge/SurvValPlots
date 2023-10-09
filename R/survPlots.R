globalVariables(c("update", "predict", "x", "lower", "upper", "na.omit"))
#' Plots for survival models
#'
#' @import pROC
#' @import survival
#' @import ggplot2
#' @import prodlim
#' @importFrom scales oob_keep
#' @importFrom ggpubr ggarrange
#'
#' @param model CoxPH model object
#' @param time Time point to show performance at
#' @param df Data contained in data frame
#' @param eventVar Name of the event variable in df
#' @param timeVar Name of the time variable in df
#' @param plotType One of "ROC", "Calibration", or "Decision"
#' @param xlim Limit of the x-axis for "Calibration" or "Decision"
#' @param ylim Limit of the y-axis for "Calibration"
#' @param colour Colour of the plot line
#' @param harm Harm of the model for the decision curve
#'
#' @return GGPlot object
#' @export
#'
#' @examples
#'
#' # Example using Rotterdam breast cancer data
#' library(survival)
#' mod <- coxph(Surv(rtime, recur) ~ age+meno+grade+size, data = rotterdam)
#' head(rotterdam)
#' survPlots(
#'   model = mod,
#'   time = 2000,
#'   df = rotterdam,
#'   eventVar = "recur",
#'   timeVar = "rtime",
#'   plotType="ROC")
#'
#' survPlots(
#'   model = mod,
#'   time = 2000,
#'   df = rotterdam,
#'   eventVar = "recur",
#'   timeVar = "rtime",
#'   plotType="Calibration")
#'
#' survPlots(
#'   model = mod,
#'   time = 2000,
#'   df = rotterdam,
#'   eventVar = "recur",
#'   timeVar = "rtime",
#'   plotType="Decision")
#'
#' survPlots(
#'   model = mod,
#'   time = 2000,
#'   df = rotterdam,
#'   eventVar = "recur",
#'   timeVar = "rtime",
#'   plotType="All")
survPlots= function(model,
                    time,
                    df,
                    eventVar = "Event",
                    timeVar="Time",
                    plotType="ROC",
                    xlim=c(0,1),
                    ylim=c(0,1),
                    colour="black",
                    harm=0){

  df[as.character(model$terms[[2]][[2]])] = df[timeVar]
  df[as.character(model$terms[[2]][[3]])] = df[eventVar]

  timeVar = as.character(model$terms[[2]][[2]])
  eventVar = as.character(model$terms[[2]][[3]])

  # Estimating the expected event status for right censored patients
  # Based on code from `pec`
  margForm <- update(model$formula, paste(".~1"))
  margFit <- prodlim(margForm, data = df)
  obs <- jackknife(margFit, times = time)

  df[timeVar] = time # Set the time here for the prediction
  pred = exp(-predict(model, newdata=df, type="expected")) # Survival probability

  # Discrimination

  if (plotType=="ROC"){
    g = roc.survvalplots(obs, pred, colour)
  }

  # Calibration

  if (plotType=="Calibration"){
    g = cal.survvalplots(obs, pred, xlim, ylim)
  }

  # Decision

  if (plotType=="Decision"){
    g = dec.survvalplots(obs, pred, model, df, time, xlim, harm)

  }

  # All

  if (plotType == "All"){
    g1 = roc.survvalplots(obs, pred, colour)
    g2 = cal.survvalplots(obs, pred, xlim, ylim)
    g3 = dec.survvalplots(obs, pred, model, df, time, xlim, harm)

    g = ggarrange(
                  ggarrange(g1, g2,ncol=2,nrow=1),
                  g3,
              ncol=1, nrow=2)
  }

  g
}


roc.survvalplots = function(obs, pred, colour){
  rocObj = roc(as.factor(round(obs)), pred) # The ROC object
  dat.ci <- ci.se(rocObj, specificities = seq(0, 1, 0.01)) # Confidence intervals
  dat.ci <- data.frame(x = as.numeric(rownames(dat.ci)), lower=dat.ci[,1], upper=dat.ci[,3])

  # This displays the AUROC on the plot
  rocText = paste0("AUROC: ",round(rocObj$auc,3),
                   "\n DeLong 95% CI: (",
                   round(ci(rocObj)[1],3),
                   ", ",
                   round(ci(rocObj)[3],3),
                   ")"
  )

  # GGplot object
  g = ggroc(rocObj, "size"=1.2, colour=colour) +
    labs(x="1-Specificity (FPR)", y="Sensitivity (TPR)", color="Model")+
    ylim(0,1) +
    xlim(1,0)+
    theme(legend.position="none")+
    theme_minimal()+
    geom_ribbon(
      data = dat.ci,
      aes(x = x, ymin = lower, ymax = upper),
      alpha = 0.2,
      inherit.aes = F) +
    geom_abline(intercept = 1, slope=1, lwd=1, lty=2, col="red")+
    annotate("text", x = 0.5, y = 0.1, label = rocText, size=3)+
    theme_light()
}


cal.survvalplots = function(obs, pred, xlim, ylim){
  cal_df = data.frame(obs, pred) # Merge expected observations and probabilities
  cal_df = na.omit(cal_df) # Remove NAs

  # GGplot object
  g = ggplot(cal_df, aes(x = pred, y = obs)) + # Plot the curve
    labs(x="Expected", y="Observed") + # Labels
    geom_abline(intercept = 0, slope=1, lwd=1, lty=2, col="red") + # Diagonal line
    scale_y_continuous(limits = ylim, oob = scales::oob_keep)+ # Set the y-axis
    scale_x_continuous(limits = xlim, oob = scales::oob_keep)+ # Set the y-axis
    theme_light()+  # Set the theme
    geom_smooth(method="loess")
}


dec.survvalplots = function(obs, pred, model, df, time, xlim, harm){
  pt = seq(0, 1, 0.02) # A range of probability thresholds
  obs2 = factor(round(obs), levels=c(0,1))

  # Model curve
  NB = pt # A vector of net benefits the same length as pt
  for (i in 1:length(pt)){ # For each threshold
    pred2 = as.numeric(pred>=pt[i]) # Set any predictions greater than the threshold to 1
    pred2 = factor(pred2, levels=c(0,1)) # Convert to predictor
    FP = table(pred2, obs2)[2]/length(pred2) # False positives divided by total
    TP = table(pred2, obs2)[4]/length(pred2) # True positives divided by total
    NB[i] = TP - FP*(pt[i]/(1-pt[i]))-harm # Net benefit
  }

  # Treat all curve
  NB1=pt # A vector of net benefits the same length as pt
  for (i in 1:length(pt)){ # For each threshold
    pred2 = rep(1, length(pred)) # Set all predictions to 1
    pred2 = factor(pred2, levels=c(0,1)) # For each threshold
    FP = table(pred2, obs2)[2]/length(pred2) # False positives divided by total
    TP = table(pred2, obs2)[4]/length(pred2) # True positives divided by total
    NB1[i] = TP - FP*(pt[i]/(1-pt[i]))-harm # Net benefit
  }

  # Treat all curve
  NB0=pt # A vector of net benefits the same length as pt
  for (i in 1:length(pt)){ # For each threshold
    pred2 = rep(0, length(pred)) # Set all predictions to 0
    pred2 = factor(pred2, levels=c(0,1)) # For each threshold
    FP = table(pred2, obs2)[2]/length(pred2) # False positives divided by total
    TP = table(pred2, obs2)[4]/length(pred2) # True positives divided by total
    NB0[i] = TP - FP*(pt[i]/(1-pt[i]))-harm # Net benefit
  }

  # Data frame of the net benefit
  nb_df = data.frame("Threshold" = pt,
                     "Model" = NB,
                     "Treat all" = NB1,
                     "Treat none" = NB0)

  # GGplot object
  g = ggplot(data=nb_df) + # Plot the curve
    geom_line(aes(x=pt, y=NB, color="Model"), size=1)+ # Net benefit of the model
    labs(x="Probability threshold", y="Net benefit") + # Labels
    geom_line(aes(x=pt, y=NB1, color="Treat all"), size=1)+ # Net benefit of treat all
    geom_line(aes(x=pt, y=NB0, color="Treat none"), size=1)+ # Net benefit of treat none
    ylim(-0.1,NA)+ # Show just below the x-axis
    xlim(xlim[1], xlim[2])+ # Set the x-axis limits
    theme_light()+  # Set the theme
    labs(color='') # Hide the legend title
}
