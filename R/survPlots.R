globalVariables(c("update", "predict", "x", "lower", "upper", "na.omit"))
#' Plots for survival models
#'
#' @import pROC
#' @import survival
#' @import ggplot2
#' @import prodlim
#' @import dcurves
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
#' @param netInterventions True or False to show net interventions for "Decision"
#' @param colour Colour of the plot line
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
survPlots= function(model, time, df,
                    eventVar = "Event", timeVar="Time",
                    plotType="ROC",
                    xlim=c(0,1), ylim=c(0,1),
                    netInterventions=F,
                    colour="black"){

  # Estimating the expected event status for right censored patients
  # Based on code from `pec`
  margForm <- update(model$formula, paste(".~1"))
  margFit <- prodlim(margForm, data = df)
  obs <- jackknife(margFit, times = time)

  df[timeVar] = time # Set the time here for the prediction
  pred = exp(-predict(model, newdata=df, type="expected")) # Survival probability

  # Discrimination

  if (plotType=="ROC"){
    rocObj = roc(as.factor(round(obs)), pred)
    dat.ci <- ci.se(rocObj, specificities = seq(0, 1, 0.01))
    dat.ci <- data.frame(x = as.numeric(rownames(dat.ci)), lower=dat.ci[,1], upper=dat.ci[,3])

    rocText = paste0("AUROC: ",round(rocObj$auc,3),
                     "\n DeLong 95% CI: (",
                     round(ci(rocObj)[1],3),
                     ", ",
                     round(ci(rocObj)[3],3),
                     ")"
    )

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

  # Calibration

  if (plotType=="Calibration"){
    cal_df = data.frame(obs, pred, event2 = abs(1-df[eventVar])) # Merge expected observations and probabilities
    cal_df = na.omit(cal_df)

    g = ggplot(cal_df, aes(x = pred, y = obs)) + # Plot the curve
      labs(x="Expected", y="Observed") + # Labels
      geom_abline(intercept = 0, slope=1, lwd=1, lty=2, col="red") + # Diagonal line
      scale_y_continuous(limits = ylim, oob = scales::oob_keep)+ # Set the y-axis
      scale_x_continuous(limits = xlim, oob = scales::oob_keep)+ # Set the y-axis
      theme_light()+  # Set the theme
      geom_smooth(method="loess")
  }

  if (plotType=="Decision"){
    df$pred = 1-pred
    decForm = update(model$formula, paste(".~pred"))
    g = dca(decForm,
            data = df,
            time = time,
            thresholds = seq(xlim[1], xlim[2], (xlim[2]-xlim[1])/100),
            label = list(pred = "Prediction Model")
    )
    if (netInterventions){
      g = net_intervention_avoided(g)
    }
    plot(g, smooth = TRUE)
  }

  if (plotType == "All"){
    rocObj = roc(as.factor(round(obs)), pred)
    dat.ci <- ci.se(rocObj, specificities = seq(0, 1, 0.01))
    dat.ci <- data.frame(x = as.numeric(rownames(dat.ci)), lower=dat.ci[,1], upper=dat.ci[,3])

    rocText = paste0("AUROC: ",round(rocObj$auc,3),
                     "\n DeLong 95% CI: (",
                     round(ci(rocObj)[1],3),
                     ", ",
                     round(ci(rocObj)[3],3),
                     ")"
    )

    g1 = ggroc(rocObj, "size"=1.2, colour=colour) +
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

    cal_df = data.frame(obs, pred, event2 = abs(1-df[eventVar])) # Merge expected observations and probabilities
    cal_df = na.omit(cal_df)

    g2 = ggplot(cal_df, aes(x = pred, y = obs)) + # Plot the curve
      labs(x="Expected", y="Observed") + # Labels
      geom_abline(intercept = 0, slope=1, lwd=1, lty=2, col="red") + # Diagonal line
      scale_y_continuous(limits = ylim, oob = scales::oob_keep)+ # Set the y-axis
      scale_x_continuous(limits = xlim, oob = scales::oob_keep)+ # Set the y-axis
      theme_light()+  # Set the theme
      geom_smooth(method="loess")


    df$pred = 1-pred
    decForm = update(model$formula, paste(".~pred"))
    g3 = dca(decForm,
            data = df,
            time = time,
            thresholds = seq(xlim[1], xlim[2], (xlim[2]-xlim[1])/100),
            label = list(pred = "Prediction Model")
    )
    if (netInterventions){
      g3 = net_intervention_avoided(g)
    }
    g3 = plot(g3, smooth = TRUE)

    g = ggarrange(
                  ggarrange(g1, g2,ncol=2,nrow=1),
                  g3,
              ncol=1, nrow=2)
  }

  g
}
