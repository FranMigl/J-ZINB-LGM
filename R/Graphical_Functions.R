# here we can find the graphical functions used to plot the results

library(ggplot2)

roc <- function(df){
  roc_plot <- ggplot(df, aes(x=FPR, y = TPR)) +
    geom_point(size = 1) +
    geom_point(data = df[which.max(df$J),], aes(x = FPR, y = TPR), colour = "blue", size = 4) +
    geom_line() + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = paste('ROC CURVE,', comment(df)), x = "False Positive Rate", y = "True Positive Rate", subtitle = paste('Best Youden Lambda =', round(df$Lambda[which.max(df$J)],3)), 
         caption = paste("N.Edges best λ = ", df$N.edges[which.max(df$J)], " / ", round(df$NZ[1]/2))) +
    geom_ribbon(aes(ymin = 0, ymax = TPR), fill = "darkorange", alpha = 0.5) +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_line(colour = "gray90"),
          axis.title = element_text(size=14, face = "bold"),
          plot.subtitle=element_text(color="darkslateblue")) +
    ylim(0,1)
  return(roc_plot)
  
}

prc <- function(df){
  prc_plot <- ggplot(df, aes(x = TPR, y = Precision)) +
    ylim(0,1) +
    geom_point(size = 1) +
    geom_point(data = df[which.max(df$J),], aes(x = TPR, y = Precision), colour = "brown1", size = 4) +
    geom_line() +
    geom_point(size = 1)+
    labs(title = paste("Precision Recall Curve,",comment(df)), x = "Recall (TPR)", y = "Precision", subtitle = paste('Best Youden Lambda =',round(df$Lambda[which.max(df$J)],3)), 
         caption = paste("N.Edges best λ = ", df$N.edges[which.max(df$J)], " / ", round(df$NZ[1]/2))) +
    geom_ribbon(aes(ymin = 0, ymax = Precision), fill = "darkslateblue", alpha = 0.5) +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_line(colour = "gray90"),
          axis.title = element_text(size=14, face = "bold"),
          plot.subtitle=element_text(color="darkred")) +
    ylim(0,1)
  return(prc_plot)
}

Accuracy <- function(df){
  acc_plot <- ggplot(df, aes(x = Lambda, y = Balanced_Accuracy)) +
    geom_point(size = 1) +
    geom_point(data = df[which.max(df$Balanced_Accuracy),], aes(x = df$Lambda[which.max(df$Balanced_Accuracy)], y = Balanced_Accuracy), colour = "brown1", size = 3) +
    geom_line() + 
    geom_point() +
    labs(title = paste("Balanced Accuracy with Lambda,", comment(df)), x = "Increasing Lambda", y = "Balanced Accuracy", subtitle = paste("Best Balanced Accuracy Lambda = ", round(df$Lambda[which.max(df$Balanced_Accuracy)],3)), 
         caption = paste("N.Edges best λ = ", df$N.edges[which.max(df$Balanced_Accuracy)], " / ", round(df$NZ[1]/2))) +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_line(colour = "gray90"),
          axis.title = element_text(size=14, face = "bold"),
          plot.subtitle=element_text(color="darkred")) +
    ylim(0,1)
  return(acc_plot)
}

F1_stat <- function(df){
  F1_plot <- ggplot(df, aes(x = Lambda, y = F1)) +
    geom_point(size = 1) +
    geom_point(data = df[which.max(df$F1),], aes(x = df$Lambda[which.max(df$F1)], y = F1), colour = "brown1", size = 3) +
    geom_line() + 
    geom_point() +
    labs(title = paste("F1 Statistic with Lambda,", comment(df)), x = "Increasing Lambda", y = "F1", subtitle = paste('Best F1 Lambda =',round(df$Lambda[which.max(df$F1)],3)), 
         caption = paste("N.Edges best λ = ", df$N.edges[which.max(df$F1)], " / ", round(df$NZ[1]/2))) +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_line(colour = "gray90"),
          axis.title = element_text(size=14, face = "bold"),
          plot.subtitle=element_text(color="darkred")) +
    ylim(0,1)
  return(F1_plot)
}