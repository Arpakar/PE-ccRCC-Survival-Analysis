library(survival)
library(survminer)
kmsurvo<-Surv(combined$Rectime,combined$Censor.Status)
sfit <- survfit(Surv(combined$Rectime,combined$Censor.Status)~RAP2B, data=combined)
ggsurvplot(sfit, pval=TRUE, 
           legend.labs=c("Normal", "Altered"), legend=c(.75,.75),  
           title="RAP2B")
ggsurvplot(sfit, legend = c(0.2, 0.2))
plot(sfit)

ggsurv <- ggsurvplot(
  sfit,
  size = 1.25,
  palette = c("forestgreen", "firebrick"),
  pval = 0.043125,
  pval.size = 5,
  legend.title = "",
  legend.labs = c("Normal", "Altered"),
  xlab = "Time in days",
  font.x = c(16, "italic", "black"),
  font.y = c(16, "italic", "black"),
  legend = c(.75, .75),
  title = "RAP2B",
  main = "Survival curve",
  font.main = c(20, "darkblue"),
  font.tickslab = c(16, "plain", "black")
)

ggsurv$plot <- ggsurv$plot +
  theme(
    plot.title = element_text(hjust = 0.5),   # Center the title
    legend.text = element_text(size = 18, color = "black", face = "bold")
  )

ggsurv

