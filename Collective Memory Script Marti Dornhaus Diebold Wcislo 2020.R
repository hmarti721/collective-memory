setwd("C:/Users/hanna/Documents/Dissertation/Collective Memory paper for BEAS/Data analysis")
allTreatments_glmm <- read.csv(file = "colmemory18June.csv", header = TRUE, sep = ',')
head(allTreatments_glmm)

library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(emmeans)

#Remove Set V (both leaf types removed with forceps) for main analysis
allTreatments_glmm%>%
  filter(Set!= "V")-> no_double_glmm

# Create 2 datasets for time points when naive ants are 80-99% and 100% of their colonies, respectively
turnover80_99<-filter(no_double_glmm, Prop.new>0.799 & Prop.new<1.0)
turnover100<-filter(no_double_glmm, Prop.new == 1.0)

# Do cycloheximide treated colonies still take fewer leaves of the previously treated type when all experienced ants are dead?
mod100<- glmer.nb(PT_taken~
                    Treatment+
                    Leaf_removal+
                    offset(ln_LeavesTaken)+
                    (1|Week)+
                    (1|Colony)+
                    (1|Experimental_leaf)+
                    (1|Queen), data = turnover100)
summary(mod100)
## No effect of cycloheximide, but there is an effect of leaf removal pattern

# Do cycloheximide treated colonies still take fewer leaves of the previously treated type when colonies are 80-99% naive ants?
mod80_99.new<- glmer.nb(PT_taken~
                      Treatment+
                      Leaf_removal+
                      offset(ln_LeavesTaken)+
                      (1|Week)+
                      (1|Colony)+
                        (1|Experimental_leaf)+
                      (1|Queen), data = filter(turnover80_99, Agegroup=="new"))
summary(mod80_99.new)


turnover80_99n<-filter(turnover80_99, Agegroup=="new")
summary(turnover80_99n$PT_taken)
emmeans(mod80_99.new, "Treatment", type="response")
emmeans(mod80_99.new, "Leaf_removal", type="response")
# Yes, memory retained when colonies are 80-99% naive ants

mod80_99.old<- glmer.nb(PT_taken~
                          Treatment+
                          Leaf_removal+
                          offset(ln_LeavesTaken)+
                          (1|Week)+
                          (1|Colony)+
                          (1|Experimental_leaf)+
                          (1|Queen), data = filter(turnover80_99, Agegroup=="old"))
summary(mod80_99.old)

turnover80_99o<-filter(turnover80_99, Agegroup=="old")
summary(turnover80_99o$PT_taken)
emmeans(mod80_99.old, "Treatment", type="response")
emmeans(mod80_99.old, "Leaf_removal", type="response")

##Overall, are do colonies remember cycloheximide treatment (Treatment)? Do they relearn that prev treated leaf type is beneficial when they can incorporate it (Leaf_removal)?
mod.all <- glmer.nb(PT_taken~
                      Treatment+
                      Leaf_removal+
                      offset(ln_LeavesTaken)+
                      (1|Week)+
                      (1|Colony)+
                      (1|Experimental_leaf)+
                      (1|Queen), data = no_double_glmm)
summary(mod.all)
emmeans(mod.all, "Treatment", type="response")
emmeans(mod.all, "Leaf_removal", type="response")

##Overall, are do naive ants conform to preferences of experienced nestmates (Treatment)? Do they relearn that prev treated leaf type is beneficial when they can incorporate it (Leaf_removal)?
mod.new <- glmer.nb(PT_taken~
                      Treatment+
                      Leaf_removal+
                      offset(ln_LeavesTaken)+
                      (1|Week)+
                      (1| Colony)+
                      (1|Experimental_leaf)+
                      (1|Queen), data = filter(no_double_glmm, Agegroup=="new"))
summary(mod.new)
emmeans(mod.new, "Treatment", by="Leaf_removal", type="response")

##Overall, are do experienced ants remember cycloheximide treatment (Treatment)? Do they relearn that prev treated leaf type is beneficial when they can incorporate it (Leaf_removal)?
mod.old <- glmer.nb(PT_taken~
                      Treatment+
                      Leaf_removal+
                      offset(ln_LeavesTaken)+ 
                      (1|Week)+
                      (1| Colony)+
                      (1|Experimental_leaf)+
                      (1|Queen), data = filter(no_double_glmm, Agegroup=="old"))
summary(mod.old)
emmeans(mod.old, "Treatment", by="Leaf_removal", type="response")



####### Calculate the effect of past leaf incorporation on leaf choice
mergeEorN1 <- read.csv(file = "mergeEorN1.csv", header = TRUE, sep = ',')
leaves_incorp <- read.csv(file = "pastleafincorp.csv", header = TRUE, sep = ',')

## Calculate the number of leaves of each type incorporated into the garden based on the number of each leaf type taken - the number of each type removed from the colony after 2-6 hours had passed
leaves_incorp%>%
  mutate(LS_taken=Leaves_taken-IB_taken)%>%
  mutate(LS_incorp=LS_taken-LS.removed, IB_incorp=IB_taken-IB.removed)%>%
  mutate(total_incorp=LS_incorp+IB_incorp)%>%
  mutate(PropIB_incorp_lastweek=IB_incorp/total_incorp)%>%
  group_by(Colony)%>%
  arrange(Colony,Week)%>%
  mutate(LS_incorp_ever=cumsum(LS_incorp), IB_incorp_ever=cumsum(IB_incorp), total_incorp_ever=cumsum(total_incorp))->leaves_incorp1

###NEXT: select only the necessary columns, make the time shift variable, make proportion previously treated incorporated ever
# remerge the cumsum columns onto a df that has data separated by age group
leaves_incorp1%>%
  select(Colony, Week, LS_incorp_ever, IB_incorp_ever, total_incorp_ever, Treatment, PropIB_incorp_lastweek)%>%
  mutate(Prop_PT_incorp_ever=case_when(Treatment=="CHX IB"|Treatment=="Sham IB"~ IB_incorp_ever/total_incorp_ever,
                                       Treatment=="CHX LS"| Treatment=="Sham LS" ~ LS_incorp_ever/total_incorp_ever),
         PT_incorp_ever = Prop_PT_incorp_ever*total_incorp_ever,
         PropPT_incrop_lastweek= case_when(Treatment=="CHX IB"|Treatment=="Sham IB"~ PropIB_incorp_lastweek,
                                           Treatment=="CHX LS"| Treatment=="Sham LS" ~ 1-PropIB_incorp_lastweek))->leaves_incorp2


#Then make a time plus 1 column to shift the cumulative column so that I can compare the leaves taken at time x with leaves incorporated before time x
leaves_incorp2%>%
  mutate(tplus1=case_when(Week==0~1,
                          Week >0 ~ Week+2))%>%
  rename(Week_actual = Week)%>%
  rename(Week=tplus1)->leaves_incorp3

cumsum<-merge(mergeEorN1, leaves_incorp3, by=c("Colony", "Week", "Treatment"))

cumsum%>%
  mutate(EorN=case_when(Age=="new"~"Naive",
                        Age=="old"~"Experienced"))->cumsum
###Do we still see an effect of treatment after accounting for the effect of leaves incorporated
cumsum %>%
  mutate(ln_LeavesTaken = log(Leaves_taken.a),
         PT_taken = PT_prop.a*Leaves_taken.a)-> cumsum_glmm



#Test for effect of past leaf incorporation (PT_incorp_ever) and cycloheximide treatment (Treatment_b) on leaf choices of old, experienced ants
# These results are in table 3
cumsum_glmm%>%na.omit()->cumsum_glmm
mod.oldincorp <- glmer.nb(PT_taken~
                      Treatment_b+
                      PT_incorp_ever+
                      offset(ln_LeavesTaken)+
                      (1|Week)+
                      (1| Colony), data = filter(cumsum_glmm, Age=="old"))
summary(mod.oldincorp)

#Calculate effect sizes as estimated marginal means
emmeans(mod.oldincorp, "Treatment_b", type="response")
emmeans(mod.oldincorp, "PT_incorp_ever", type="response", cov.reduce = range)

#Test for effect of past leaf incorporation (PT_incorp_ever) and cycloheximide treatment (Treatment_b) on leaf choices of new, naive ants
mod.newincorp <- glmer.nb(PT_taken~
                      Treatment_b+
                      PT_incorp_ever+
                      offset(ln_LeavesTaken)+
                      (1|Week)+
                      (1| Colony), data = filter(cumsum_glmm, Age=="new"))
summary(mod.newincorp)

#Calculate effect sizes as estimated marginal means
emmeans(mod.newincorp, "Treatment_b", type="response")
emmeans(mod.newincorp, "PT_incorp_ever", type="response", cov.reduce = range)



#####################################################################################
# Test the effects of outbound contacts, presence of tiny ants, and previous IB attempts on  leaf choice
# Results for Table 4

data.full <- read.csv(file = "videoscoringdata.csv", header = TRUE, sep = ',')

#Logistic regression
mod.control <- glmer(Choice_bin~
                       TouchIBtiny+
                       Outbound_propIB+
                       prev_IB_attempts+
                       (1|Week)+
                       (1|Colony), data = filter(data.full, Treatment=="Control"), family=binomial)
summary(mod.control)



mod.chx<-  glmer(Choice_bin~
                   TouchIBtiny+
                   Outbound_propIB+
                   prev_IB_attempts+
                   (1|Week)+
                   (1|Colony), data =filter(data.full, Treatment=="CHX"), family=binomial)
summary(mod.chx)



#Make column called 'leaf type' 
data.full%>%
  gather(key="leaf_type", value = "tinies", TouchIBtiny, TouchLStiny)->data.full1

#Are there more tinies on IB in CHX treated colonies compared to IB in Control colonies?
wilcox.test(tinies~Treatment, data=filter(data.full1, leaf_type=="TouchIBtiny"))

#Are there more tinies on IB compared to LS in CHX treated colonies?
wilcox.test(tinies~leaf_type, data=filter(data.full1, Treatment=="CHX"))



data.full1%>%
  mutate(Choice_Condition=case_when(Treatment=="CHX"&Choice_bin==0~"LS in Treated Colony",
                                    Treatment=="CHX"&Choice_bin==1~"IB in Treated Colony",
                                    Treatment=="Control"&Choice_bin==0~"LS in Untreated Colony",
                                    Treatment=="Control"&Choice_bin==1~"IB in Untreated Colony",))-> data.full1


StatN <- ggproto("StatN", Stat,
                 required_aes = c("x", "y"), 
                 compute_group = function(data, scales) {
                   y <- data$y
                   y <- y[!is.na(y)]
                   n <- length(y)
                   data.frame(x = data$x[1], y = -0.1, label = n)
                 }
)
stat_n <- function(mapping = NULL, data = NULL, geom = "text", 
                   position = "identity", inherit.aes = TRUE, show.legend = NA, 
                   na.rm = FALSE, ...) {
  ggplot2::layer(stat = StatN, mapping = mapping, data = data, geom = geom, 
                 position = position, inherit.aes = inherit.aes, show.legend = show.legend, 
                 params = list(na.rm = na.rm, ...))
}

ggplot(data.full1, aes(x=Choice_Condition, y= tinies))+
  geom_boxplot()+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  stat_n()


####################### Do ants carrying IB in treated colonies recieve more inbound contacts?


contacts <- read.csv(file = "contacts.csv", header = TRUE, sep = ',')


ggplot(contacts, aes(x=Choice_Condition, y=inbound.antennation))+
  geom_boxplot()+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("# of antennations on inbound journey")+
  xlab("Leaf chosen by forager")+
  stat_n(position=position_nudge(y=-0.5))


############ Do ants carrying IB in cycloheximide treated colonies receive a different amount of antennations than ants carrying LS in treated colonies?
wilcox.test(inbound.antennation~Choice, data=filter(contacts, Treatment=="CHX"))

######### Do ants carrying IB in CHX colonies receive a different amount of antennations than ants carrying IB in control colonies?
wilcox.test(inbound.antennation~Treatment, data=filter(contacts, Choice=="IB"))







##### Generate figures

#Define color palette
#This color palette is for the paper figures, can be printed in greyscale
GreycompPalette<-c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4", "grey50")


# Fig 3a
ggplot(no_double_glmm, aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  scale_color_manual(values=GreycompPalette)+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choice")+
  xlab("Worker Turnover")+
  ylim(0,1.0)

### Fig 3b effect of cycloheximide treatment when colonies are 100% naive ants
ggplot(turnover100, aes(x=Set, y=Prop_PT, fill=Set))+
  geom_boxplot()+
  scale_fill_manual(values=GreycompPalette)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choice at 100% Turnover")+
  xlab("Treatment")+
  stat_n()

# Fig 4a
ggplot(filter(no_double_glmm, Agegroup=="old"), aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  scale_color_manual(values=GreycompPalette)+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Experienced Ants")+
  xlab("Worker Turnover")+
  ylim(0,1.0)

# Fig 4b boxplots, mean leaf choice for experienced ants in treatments I-IV
ggplot(filter(no_double_glmm, Agegroup=="old"), aes(x=Set, y=Prop_PT, fill=Set))+
  scale_fill_manual(values=GreycompPalette)+
  geom_boxplot()+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Experienced Ants")+
  xlab("Treatment")+
  stat_n()+
  ylim(-0.1,1.0)

# Fig 5a
ggplot(filter(no_double_glmm, Agegroup=="new"), aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  scale_color_manual(values=GreycompPalette)+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choice of Naive Ants")+
  xlab("Worker Turnover")+
  ylim(0,1.0)

# Fig 5b boxplots, mean leaf choice for naive ants in treatments I-IV
ggplot(filter(no_double_glmm, Agegroup=="new"), aes(x=Set, y=Prop_PT, fill=Set))+
  scale_fill_manual(values=GreycompPalette)+
  geom_boxplot()+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Naive Ants")+
  xlab("Treatment")+
  stat_n()+
  ylim(-0.1,1.0)

# Fig for ppt 80-99% Naive only
ggplot(filter(turnover80_99, Agegroup=="new"), aes(x=Set, y=Prop_PT, fill=Set))+
  scale_fill_manual(values=GreycompPalette)+
  geom_boxplot()+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Naive Ants")+
  xlab("Treatment")+
  stat_n()+
  ylim(-0.1,1.0)

# Fig for ppt 80-99% Naive, Exp only
ggplot(filter(turnover80_99, Agegroup=="old"), aes(x=Set, y=Prop_PT, fill=Set))+
  scale_fill_manual(values=GreycompPalette)+
  geom_boxplot()+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Experienced Ants")+
  xlab("Treatment")+
  stat_n()+
  ylim(-0.1,1.0)

# Figs for ppt building treatments, exp and naive
ggplot(filter(no_double_glmm, Set=="I"), aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  scale_color_manual(values=GreycompPalette)+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choice")+
  xlab("Worker Turnover")+
  ylim(0,1.0)

# Figs for ppt building treatments, exp and naive
ggplot(filter(no_double_glmm, Leaf_removal=="IB removal"), aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  scale_color_manual(values=GreycompPalette)+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choice")+
  xlab("Worker Turnover")+
  ylim(0,1.0)
GreycompPalette<-c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4", "grey50")

# Figs for ppt building treatments, exp and naive
ggplot(filter(no_double_glmm, Set=="III"), aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  scale_color_manual(values="#41b6c4")+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choice")+
  xlab("Worker Turnover")+
  ylim(0,1.0)

# Figs for ppt building treatments, exp and naive
ggplot(filter(no_double_glmm, Leaf_removal=="no removal"), aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  scale_color_manual(values=c("#41b6c4","#a1dab4"))+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choice")+
  xlab("Worker Turnover")+
  ylim(0,1.0)

# Figs for ppt building treatments, exp and naive
ggplot(filter(cumsum_glmm, Age=="old"), aes(x=Prop_PT_incorp_ever, y=PT_prop.a, color=Treatment_b))+
  geom_point(size=3)+
  scale_color_manual(values=c("#41b6c4","#a1dab4"))+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf choices of experienced ants")+
  xlab("Past leaf incorporation")+
  ylim(0,1.0)


#Supplemental fig 1a
ggplot(filter(allTreatments_glmm, Agegroup=="old"), aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Experienced Ants")+
  scale_color_manual(values=GreycompPalette, name="Treatment", labels = c("I, Treated, one leaf type removed", "II, Treated, no leaves removed","III, Untreated, one leaf type removed", "IV, Untreated, no leaves removed", "V, Untreated, both leaves removed"))+
  xlab("Worker Turnover")+
  ylim(0,1.0)





### For supplementary materials
##Glmms with all treatments, including the double leaf removal treatment
mod.new5 <- glmer.nb(PT_taken~
                      Treatment+
                      Leaf_removal+
                      offset(ln_LeavesTaken)+
                      (1|Week)+
                      (1| Colony)+
                       (1|Year)+
                      (1| Queen), data = filter(allTreatments_glmm, Agegroup=="new"))
summary(mod.new5)

mod.old5 <- glmer.nb(PT_taken~
                      Treatment+
                      Leaf_removal+
                      offset(ln_LeavesTaken)+ 
                      (1|Week)+
                      (1| Colony)+
                       (1|Year)+
                       (1| Queen), data = filter(allTreatments_glmm, Agegroup=="old"))
summary(mod.old5)

# Supplemental fig 1b
ggplot(filter(allTreatments, Agegroup=="old"), aes(x=Set, y=Prop_PT, fill=Set))+
  scale_fill_manual(values=GreycompPalette)+
  geom_boxplot()+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Experienced Ants")+
  xlab("Treatment")+
  stat_n()+
  ylim(-0.1,1.0)

#Supplemental fig 2a
ggplot(filter(allTreatments, Agegroup=="new"), aes(x=Prop.new, y=Prop_PT, color=Set))+
  geom_point()+
  scale_color_manual(values=GreycompPalette)+
  geom_smooth(method='lm',formula=y~x, se =FALSE,aes(weight=Leaves_taken), size=3)+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Naive Ants")+
  xlab("Worker Turnover")+
  ylim(0,1.0)

# Supplemental fig2b
ggplot(filter(allTreatments, Agegroup=="new"), aes(x=Set, y=Prop_PT, fill=Set))+
  geom_boxplot()+
  theme_bw(base_size=20)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Leaf Choices of Naive Ants")+
  xlab("Treatment")+
  stat_n()+
  scale_fill_manual(values=GreycompPalette)+
  #scale_fill_manual(name="Treatment", labels = c("I, Treated, one leaf type removed", "II, Treated, no leaves removed","III, Untreated, one leaf type removed", "IV, Untreated, no leaves removed", "V, Untreated, both leaves removed"))+
  ylim(-0.1,1)




