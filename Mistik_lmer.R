#####################################
###
### Mistik models fit with lmer()
###
###
#####################################


library(arm)
library(lsmeans)
library(car)
library(plyr)
library(reshape2)
library(xtable)
library(ggplot2)
library(pander)

source('mistik_load.R') # Bring in cleaned mistik objects

############################
### Mixed-effects models ###
############################

### Leaf Area ###

rplot<-function(model){
  plot(resid(model)~fitted(model),pch=20,cex=0.8)
  abline(0,0,lty=4)
}

la.mer<-lmer(log(sla)~species*gappos+(1|vegtrt),data=LA,REML=T)

rplot(la.mer)
qqnorm(resid(la.mer));qqline(resid(la.mer)) # Normality of error looks OK
plot(log(LA$sla)~fitted(la.mer));abline(0,1)

la.aov<-Anova(la.mer,type=2,test.statistic='F')
la.aov

row.names(la.aov)<-c('S','G','S G')

print(xtable(la.aov,
             caption='Results from Analysis of Deviance on differences in 
leaf area per unit wet weight for species and gap position. 
The standard deviation of the modeled random effect for transect was 6.0441, 
and the standard deviation of residual error was 21.6106'),type='html')
      


la.lsm<-lsmeans(la.mer,specs='pairwise~species*gappos')

la.plot<-ggplot(ldply(la.lsm[1])[-1],
                aes(x=gappos,y=exp(lsmean),
                    ymin=ifelse(exp(lsmean-SE)>0,exp(lsmean-SE),0),
                    ymax=exp(lsmean+SE),
                    fill=species))+
  geom_bar(stat='identity',position='dodge')+
  geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
  scale_fill_grey()+
  xlab("Gap Position")+
  ylab(expression(""*cm^2*"/g"))+
  theme_bw()+
  theme(axis.line=element_line(color='black'),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank(),
       panel.border=element_blank(),
       panel.background=element_blank())
la.plot

### Photosynthesis ###

photo.mer<-lmer(sqrt(photo)~light*species*gappos*vegtrt+(1|block/transect/gappos/treenum),
            data=mistik,REML=T)

rplot(photo.mer)
qqnorm(resid(photo.mer));qqline(resid(photo.mer))
qplot(resid(photo.mer),binwidth=0.1)
plot(sqrt(mistik$photo)~fitted(photo.mer));abline(0,1)

photo.aov<-Anova(photo.mer,type=2,test.statistic='F')
photo.aov
row.names(photo.aov)<-c('L','S','G','V','L S','L G','S G','L V','S V','G V','L S G','L S V','L G V','S G V','L S G V')

photo.lsm<-lsmeans(photo.mer,specs='pairwise~light*species*gappos|light:vegtrt:species')

photo.lsm
photo.plot<-ggplot(ldply(photo.lsm[1])[-1],aes(x=gappos,y=I(lsmean^2),
                                 ymin=ifelse(I((lsmean-SE)^2)>0,I((lsmean-SE)^2),0),
                                 ymax=I((lsmean+SE)^2),
                                 fill=species))+
  geom_bar(stat='identity',position='dodge')+
  geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
  facet_grid(.+light~vegtrt)+scale_fill_grey()+
  xlab("Gap Position")+ylim(0,25)+
  ylab(expression(" A ("*mu*"mol "*CO[2]*" "*m^-2*" "*s^-1*")"))+
  theme_bw()+
  theme(axis.line=element_line(color='black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color='black'),
        panel.background=element_blank())
photo.plot
# photo.plot2<-ggplot(ldply(photo.lsm[1])[-1],aes(x=gappos,y=I(lsmean^2),
#                                                ymin=ifelse(I((lower.CL)^2)>0,I((lower.CL)^2),0),
#                                                ymax=I((upper.CL)^2),
#                                                fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   facet_grid(.+light~vegtrt)+scale_fill_grey()+
#   xlab("Gap Position")+ylim(0,40)+
#   ylab(expression(" A ("*mu*"mol "*CO[2]*" "*m^-2*" "*s^-1*")"))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# photo.plot2


### Conductance ###

cond.mer<-lmer(sqrt(cond)~light*species*gappos*vegtrt+(1|block/transect/gappos/treenum),
               data=mistik,REML=T)
rplot(cond.mer) # may be a few outliers, but not sure whether to toss them or not...

qqnorm(resid(cond.mer));qqline(resid(cond.mer))
qplot(resid(cond.mer),binwidth=0.01)

cond.aov<-Anova(cond.mer,type=2,test.statistic='F')
cond.aov
row.names(cond.aov)<-row.names(photo.aov)

cond.lsm<-lsmeans(cond.mer,specs='pairwise~light*species*gappos|light:vegtrt:species')
cond.lsm

cond.plot<-ggplot(ldply(cond.lsm[1])[-1],aes(x=gappos,y=I(lsmean^2),
                                 ymin=ifelse(I((lsmean-SE)^2)>0,I((lsmean-SE)^2),0),
                                 ymax=I((lsmean+SE)^2),
                                 fill=species))+
  geom_bar(stat='identity',position='dodge')+
  geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
  facet_grid(.+light~vegtrt)+scale_fill_grey()+
  xlab("Gap Position")+ylim(0,0.5)+
  ylab(expression(""*g[sw]*" ("*mu*"mol "*CO[2]*" "*m^-2*" "*s^-1*")"))+
  theme_bw()+
  theme(axis.line=element_line(color='black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color='black'),
        panel.background=element_blank())
cond.plot

### Transpiration ###

trmmol.mer<-lmer(sqrt(trmmol)~light*species*gappos*vegtrt+(1|block/transect/gappos/treenum),
                 data=mistik,REML=T)

rplot(trmmol.mer)
qqnorm(resid(trmmol.mer));qqline(resid(trmmol.mer))
qplot(resid(trmmol.mer),binwidth=0.05)
plot(sqrt(mistik$trmmol)~fitted(trmmol.mer));abline(0,1)

trmmol.aov<-Anova(trmmol.mer,type=2,test.statistic='F')
trmmol.aov
row.names(trmmol.aov)<-row.names(photo.aov)


trmmol.mer2<-update(trmmol.mer,.~.-vegtrt-light:vegtrt-species:vegtrt-
                      gappos:vegtrt-light:species:vegtrt-light:gappos:vegtrt-
                      species:gappos:vegtrt-light:species:gappos:vegtrt)

trmmol.aov2<-Anova(trmmol.mer2,type=2,test.statistic='F')
trmmol.aov2

trmmol.lsm<-lsmeans(trmmol.mer2,specs='pairwise~light*species*gappos')
trmmol.lsm

trmmol.plot<-ggplot(ldply(trmmol.lsm[1])[-1],aes(x=gappos,y=I(lsmean^2),
                                 ymin=ifelse(I((lsmean-SE)^2)>0,I((lsmean-SE)^2),0),
                                 ymax=I((lsmean+SE)^2),
                                 fill=species))+
  geom_bar(stat='identity',position='dodge')+
  geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
  facet_grid(.~light)+scale_fill_grey()+
  xlab("Gap Position")+ylim(0,5)+
  ylab(expression(""*mu*"mol "*H[2]*"O "*m^2*" "*s^-1*""))+
  theme_bw()+
  theme(axis.line=element_line(color='black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color='black'),
        panel.background=element_blank())
trmmol.plot

# 
# # Ditch WUE? looks a little too crazy...
# ### WUE ###
# 
# wue.mer<-lmer(log(wue)~light*species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#               data=mistik,REML=T)
# rplot(wue.mer) # outliers over 2.5?
# qqnorm(resid(wue.mer));qqline(resid(wue.mer))
# qplot(resid(wue.mer),binwidth=0.1)
# plot(log(mistik$wue)~fitted(wue.mer));abline(0,1)
# 
# wue.aov<-Anova(wue.mer,type=2,test.statistic='F')
# wue.aov
# print(xtable(wue.aov,include.rownames=F),file='wue.tex',floating=F)
# 
# wue.lsm<-ldply(lsmeans(wue.mer,
#                          specs='pairwise~light*species*gappos|light:vegtrt:species')[1])[-1]
# wue.lsm
# wue.plot<-ggplot(wue.lsm,aes(x=gappos,y=exp(lsmean),
#                                  ymin=ifelse(exp(lsmean-SE)>0,exp(lsmean-SE),0),
#                                  ymax=exp(lsmean+SE),
#                                  fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   facet_grid(.+light~vegtrt)+scale_fill_grey()+
#   xlab("Gap Position")+ylim(0,40)+
#   ylab(expression(""*mu*"mol "*CO[2]*"/"*mu*"mol "*H[2]*"O"))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# wue.plot


### PPI ###


m1200<-mistik[mistik$PARi>1100,]
m50<-mistik[mistik$PARi<60,]
temp1<-ddply(m1200[,-c(3,5,37)],.(species,block,transect,vegtrt,gappos,treenum),
             summarize,photo1200=mean(photo))
temp2<-ddply(m50[-c(3,5,37)],.(species,block,transect,vegtrt,gappos,treenum),
             summarize,photo50=mean(photo))
plastic<-merge(temp1,temp2,by=c('block','transect','vegtrt','gappos','species','treenum'))
head(plastic)
plastic$plast<-(plastic$photo1200-plastic$photo50)/plastic$photo1200

plastic<-plastic[plastic$plast>0,]
plast.gg<-ggplot(plastic,aes(x=gappos,y=plast))+
  ggtitle("Photosynthetic Plasticity")+
  xlab("Gap Position")+
  ylab("PPI")

plast.gg+geom_point()+
  facet_grid(.~species*vegtrt)

plast.gg+geom_boxplot()+
  facet_grid(.~vegtrt)+
  aes(fill=species)

qplot(plast,data=plastic,binwidth=0.05)+facet_wrap(~species)
qplot(logit(plast),data=plastic,binwidth=0.3)+facet_wrap(~species)

plast.mer<-lmer(logit(plast)~species*gappos*vegtrt+(1|block/transect/gappos),
           data=plastic,REML=T)
rplot(plast.mer)
qqnorm(resid(plast.mer));qqline(resid(plast.mer)) # Normal-ish

plot(logit(plastic$plast)~fitted(plast.mer));abline(0,1) # Fit looks pretty good.

plast.aov<-Anova(plast.mer,type=2,test.statistic='F')
plast.aov # Nothing significant, do no further testing

table2<-cbind(photo.aov,cond.aov,trmmol.aov)
table2

# Random effects table
ran<-as.data.frame(rbind(sigma.hat(photo.mer)$sigma,
           sigma.hat(cond.mer)$sigma,
           sigma.hat(trmmol.mer)$sigma))[,c(5,4,3,2,1)]

names(ran)<-c('Block','Transect','Plot','Tree','Residual')

row.names(ran)<-c('Photosynthesis','Conductance','Transpiration')
ran # These are still standard deviations, not variances

# #######################
# ### Light submodels ###
# #######################
# 
# ### Photosynthesis 1200 ###
# 
# mer1<-lmer(photo~species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#                      data=m1200,REML=T)
# 
# rplot(mer1)
# qqnorm(resid(mer1));qqline(resid(mer1))
# qplot(resid(mer1,binwidth=0.1))
# abs(resid(mer1))>5
# m1200<-m1200[-c(357),] # Remove outlier
# 
# mer2<-update(mer1,data=m1200)
# rplot(mer2)
# qqnorm(resid(mer2));qqline(resid(mer1))
# qplot(resid(mer2),binwidth=0.1)
# 
# photo.1200.mer<-mer2
# 
# photo.1200.aov<-Anova(photo.1200.mer,type=2,test.statistic='F')
# photo.1200.aov
# 
# plot(m1200$photo~fitted(photo.1200.mer)); abline(0,1)
# 
# summary(photo.1200.mer)
# 
# photo.1200.lsm<-ldply(lsmeans(photo.1200.mer,specs='pairwise~species*gappos|vegtrt:species')[1])[-1]
# 
# photo.1200.plot<-ggplot(photo.1200.lsm,aes(x=gappos,y=lsmean,
#                                            ymin=ifelse(lsmean-SE>0,lsmean-SE,0),
#                                            ymax=lsmean+SE,
#                                            fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   facet_grid(.~vegtrt)+scale_fill_grey()+
#   xlab("Gap Position")+ylim(0,25)+
#   ylab(expression(""*mu*"mol "*CO[2]*""))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# 
# 
# photo.1200.plot
# 
# 
# ### Photosynthesis 50 ###
# 
# photo.50.mer<-lmer(photo~species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#            data=m50) 
# rplot(photo.50.mer)# log-transforming this model is equally bad
# qplot(resid(photo.50.mer),binwidth=0.1)
# qqnorm(resid(photo.50.mer));qqline(resid(photo.50.mer)) # pretty good
# 
# photo.50.aov<-Anova(photo.50.mer,type=2,test.statistic='F')
# photo.50.aov
# 
# photo.50.lsm<-ldply(lsmeans(photo.50.mer,specs='pairwise~gappos*species|vegtrt:species')[1])[-1]
# 
# photo.50.plot<-ggplot(photo.50.lsm,aes(x=gappos,y=lsmean,
#                                        ymin=ifelse(lsmean-SE>0,lsmean-SE,0),
#                                        ymax=lsmean+SE,
#                                        fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   facet_grid(.~vegtrt)+scale_fill_grey()+
#   xlab("Gap Position")+ylim(0,5)+
#   ylab(expression(""*mu*"mol "*CO[2]*""))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# 
# photo.50.plot
# 
# 
# ### Conductance 1200 ###
# 
# cond.1200.mer<-lmer(cond~species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#            data=m1200,REML=T)
# 
# rplot(cond.1200.mer)
# qqnorm(resid(cond.1200.mer));qqline(resid(cond.1200.mer))
# qplot(resid(cond.1200.mer),binwidth=0.01)
# abs(resid(cond.1200.mer))>0.15
# m1200<-m1200[-c(55,57,66,67,126),] # Remove outliers
# 
# mer2<-update(cond.1200.mer,data=m1200)
# rplot(mer2)
# qqnorm(resid(mer2));qqline(resid(mer2))
# qplot(resid(mer2),binwidth=0.01)
# 
# cond.1200.mer<-mer2
# 
# cond.1200.aov<-Anova(cond.1200.mer,type=2,test.statistic='F')
# cond.1200.aov
# 
# plot(m1200$cond~fitted(cond.1200.mer)); abline(0,1)
# 
# summary(cond.1200.mer)
# 
# cond.1200.lsm<-ldply(lsmeans(cond.1200.mer,specs='pairwise~species*gappos|vegtrt:species')[1])[-1]
# 
# cond.1200.plot<-ggplot(cond.1200.lsm,aes(x=gappos,y=lsmean,
#                                          ymin=ifelse(lsmean-SE>0,lsmean-SE,0),
#                                          ymax=lsmean+SE,
#                                          fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   facet_grid(.~vegtrt)+scale_fill_grey()+
#   xlab("Gap Position")+ylim(0,1)+
#   ylab(expression(""*mu*"mol "*H[2]*"O"))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# 
# 
# cond.1200.plot
# 
# 
# ### Conductance 50 ###
# 
# cond.50.mer<-lmer(cond~species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#                    data=m50) 
# rplot(cond.50.mer)# log-transforming this model is equally bad
# qplot(resid(cond.50.mer),binwidth=0.01)
# abs(resid(cond.50.mer))>0.15
# m50<-m50[-c(332,67),]
# 
# qqnorm(resid(cond.50.mer));qqline(resid(cond.50.mer)) # pretty good
# 
# cond.50.aov<-Anova(cond.50.mer,type=2,test.statistic='F')
# cond.50.aov
# 
# cond.50.lsm<-ldply(lsmeans(cond.50.mer,specs='pairwise~gappos*species|vegtrt:species')[1])[-1]
# 
# cond.50.plot<-ggplot(cond.50.lsm,aes(x=gappos,y=lsmean,
#                                      ymin=ifelse(lsmean-SE>0,lsmean-SE,0),
#                                      ymax=lsmean+SE,
#                                      fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   facet_grid(.~vegtrt)+scale_fill_grey()+
#   xlab("Gap Position")+ylim(0,1)+
#   ylab(expression(""*mu*"mol "*CO[2]*""))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# 
# cond.50.plot
# 
# 
# 
# ### Transpiration 1200 ###
# 
# tr.1200.mer<-lmer(trmmol~species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#                   data=m1200)
# 
# rplot(tr.1200.mer) # transforming causes similar problems
# qqnorm(resid(tr.1200.mer));qqline(resid(tr.1200.mer))
# qplot(resid(tr.1200.mer))
# plot(m1200$trmmol~fitted(tr.1200.mer));abline(0,1)
# 
# trmmol.1200.aov<-Anova(tr.1200.mer,type=2,test.statistic='F')
# trmmol.1200.aov
# 
# tr.1200.lsm<-lsmeans(update(tr.1200.mer,.~.-species:gappos:vegtrt-gappos:vegtrt-species:vegtrt-vegtrt),
#                      specs='pairwise~species:gappos')
# 
# tr.1200.plot<-ggplot(ldply(tr.1200.lsm[1])[-1],aes(x=gappos,y=lsmean,
#                                        ymin=ifelse(lsmean-SE>0,lsmean-SE,0),
#                                        ymax=lsmean+SE,
#                                        fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   scale_fill_grey()+xlab('Gap Position')+ylim(0,4)+
#   ylab(expression(""*mu*"mol "*H[2]*"O "*m^2*" "*s^-1*""))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# 
# tr.1200.plot
# 
# 
# ### Transpiration 50 ###
# 
# tr.50.mer<-lmer(trmmol~species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#                   data=m50)
# 
# 
# rplot(tr.50.mer)
# qqnorm(resid(tr.50.mer));qqline(resid(tr.50.mer))
# qplot(resid(tr.50.mer))
# plot(m50$trmmol~fitted(tr.50.mer));abline(0,1)
# 
# trmmol.50.aov<-Anova(tr.50.mer,type=2,test.statistic='F')
# trmmol.50.aov
# 
# 
# tr.50.lsm<-lsmeans(update(tr.50.mer,.~.-species:gappos:vegtrt-gappos:vegtrt-species:vegtrt-vegtrt),
#                    specs='pairwise~species:gappos')
# 
# tr.50.plot<-ggplot(ldply(tr.50.lsm[1])[-1],
#                    aes(x=gappos,y=lsmean,
#                        ymin=ifelse(lsmean-SE>0,lsmean-SE,0),
#                        ymax=lsmean+SE,
#                        fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   scale_fill_grey()+
#   xlab("Gap Position")+ylim(0,4)+
#   ylab(expression(""*mu*"mol "*H[2]*"O "*m^2*" "*s^-1*""))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# 
# tr.50.plot
# 
# 
# ### WUE 1200 ###
# 
# wue.1200.mer<-lmer(log(wue)~species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#                   data=m1200)
# 
# rplot(wue.1200.mer)
# qqnorm(resid(wue.1200.mer));qqline(resid(wue.1200.mer))
# qplot(resid(wue.1200.mer),binwidth=0.1)
# plot(log(m1200$wue)~fitted(wue.1200.mer));abline(0,1)
# 
# wue.1200.aov<-Anova(wue.1200.mer,type=2,test.statistic='F')
# wue.1200.aov
# 
# wue.1200.aov2<-Anova(update(wue.1200.mer,.~.-species:gappos:vegtrt-gappos:vegtrt-species:vegtrt-vegtrt),
#                             type=2,test.statistic='F')
# wue.1200.aov2
# 
# wue.1200.lsm<-lsmeans(update(wue.1200.mer,.~.-species:gappos:vegtrt-gappos:vegtrt-species:vegtrt-vegtrt),
#                              specs='pairwise~species:gappos')
# 
# wue.1200.plot<-ggplot(ldply(wue.1200.lsm[1])[-1],
#                       aes(x=gappos,y=exp(lsmean),
#                           ymin=ifelse(exp(lsmean-SE)>0,exp(lsmean-SE),0),
#                           ymax=exp(lsmean+SE),
#                           fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   scale_fill_grey()+xlab('Gap Position')+
#   ylab(expression(""*mu*"mol "*CO[2]*"/"*mu*"mol "*H[2]*"O"))+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# 
# 
# wue.1200.plot
# 
# ### WUE 50 ###
# 
# wue.50.mer<-lmer(log(wue)~species*gappos*vegtrt+(1|block/transect/gappos/treenum),
#                  data=m50)
# 
# rplot(wue.50.mer)
# qqnorm(resid(wue.50.mer));qqline(resid(wue.50.mer))
# plot(log(m50$wue)~fitted(wue.50.mer));abline(0,1)
# qplot(resid(wue.50.mer,binwidth=0.1))
# summary(wue.50.mer)
# 
# wue.50.aov<-Anova(wue.50.mer,type=2,test.statistic='F')
# wue.50.aov
# 
# wue.50.aov2<-Anova(update(wue.50.mer,.~.-species:gappos:vegtrt-gappos:vegtrt-species:vegtrt-vegtrt),
#                    type=2,tes.statistic='F')
# wue.50.aov2
# 
# # Plot Multiple comparisons
# # wue.50.lsm<-lsmeans(wue.50.mer,specs='pairwise~species*gappos*vegtrt|vegtrt')
# wue.50.lsm<-lsmeans(update(wue.50.mer,.~.-species:gappos:vegtrt-gappos:vegtrt-species:vegtrt-vegtrt),
#                       specs='pairwise~species:gappos')
# 
# wue.50.plot<-ggplot(ldply(wue.50.lsm[1])[-1],
#                     aes(x=gappos,y=exp(lsmean),
#                         ymin=ifelse(exp(lsmean-SE)>0,exp(lsmean-SE),0),
#                         ymax=exp(lsmean+SE),
#                         fill=species))+
#   geom_bar(stat='identity',position='dodge')+
#   geom_errorbar(width=0.4,position=position_dodge(width=0.9))+
#   scale_fill_grey()+xlab('Gap Position')+ylim(0,5)+
#   ylab(expression(""*mu*"mol "*CO[2]*"/"*mu*"mol "*H[2]*"O"))+
# #   facet_wrap(~vegtrt)+
#   theme_bw()+
#   theme(axis.line=element_line(color='black'),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_rect(color='black'),
#         panel.background=element_blank())
# wue.50.plot
# 

