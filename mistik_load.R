##################################
#### Load and clean Mistik data
##################################

rm(list=ls())

library(ggplot2)
library(car)

###################################################
### Define data objects and remove logging errors
###################################################

### Gas Exchange ###

plots<-read.table('./Mistik_Data/plots.txt',sep='\t',header=T)
head(plots)

trees<-read.table('./Mistik_Data/trees.txt',sep='\t',header=T)
head(trees)

obs<-read.table('./Mistik_Data/obs_raw.txt',sep='\t',header=T)
head(obs)

obs.trees<-merge(obs,trees,by='treenum')
mistik.load<-merge(obs.trees,plots,by='transect')
summary(mistik.load)

mistik<-mistik.load[mistik.load$photo>0&mistik.load$cond>0& # Delete neg vals, ALIN
                      mistik.load$trmmol>0&mistik.load$species!='ALIN',]

mistik$light<-factor(ifelse(mistik$PARi>1100,'1200 mmol PPFD','50 mmol PPFD'))
mistik$gap<-factor(mistik$gap,ordered=F)
mistik$transect<-factor(mistik$transect,ordered=F)
mistik$species<-factor(mistik$species)
mistik$treenum<-factor(mistik$treenum,ordered=F)
mistik$ind<-factor(mistik$ind,ordered=F)
mistik$gappos<-factor(mistik$gappos,levels=c('SU','SE','C','NE','NU'))
mistik$vegtrt<-factor(ifelse(mistik$vegtrt=='r','Brush Saw','Meri-Crusher'))
mistik$photo<-mistik$photo*(mistik$Area/6) # Correct for leaf area
mistik$cond<-mistik$cond*(mistik$Area/6)
mistik$trmmol<-mistik$trmmol*(mistik$Area/6)

# Fix date and time
mistik$date.time<-paste(substring(mistik$date,1,9),substring(mistik$HMS,12))
mistik$date<-strptime(mistik$date.time,format='%m/%d/%Y %H:%M:%S',tz='Canada/Saskatchewan')
mistik$cont.date<-as.numeric(mistik$date)+mistik$time
mistik<-mistik[,-c(4,38)] # Drop unused date and time columns
head(mistik)

mistik.raw<-mistik
row.names(mistik)<-1:nrow(mistik)

### Leaf Area ###

LA<-read.table("./Mistik_data/leaf_area.txt", sep="\t", header=T, row.names=1)

LA$gappos<-recode(LA$gappos, "c('SU','NU')='U'")
LA$gappos<-factor(LA$gappos, levels=c("U","SE","C","NE"))

LA$sla<-LA$area/LA$weight

summary(LA)

# Examine data distribution, clean data


#################################################
### Exploratory analysis and removal of outliers
#################################################

### Leaf Area ###

qplot(x=sla,data=LA,binwidth=10)+facet_wrap(~species)
LA[LA$sla>200,]

LA.gg<-ggplot(LA, aes(factor(gappos),y=sla))+
  ylab("SLA")+
  xlab("Gap Position")
LA.gg+geom_boxplot()+facet_grid(.~species)
LA.gg+geom_point()+facet_grid(.~species)
LA<-LA[-c(103),] # Remove outlier

### Photosynthesis ###

qplot(x=photo,data=mistik,binwidth=1)+facet_wrap(~species+light)
# Looks like a gamma distribution?

photo.gg<-ggplot(mistik,aes(factor(gappos),y=photo))+
  ggtitle(expression(""*A[net]*""))+
  xlab("Gap Position")+
  ylab(expression(""*mu*"mol "*CO[2]*" "*m^-2*" "*s^-1*""))

photo.gg+geom_point()+
  facet_grid(.+light~species*vegtrt)
photo.gg+geom_boxplot()+
  facet_grid(.+light~species)+
  aes(fill=factor(vegtrt))

### Stomatal Conductance ###

qplot(x=cond,data=mistik,binwidth=0.03)+facet_wrap(~species+light)
# Gamma distribution again

cond.gg<-ggplot(mistik,aes(x=gappos,y=photo))+
  ggtitle('Conductance')+
  xlab('gap position')+
  ylab(expression(""*mu*"mol "*H[2]*"O "*m^-2*" "*s^-1*""))
cond.gg+geom_point()+
  facet_grid(.+light~species*vegtrt)
cond.gg+geom_boxplot()+
  facet_grid(.+light~species)+
  aes(fill=vegtrt)

### Transpiration ###

qplot(x=trmmol,data=mistik,binwidth=0.2)+facet_wrap(~species+light)
# Gamma distribution?

trmmol.gg<-ggplot(mistik,aes(factor(gappos),y=trmmol))+
  ggtitle('E')+
  xlab('Gap Position')+
  ylab(expression(""*mu*"mol "*H[2]*"0 "*m^-2*" "*s^-1*""))

trmmol.gg+geom_point()+
  facet_grid(.+light~species*vegtrt)
trmmol.gg+geom_boxplot()+
  facet_grid(.+light~species)+
  aes(fill=factor(vegtrt))

