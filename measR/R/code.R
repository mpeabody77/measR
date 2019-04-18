####################################################################################
#                                                                                  #
#                           Read in DATA                                           #
#                                                                                  #
####################################################################################

data<-read.csv("Exam1_All.csv", na.strings=c("","NA"))

####################################################################################
#                                                                                  #
#                           Model Specifications                                   #
#                                                                                  #
####################################################################################


#     measR(data,
#           model=
#           iselect=,
#           pselect=,
#           iaf=, 
#           paf=, 
#           idf=, 
#           pdf=,
#           item.start=, 
#           n.items=, 
#           n.people=, 
#           converg=, 
#           max.iter=)



item.start<-3
n.items<-18
n.people<-35
converg<-abs(.01)
#add max.iter?

demographics<-data[, 1:item.start-1]
demographics$entry<-as.numeric(row.names(demographics))
data2<-data[,item.start:((item.start-1)+n.items)]
#data2<-colnames()



####################################################################################
#                                                                                  #
#                           Remove Extreme Values                                  #
#                                                                                  #
####################################################################################

repeat {
  person.raw<-apply(data2, 1, sum, na.rm=TRUE)
  person.admin<-apply(data2, 1, function(x){sum(!is.na(x))})
  person.prop<-person.raw/person.admin
  data2<-data2[c(person.prop!=1 & person.prop!=0),]
  
  item.raw<-apply(data2, 2, sum, na.rm=TRUE)
  item.admin<-apply(data2, 2, function(x){sum(!is.na(x))})
  item.prop<-item.raw/item.admin
  data2<-data2[, c(item.prop!=1 & item.prop!=0)]
  
  person.raw<-apply(data2, 1, sum, na.rm=TRUE)
  person.admin<-apply(data2, 1, function(x){sum(!is.na(x))})
  person.prop<-person.raw/person.admin
  item.raw<-apply(data2, 2, sum, na.rm=TRUE)
  item.admin<-apply(data2, 2, function(x){sum(!is.na(x))})
  item.prop<-item.raw/item.admin
  
  if (length(which(person.prop==1 | person.prop==0 |
                 item.prop==1 | item.prop==0))==0)  {
                   break
                 }
}

####################################################################################
#                                                                                  #
#                                  JMLE Iteration                                  #
#                                                                                  #
####################################################################################

person.ability<-log(person.prop/(1-person.prop))
item.difficulty<-log((1-item.prop)/item.prop)
avg.item.difficulty<-mean(item.difficulty)
adj.item.difficulty<-item.difficulty-avg.item.difficulty


repeat { 
  expected.values<-function(x,y){(exp(x-y)/(1+exp(x-y) ) ) }
  Exp.Val<-outer(person.ability, adj.item.difficulty, expected.values)
  exp.var<-function(x){x*(1-x)}
  var<-exp.var(Exp.Val)
  
  person.var<-apply(var, 1, sum,na.rm=TRUE)
  person.var<-person.var*-1
  item.var<-apply(var, 2, sum, na.rm=TRUE)
  item.var<-item.var*-1
  
  res<-data2-Exp.Val
  #res<-data2[,(1:ncol(data2))]-Exp.Val
  person.sum.of.residuals<-apply(res, 1, sum, na.rm=TRUE)
  item.sum.of.residuals<-apply(res, 2, sum, na.rm=TRUE)
  item.sum.of.residuals<-(-1*item.sum.of.residuals)
  
  per.square.resid<- sapply(person.sum.of.residuals, function(x) x^2)
  per.sum.of.square.resid<-sum(per.square.resid)
  item.square.resid<- sapply(item.sum.of.residuals, function(x) x^2)
  item.sum.of.square.resid<-sum(item.square.resid)
  
  new.ability<-person.ability-person.sum.of.residuals/person.var
  new.difficulty<-adj.item.difficulty-item.sum.of.residuals/item.var
  avg.new.difficulty<-mean(new.difficulty)
  adj.new.difficulty<-new.difficulty-avg.new.difficulty
  
  person.ability<-new.ability
  adj.item.difficulty<-adj.new.difficulty
  
  if(per.sum.of.square.resid<=converg) {
    break
  }
  }


####################################################################################
#                                                                                  #
#                               FIT STATISTICS                                     #
#                                                                                  #
####################################################################################

#OUTFIT MNSQ
Zscore<- res / (var^(1/2))
Zsquared<-(res^2) / var
per.sum.sq.resid<-apply(Zsquared, 1, sum, na.rm=TRUE)
item.sum.sq.resid<-apply(Zsquared, 2, sum, na.rm=TRUE)
n.meas.person<-apply(Zsquared, 2, length)
n.meas.item<-apply(Zsquared, 1, length)
Person.Outfit.mnsq<-per.sum.sq.resid/n.meas.item
Item.Outfit.mnsq<-item.sum.sq.resid/n.meas.person

#INFIT MNSQ
In<-Zsquared*var
Per.SumIn<-apply(In, 1, sum, na.rm=TRUE)
Per.SumVar<-apply(var, 1, sum, na.rm=TRUE)
Per.Infit.mnsq<-Per.SumIn/Per.SumVar
Item.SumIn<-apply(In, 2, sum, na.rm=TRUE)
Item.SumVar<-apply(var, 2, sum, na.rm=TRUE)
Item.Infit.mnsq<-Item.SumIn/Item.SumVar

#ZSTD from p100 Wright & Masters
#Eni is the Expected mean of Xni
Eni<-(0*(1-Exp.Val)) + (1*Exp.Val)
#Wni is the Variance of Xni
Wni<-(((0-Eni)^2)*(1-Exp.Val)) + (((1-Eni)^2) * Exp.Val)
#Cni is the Kurtosis of Xni
Cni<-(((0-Eni)^4)*(1-Exp.Val)) + (((1-Eni)^4) * Exp.Val)

#Person INFIT ZSTD
Per.sum.Wni<-apply(Wni, 1, sum, na.rm=TRUE)
Per.In.Zstd.var<-(Cni-Wni^2)/Per.sum.Wni^2
sum.Per.In.Zstd.var<-apply(Per.In.Zstd.var, 1, sum, na.rm=TRUE)

per.in.qi<-sqrt(sum.Per.In.Zstd.var)
per.in.qi[per.in.qi>sqrt(2)]<-sqrt(2)
Per.In.ZSTD<-(Per.Infit.mnsq^(1/3)-1)*(3/per.in.qi)+(per.in.qi/3)

#Person OUTFIT ZSTD 
Per.Out.Wni2<-Wni^2
#2.  the contribution to the model variance of the OUTFIT MEANSQ 
#    by one observation is not allowed to overwhelm the contributions 
#    of the other observations. Wni^2 is not allowed to be less than 0.00001
Per.Out.Wni2[Per.Out.Wni2<.00001]<-.00001
Per.CniWni<-(Cni/Per.Out.Wni2)/n.meas.item^2
sum.Per.Out.Zstd.CniWni<-apply(Per.CniWni, 1, sum, na.rm=TRUE) #Sum CniWni
sum.Per.Out.Zstd.var<-sum.Per.Out.Zstd.CniWni-(1/n.meas.item)  #SumVar / df
#1. the lower limit of the degrees of freedom of the mean-square statistic 
#is set at 1. qi is not allowed to be more than sqrt(2) for both. 
per.out.qi<-sqrt(sum.Per.Out.Zstd.var)
per.out.qi[per.out.qi>sqrt(2)]<-sqrt(2)
Per.Out.ZSTD<-(Person.Outfit.mnsq^(1/3)-1)*(3/per.out.qi)+(per.out.qi/3)


#Item INFIT ZSTD 
t.Wni<-t(Wni)
t.Cni<-t(Cni)
Item.sum.Wni<-apply(t.Wni, 1, sum, na.rm=TRUE)
Item.In.Zstd.var<-(t.Cni-t.Wni^2) / Item.sum.Wni^2 #this is the problem
sum.Item.In.Zstd.var<-apply(Item.In.Zstd.var, 1, sum, na.rm=TRUE)

Item.in.qi<-sqrt(sum.Item.In.Zstd.var)
Item.in.qi[Item.in.qi>sqrt(2)]<-sqrt(2)
Item.In.ZSTD<-(Item.Infit.mnsq^(1/3)-1)*(3/Item.in.qi)+(Item.in.qi/3)


#Item OUTFIT ZSTD
Item.Out.Wni2<-Wni^2
#2.  the contribution to the model variance of the OUTFIT MEANSQ 
#    by one observation is not allowed to overwhelm the contributions 
#    of the other observations. Wni^2 is not allowed to be less than 0.00001
Item.Out.Wni2[Item.Out.Wni2<.00001]<-.00001
Item.CniWni<-(Cni/Item.Out.Wni2)/n.meas.person^2
sum.Item.Out.Zstd.CniWni<-apply(Item.CniWni, 2, sum, na.rm=TRUE) #Sum CniWni
sum.Item.Out.Zstd.var<-sum.Item.Out.Zstd.CniWni-(1/n.meas.person)  #SumVar / df
#1. the lower limit of the degrees of freedom of the mean-square statistic 
#is set at 1. qi is not allowed to be more than sqrt(2) for both. 
Item.out.qi<-sqrt(sum.Item.Out.Zstd.var)
Item.out.qi[Item.out.qi>sqrt(2)]<-sqrt(2)
Item.Out.ZSTD<-(Item.Outfit.mnsq^(1/3)-1)*(3/Item.out.qi)+(Item.out.qi/3)



#transform into data frames for IFILE & PFILE
Person.Outfit.mnsq<-as.data.frame(Person.Outfit.mnsq)
Person.Outfit.mnsq$entry<-row.names(Person.Outfit.mnsq)
Item.Outfit.mnsq<-as.data.frame(Item.Outfit.mnsq)
Item.Outfit.mnsq$item.names<-row.names(Item.Outfit.mnsq)
Per.Infit.mnsq<-as.data.frame(Per.Infit.mnsq)
Per.Infit.mnsq$entry<-row.names(Per.Infit.mnsq)
Item.Infit.mnsq<-as.data.frame(Item.Infit.mnsq)
Item.Infit.mnsq$item.names<-row.names(Item.Infit.mnsq)
Per.Out.ZSTD<-as.data.frame(Per.Out.ZSTD)
Per.Out.ZSTD$entry<-row.names(Per.Out.ZSTD)
Item.Out.ZSTD<-as.data.frame(Item.Out.ZSTD)
Item.Out.ZSTD$item.names<-row.names(Item.Out.ZSTD)
Per.In.ZSTD<-as.data.frame(Per.In.ZSTD)
Per.In.ZSTD$entry<-row.names(Per.In.ZSTD)
Item.In.ZSTD<-as.data.frame(Item.In.ZSTD)
Item.In.ZSTD$item.names<-row.names(Item.In.ZSTD)



####################################################################################
#                                                                                  #
#       Put EXTREMES back in before creating IFILE / PFILE / Summary Stats         #
#                                                                                  #
####################################################################################

#RMT Estimating Rasch (person, ability, theta) Measures with Known Dichotomous Item
#    Difficulties: Anchored Maximum Likelihood Estimation (AMLE)

extreme.matrix<-data[,item.start:((item.start-1)+n.items)]
#per.extreme<-extreme.matrix[-which(row.names(extreme.matrix) %in% row.names(data2)), ]
#itm.extreme<-extreme.matrix[, -which(names(extreme.matrix) %in% names(data2))]




extreme.person.raw<-apply(extreme.matrix, 1, sum, na.rm=TRUE)
names(extreme.person.raw)<-rownames(extreme.matrix)
extreme.person.admin<-apply(extreme.matrix, 1, function(x){sum(!is.na(x))})
names(extreme.person.admin)<-rownames(extreme.matrix)
extreme.person.prop<-extreme.person.raw/extreme.person.admin
names(extreme.person.prop)<-rownames(extreme.matrix)

extreme.item.raw<-apply(extreme.matrix, 2, sum, na.rm=TRUE)
extreme.item.admin<-apply(extreme.matrix, 2, function(x){sum(!is.na(x))})
extreme.item.prop<-extreme.item.raw/extreme.item.admin



#already have these from before...
extreme.person.ability<-log(extreme.person.prop/(1-extreme.person.prop))
extreme.item.difficulty<-log((1-extreme.item.prop)/extreme.item.prop)
extreme.item.difficulty<-ifelse(extreme.item.difficulty=="-Inf" | 
                                  extreme.item.difficulty=="Inf", 
                                NA, extreme.item.difficulty)
extreme.avg.item.difficulty<-mean(extreme.item.difficulty, na.rm=TRUE)
extreme.adj.item.difficulty<-extreme.item.difficulty-extreme.avg.item.difficulty


repeat { 
  expected.values<-function(x,y){(exp(x-y)/(1+exp(x-y) ) ) }
  extreme.Exp.Val<-outer(extreme.person.ability, extreme.adj.item.difficulty, expected.values)
  exp.var<-function(x){x*(1-x)}
  extreme.var<-exp.var(extreme.Exp.Val)
  
  extreme.person.var<-apply(extreme.var, 1, sum,na.rm=TRUE)
  extreme.person.var<-extreme.person.var*-1
  extreme.item.var<-apply(extreme.var, 2, sum, na.rm=TRUE)
  extreme.item.var<-extreme.item.var*-1
  
  extreme.res<-extreme.matrix-extreme.Exp.Val
  #res<-data2[,(1:ncol(data2))]-Exp.Val
  extreme.person.sum.of.residuals<-apply(extreme.res, 1, sum, na.rm=TRUE)
  extreme.item.sum.of.residuals<-apply(extreme.res, 2, sum, na.rm=TRUE)
  extreme.item.sum.of.residuals<-(-1*extreme.item.sum.of.residuals)
  
  extreme.per.square.resid<- sapply(extreme.person.sum.of.residuals, function(x) x^2)
  extreme.per.sum.of.square.resid<-sum(extreme.per.square.resid)
  extreme.item.square.resid<- sapply(extreme.item.sum.of.residuals, function(x) x^2)
  extreme.item.sum.of.square.resid<-sum(extreme.item.square.resid)
  
  extreme.new.ability<-extreme.person.ability-extreme.person.sum.of.residuals/
    extreme.person.var
  extreme.new.difficulty<-extreme.adj.item.difficulty-extreme.item.sum.of.residuals/
    extreme.item.var
  extreme.avg.new.difficulty<-mean(extreme.new.difficulty)
  extreme.adj.new.difficulty<-extreme.new.difficulty-extreme.avg.new.difficulty
  
  extreme.person.ability<-extreme.new.ability
  extreme.adj.item.difficulty<-extreme.adj.new.difficulty
  
  if(extreme.per.sum.of.square.resid<=converg) {
    break
  }
}


####################################################################################
#                                                                                  #
#                                   PFILE                                          #
#                                                                                  #
####################################################################################
data2<-data[,item.start:((item.start-1)+n.items)]


person.raw.score<-apply(data2, 1, sum, na.rm=TRUE)
person.raw.score<-data.frame(person.raw.score)
person.raw.score$entry<-as.numeric(row.names(person.raw.score))


person.tot.admin<-apply(data2, 1, function(x){sum(!is.na(x))})
person.tot.admin<-data.frame(person.tot.admin)
person.tot.admin$entry<-as.numeric(row.names(person.tot.admin))

person.pct<-person.raw.score$person.raw.score/person.tot.admin$person.tot.admin
person.pct<-data.frame(person.pct)
person.pct$entry<-as.numeric(row.names(person.pct))

person.meas<-data.frame(person.ability)
person.meas$entry<-as.numeric(row.names(person.meas))

#Person SE
P.SE<-Exp.Val
P.SE.plus<-P.SE*(1-P.SE)
P.SE.sum<-apply(P.SE.plus, 1, sum)
Person.Model.SE<-1/sqrt(P.SE.sum)
Person.Model.SE<-data.frame(Person.Model.SE)
Person.Model.SE$entry<-row.names(Person.Model.SE)

pfile<-merge(demographics, person.raw.score, by="entry", all=TRUE)
pfile<-merge(pfile, person.tot.admin, by="entry", all=TRUE)
pfile<-merge(pfile, person.pct, by="entry", all=TRUE)
pfile<-merge(pfile, person.meas, by="entry", all=TRUE)
pfile<-merge(pfile, Person.Model.SE, by="entry", all=TRUE)
pfile<-merge(pfile, Per.Infit.mnsq, by="entry", all=TRUE)
pfile<-merge(pfile, Person.Outfit.mnsq, by="entry", all=TRUE)
pfile<-merge(pfile, Per.In.ZSTD, by="entry", all=TRUE)
pfile<-merge(pfile, Per.Out.ZSTD, by="entry", all=TRUE)



####################################################################################
#                                                                                  #
#                                   IFILE                                          #
#                                                                                  #
####################################################################################

item.meas<-data.frame(adj.item.difficulty)
item.meas$item.names<-row.names(item.meas)

item.raw.score<-apply(data2, 2, sum, na.rm=TRUE)
item.raw.score<-data.frame(item.raw.score)
item.raw.score$item.names<-row.names(item.raw.score)

item.tot.admin<-apply(data2, 2, function(x){sum(!is.na(x))})
item.tot.admin<-data.frame(item.tot.admin)
item.tot.admin$item.names<-row.names(item.tot.admin)

item.pct<-item.raw.score$item.raw.score/item.tot.admin$item.tot.admin
item.pct<-data.frame(item.pct)
item.pct$entry<-row.names(item.pct)

item.names<-colnames(data2)
item.names<-data.frame(item.names)
item.names$entry<-as.numeric(row.names(item.names))

#Item SE
I.SE<-Exp.Val
I.SE.plus<-I.SE*(1-I.SE)
I.SE.sum<-apply(I.SE.plus, 2, sum)
Item.Model.SE<-1/sqrt(I.SE.sum)
Item.Model.SE<-data.frame(Item.Model.SE)
Item.Model.SE$item.names<-row.names(Item.Model.SE)



ifile<-merge(item.names, item.raw.score, by="item.names", all=TRUE)
ifile<-merge(ifile, item.tot.admin, by="item.names", all=TRUE)
ifile<-merge(ifile, item.pct, by="entry", all=TRUE)
ifile<-merge(ifile, item.meas, by="item.names", all=TRUE)
ifile<-merge(ifile, Item.Model.SE, by="item.names", all=TRUE)
ifile<-merge(ifile, Item.Infit.mnsq, by="item.names", all=TRUE)
ifile<-merge(ifile, Item.In.ZSTD, by="item.names", all=TRUE)
ifile<-merge(ifile, Item.Outfit.mnsq, by="item.names", all=TRUE)
ifile<-merge(ifile, Item.Out.ZSTD, by="item.names", all=TRUE)
ifile<-ifile[order(ifile$entry),]




####################################################################################
#                                                                                  #
#                        Summary Satistics - (EVERYONE)                            #
#                                                                                  #
####################################################################################



#######################PERSON SUMMARY STATISTICS TABLE############################# 

mean.person.ability<-mean(pfile$person.ability, na.rm=TRUE)
mean.person.score<-mean(pfile$person.raw.score, na.rm=TRUE)
mean.person.count<-mean(pfile$person.tot.admin, na.rm=TRUE)
mean.person.se<-mean(pfile$Person.Model.SE, na.rm=TRUE)
mean.person.in_mnsq<-mean(pfile$Per.Infit.mnsq, na.rm=TRUE)
mean.person.in_zstd<-mean(pfile$Per.In.ZSTD, na.rm=TRUE)
mean.person.out_mnsq<-mean(pfile$Person.Outfit.mnsq, na.rm=TRUE)
mean.person.out_zstd<-mean(pfile$Per.Out.ZSTD, na.rm=TRUE)

person.means<-data.frame(mean.person.ability, mean.person.score, mean.person.count,
                         mean.person.se,mean.person.in_mnsq, mean.person.in_zstd,
                         mean.person.out_mnsq, mean.person.out_zstd)

names(person.means)[names(person.means) == "mean.person.ability"] <- "MEASURE"
names(person.means)[names(person.means) == "mean.person.score"] <- "SCORE"
names(person.means)[names(person.means) == "mean.person.count"] <- "COUNT"
names(person.means)[names(person.means) == "mean.person.se"] <- "MODEL.SE"
names(person.means)[names(person.means) == "mean.person.in_mnsq"] <- "IN.MNSQ"
names(person.means)[names(person.means) == "mean.person.in_zstd"] <- "IN.ZSTD"
names(person.means)[names(person.means) == "mean.person.out_mnsq"] <- "OUT.MNSQ"
names(person.means)[names(person.means) == "mean.person.out_zstd"] <- "OUT.ZSTD"
                    

sd.person.ability<-sd(pfile$person.ability, na.rm=TRUE)
sd.person.score<-sd(pfile$person.raw.score, na.rm=TRUE)
sd.person.count<-sd(pfile$person.tot.admin, na.rm=TRUE)
sd.person.se<-sd(pfile$Person.Model.SE, na.rm=TRUE)
sd.person.in_mnsq<-sd(pfile$Per.Infit.mnsq, na.rm=TRUE)
sd.person.in_zstd<-sd(pfile$Per.In.ZSTD, na.rm=TRUE)
sd.person.out_mnsq<-sd(pfile$Person.Outfit.mnsq, na.rm=TRUE)
sd.person.out_zstd<-sd(pfile$Per.Out.ZSTD, na.rm=TRUE)

person.sds<-data.frame(sd.person.ability, sd.person.score, sd.person.count,
                         sd.person.se,sd.person.in_mnsq, sd.person.in_zstd,
                         sd.person.out_mnsq, sd.person.out_zstd)

names(person.sds)[names(person.sds) == "sd.person.ability"] <- "MEASURE"
names(person.sds)[names(person.sds) == "sd.person.score"] <- "SCORE"
names(person.sds)[names(person.sds) == "sd.person.count"] <- "COUNT"
names(person.sds)[names(person.sds) == "sd.person.se"] <- "MODEL.SE"
names(person.sds)[names(person.sds) == "sd.person.in_mnsq"] <- "IN.MNSQ"
names(person.sds)[names(person.sds) == "sd.person.in_zstd"] <- "IN.ZSTD"
names(person.sds)[names(person.sds) == "sd.person.out_mnsq"] <- "OUT.MNSQ"
names(person.sds)[names(person.sds) == "sd.person.out_zstd"] <- "OUT.ZSTD"

min.person.ability<-min(pfile$person.ability, na.rm=TRUE)
min.person.score<-min(pfile$person.raw.score, na.rm=TRUE)
min.person.count<-min(pfile$person.tot.admin, na.rm=TRUE)
min.person.se<-min(pfile$Person.Model.SE, na.rm=TRUE)
min.person.in_mnsq<-min(pfile$Per.Infit.mnsq, na.rm=TRUE)
min.person.in_zstd<-min(pfile$Per.In.ZSTD, na.rm=TRUE)
min.person.out_mnsq<-min(pfile$Person.Outfit.mnsq, na.rm=TRUE)
min.person.out_zstd<-min(pfile$Per.Out.ZSTD, na.rm=TRUE)

person.mins<-data.frame(min.person.ability, min.person.score, min.person.count,
                       min.person.se,min.person.in_mnsq, min.person.in_zstd,
                       min.person.out_mnsq, min.person.out_zstd)

names(person.mins)[names(person.mins) == "min.person.ability"] <- "MEASURE"
names(person.mins)[names(person.mins) == "min.person.score"] <- "SCORE"
names(person.mins)[names(person.mins) == "min.person.count"] <- "COUNT"
names(person.mins)[names(person.mins) == "min.person.se"] <- "MODEL.SE"
names(person.mins)[names(person.mins) == "min.person.in_mnsq"] <- "IN.MNSQ"
names(person.mins)[names(person.mins) == "min.person.in_zstd"] <- "IN.ZSTD"
names(person.mins)[names(person.mins) == "min.person.out_mnsq"] <- "OUT.MNSQ"
names(person.mins)[names(person.mins) == "min.person.out_zstd"] <- "OUT.ZSTD"

max.person.ability<-max(pfile$person.ability, na.rm=TRUE)
max.person.score<-max(pfile$person.raw.score, na.rm=TRUE)
max.person.count<-max(pfile$person.tot.admin, na.rm=TRUE)
max.person.se<-max(pfile$Person.Model.SE, na.rm=TRUE)
max.person.in_mnsq<-max(pfile$Per.Infit.mnsq, na.rm=TRUE)
max.person.in_zstd<-max(pfile$Per.In.ZSTD, na.rm=TRUE)
max.person.out_mnsq<-max(pfile$Person.Outfit.mnsq, na.rm=TRUE)
max.person.out_zstd<-max(pfile$Per.Out.ZSTD, na.rm=TRUE)

person.maxs<-data.frame(max.person.ability, max.person.score, max.person.count,
                       max.person.se,max.person.in_mnsq, max.person.in_zstd,
                       max.person.out_mnsq, max.person.out_zstd)

names(person.maxs)[names(person.maxs) == "max.person.ability"] <- "MEASURE"
names(person.maxs)[names(person.maxs) == "max.person.score"] <- "SCORE"
names(person.maxs)[names(person.maxs) == "max.person.count"] <- "COUNT"
names(person.maxs)[names(person.maxs) == "max.person.se"] <- "MODEL.SE"
names(person.maxs)[names(person.maxs) == "max.person.in_mnsq"] <- "IN.MNSQ"
names(person.maxs)[names(person.maxs) == "max.person.in_zstd"] <- "IN.ZSTD"
names(person.maxs)[names(person.maxs) == "max.person.out_mnsq"] <- "OUT.MNSQ"
names(person.maxs)[names(person.maxs) == "max.person.out_zstd"] <- "OUT.ZSTD"


Person.Summary<-rbind(person.means, person.sds, person.maxs, person.mins)
Person.Summary$STAT<-c("Mean", "SD", "Max", "Min")
Person.Summary<-Person.Summary[,c("STAT", "MEASURE", "SCORE", "COUNT", "MODEL.SE", "IN.MNSQ",
                                  "IN.ZSTD", "OUT.MNSQ", "OUT.ZSTD")]


#######################ITEM SUMMARY STATISTICS TABLE###############################


mean.item.difficulty<-mean(ifile$adj.item.difficulty, na.rm=TRUE)
mean.item.score<-mean(ifile$item.raw.score, na.rm=TRUE)
mean.item.count<-mean(ifile$item.tot.admin, na.rm=TRUE)
mean.item.se<-mean(ifile$Item.Model.SE, na.rm=TRUE)
mean.item.in_mnsq<-mean(ifile$Item.Infit.mnsq, na.rm=TRUE)
mean.item.in_zstd<-mean(ifile$Item.In.ZSTD, na.rm=TRUE)
mean.item.out_mnsq<-mean(ifile$Item.Outfit.mnsq, na.rm=TRUE)
mean.item.out_zstd<-mean(ifile$Item.Out.ZSTD, na.rm=TRUE)

item.means<-data.frame(mean.item.difficulty, mean.item.score, mean.item.count,
                         mean.item.se,mean.item.in_mnsq, mean.item.in_zstd,
                         mean.item.out_mnsq, mean.item.out_zstd)

names(item.means)[names(item.means) == "mean.item.difficulty"] <- "MEASURE"
names(item.means)[names(item.means) == "mean.item.score"] <- "SCORE"
names(item.means)[names(item.means) == "mean.item.count"] <- "COUNT"
names(item.means)[names(item.means) == "mean.item.se"] <- "MODEL.SE"
names(item.means)[names(item.means) == "mean.item.in_mnsq"] <- "IN.MNSQ"
names(item.means)[names(item.means) == "mean.item.in_zstd"] <- "IN.ZSTD"
names(item.means)[names(item.means) == "mean.item.out_mnsq"] <- "OUT.MNSQ"
names(item.means)[names(item.means) == "mean.item.out_zstd"] <- "OUT.ZSTD"

sd.item.difficulty<-sd(ifile$adj.item.difficulty, na.rm=TRUE)
sd.item.score<-sd(ifile$item.raw.score, na.rm=TRUE)
sd.item.count<-sd(ifile$item.tot.admin, na.rm=TRUE)
sd.item.se<-sd(ifile$Item.Model.SE, na.rm=TRUE)
sd.item.in_mnsq<-sd(ifile$Item.Infit.mnsq, na.rm=TRUE)
sd.item.in_zstd<-sd(ifile$Item.In.ZSTD, na.rm=TRUE)
sd.item.out_mnsq<-sd(ifile$Item.Outfit.mnsq, na.rm=TRUE)
sd.item.out_zstd<-sd(ifile$Item.Out.ZSTD, na.rm=TRUE)

item.sds<-data.frame(sd.item.difficulty, sd.item.score, sd.item.count,
                       sd.item.se,sd.item.in_mnsq, sd.item.in_zstd,
                       sd.item.out_mnsq, sd.item.out_zstd)

names(item.sds)[names(item.sds) == "sd.item.difficulty"] <- "MEASURE"
names(item.sds)[names(item.sds) == "sd.item.score"] <- "SCORE"
names(item.sds)[names(item.sds) == "sd.item.count"] <- "COUNT"
names(item.sds)[names(item.sds) == "sd.item.se"] <- "MODEL.SE"
names(item.sds)[names(item.sds) == "sd.item.in_mnsq"] <- "IN.MNSQ"
names(item.sds)[names(item.sds) == "sd.item.in_zstd"] <- "IN.ZSTD"
names(item.sds)[names(item.sds) == "sd.item.out_mnsq"] <- "OUT.MNSQ"
names(item.sds)[names(item.sds) == "sd.item.out_zstd"] <- "OUT.ZSTD"

min.item.difficulty<-min(ifile$adj.item.difficulty, na.rm=TRUE)
min.item.score<-min(ifile$item.raw.score, na.rm=TRUE)
min.item.count<-min(ifile$item.tot.admin, na.rm=TRUE)
min.item.se<-min(ifile$Item.Model.SE, na.rm=TRUE)
min.item.in_mnsq<-min(ifile$Item.Infit.mnsq, na.rm=TRUE)
min.item.in_zstd<-min(ifile$Item.In.ZSTD, na.rm=TRUE)
min.item.out_mnsq<-min(ifile$Item.Outfit.mnsq, na.rm=TRUE)
min.item.out_zstd<-min(ifile$Item.Out.ZSTD, na.rm=TRUE)

item.mins<-data.frame(min.item.difficulty, min.item.score, min.item.count,
                       min.item.se,min.item.in_mnsq, min.item.in_zstd,
                       min.item.out_mnsq, min.item.out_zstd)

names(item.mins)[names(item.mins) == "min.item.difficulty"] <- "MEASURE"
names(item.mins)[names(item.mins) == "min.item.score"] <- "SCORE"
names(item.mins)[names(item.mins) == "min.item.count"] <- "COUNT"
names(item.mins)[names(item.mins) == "min.item.se"] <- "MODEL.SE"
names(item.mins)[names(item.mins) == "min.item.in_mnsq"] <- "IN.MNSQ"
names(item.mins)[names(item.mins) == "min.item.in_zstd"] <- "IN.ZSTD"
names(item.mins)[names(item.mins) == "min.item.out_mnsq"] <- "OUT.MNSQ"
names(item.mins)[names(item.mins) == "min.item.out_zstd"] <- "OUT.ZSTD"

max.item.difficulty<-max(ifile$adj.item.difficulty, na.rm=TRUE)
max.item.score<-max(ifile$item.raw.score, na.rm=TRUE)
max.item.count<-max(ifile$item.tot.admin, na.rm=TRUE)
max.item.se<-max(ifile$Item.Model.SE, na.rm=TRUE)
max.item.in_mnsq<-max(ifile$Item.Infit.mnsq, na.rm=TRUE)
max.item.in_zstd<-max(ifile$Item.In.ZSTD, na.rm=TRUE)
max.item.out_mnsq<-max(ifile$Item.Outfit.mnsq, na.rm=TRUE)
max.item.out_zstd<-max(ifile$Item.Out.ZSTD, na.rm=TRUE)

item.maxs<-data.frame(max.item.difficulty, max.item.score, max.item.count,
                       max.item.se,max.item.in_mnsq, max.item.in_zstd,
                       max.item.out_mnsq, max.item.out_zstd)

names(item.maxs)[names(item.maxs) == "max.item.difficulty"] <- "MEASURE"
names(item.maxs)[names(item.maxs) == "max.item.score"] <- "SCORE"
names(item.maxs)[names(item.maxs) == "max.item.count"] <- "COUNT"
names(item.maxs)[names(item.maxs) == "max.item.se"] <- "MODEL.SE"
names(item.maxs)[names(item.maxs) == "max.item.in_mnsq"] <- "IN.MNSQ"
names(item.maxs)[names(item.maxs) == "max.item.in_zstd"] <- "IN.ZSTD"
names(item.maxs)[names(item.maxs) == "max.item.out_mnsq"] <- "OUT.MNSQ"
names(item.maxs)[names(item.maxs) == "max.item.out_zstd"] <- "OUT.ZSTD"



Item.Summary<-rbind(item.means, item.sds, item.maxs, item.mins)
Item.Summary$STAT<-c("Mean", "SD", "Max", "Min")
Item.Summary<-Item.Summary[,c("STAT", "MEASURE", "SCORE", "COUNT", "MODEL.SE", "IN.MNSQ",
                                  "IN.ZSTD", "OUT.MNSQ", "OUT.ZSTD")]


####################################################################################
#                                                                                  #
#                            Reliability Table - (EVERYONE)                        #
#                                                                                  #
####################################################################################


# Person Reliability


#True SD is slightly inflated, leading to inflated Sep and Rel.
## Maybe it's not since I haven't added back the extreme scores yet.  
##    Add extremes and then re-compute.

Model.RMSE.per<-sqrt((sum(pfile$Person.Model.SE^2, na.rm=TRUE))/
                      length(na.omit(pfile$Person.Model.SE)))
Model.True.SD.per<- sqrt(sd(pfile$person.ability, na.rm=TRUE)^2-(Model.RMSE.per^2)  )
Model.Separation.per<-Model.True.SD.per / Model.RMSE.per
Model.Reliability.per<- Model.Separation.per^2/(1+Model.Separation.per^2)
Model.Strata.per<-(4*Model.Separation.per+1)/3


Real.RMSE.per<-(sum(pfile$Person.Model.SE^2  * 
                             sqrt(ifelse(pfile$Per.Infit.mnsq>=1, pfile$Per.Infit.mnsq, 1)), 
                           na.rm=TRUE))/length(na.omit(pfile$Person.Model.SE))
Real.True.SD.per<- sqrt(sd(pfile$person.ability, na.rm=TRUE)^2-(Real.RMSE.per^2)  )
Real.Separation.per<-Real.True.SD.per / Real.RMSE.per
Real.Reliability.per<- Real.Separation.per^2/(1+Real.Separation.per^2)
Real.Strata.per<-(4*Real.Separation.per+1)/3





# Item Reliability

REAL
-RMSE
-True.SD
-Separation
-Reliability

MODEL
-RMSE
-True.SD
-Separation
-Reliability




  