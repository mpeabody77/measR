

####################################################################################
#                                                                                  #
#                               Person - Item Map                                  #
#                                                                                  #
####################################################################################


require(plyr)
item.map.person<-subset(pfile, select=c(person.ability, PERSON))
item.map.person$TYPE<-"PERSON"
item.map.person <- rename(item.map.person, replace = c("person.ability" = "MEASURE", 
                                                       "PERSON" = "NAME"))
item.map.item<-subset(ifile, select=c(adj.item.difficulty, item.names))
item.map.item$TYPE<-"ITEM"
item.map.item <- rename(item.map.item, replace = c("adj.item.difficulty" = "MEASURE",
                                                   "item.names" = "NAME"))
item.map<-rbind(item.map.person, item.map.item)

#for smaller test use the dotplot
#for bigger tests use a density plot
#maybe print both and see what people prefer

require(ggplot2)
ggplot(item.map, aes(x=factor(TYPE), y=MEASURE, fill=TYPE))+
  geom_dotplot(binaxis = "y",
               binwidth=.3,
               position = "dodge",
               binpositions="all")+
  theme_classic()+
  theme(axis.title.x=element_blank())+
  theme(legend.position="none")

#this is for plotting text...too busy when labels overlap
#ggplot(item.map, aes(x=factor(TYPE), y=MEASURE, label=NAME))+
#  geom_label(position=position_jitter(height = .5, width = .5))



set.seed(1234)
x <- rnorm(n=10000, mean=500, sd=100)
y <- rnorm(n=200, mean=400, sd=110)
df.x<-as.data.frame(x)
df.y<-as.data.frame(y)
df.x$TYPE<-"PERSON"
df.y$TYPE<-"ITEM"
df.x<-rename(df.x, replace=c("x" = "MEASURE"))
df.y<-rename(df.y, replace=c("y" = "MEASURE"))
item.map2<-rbind(df.x, df.y)

ggplot(item.map2, aes(MEASURE, fill=TYPE))+
  geom_density(alpha=7/10)+
  theme_classic()+
  theme(axis.title.x=element_blank())+
  theme(legend.position =c(0.8, 0.8),legend.title=element_blank())
