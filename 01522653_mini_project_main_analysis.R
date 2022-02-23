#clean workspace
rm(list=ls())
dev.off()
#packages needed and citation for reference
require(dplyr)
citation(package = "dplyr")
require(tidyr)
citation(package = "tidyr")
library(ape)
citation(package = "ape")
library(phytools)
citation(package = "phytools")
library(caper)
citation(package = "caper")
library(ggplot2)
citation(package = "ggplot2")
library(mgcv)
citation(package = "mgcv")
library(tidymv)
citation(package = "tidymv")

#import data 
full_data_snakes<- read.csv("Downloads/snaketraits/Snake_HR_wSex.csv")
str(full_data_snakes)
##wrangle data to compare sex dimorphism across in mass and  SVL across latitudes
#select columns required
snakes_sub<- dplyr::select(full_data_snakes, Study_Species:SiteLatitude, Sex,SVL)
#take out NAs and species with data on only one sex
snakes<- snakes_sub%>%drop_na(SVL)%>%filter(Tree_taxon!="Crotalus_atrox", Tree_taxon!="Morelia_viridis", Tree_taxon!="Vipera_latastei")
summary(snakes)
#get absolute latitudes for sites
snakes$abs_latitude <- abs(snakes$SiteLatitude)
#get mean mass and svl for each sex (get latitudes too they wont change because depend on site of collection for each species)
means_snakes<- snakes%>%group_by(Tree_taxon, Sex)%>%summarise(mean_SVL=mean(SVL), mean_lat=mean(abs_latitude))
#wrangle to have one row per species with means for each sex 
Summary <- means_snakes %>% group_by(Tree_taxon)%>%pivot_wider(names_from =  Sex, values_from = c(mean_SVL, mean_lat))
Summary
summary(Summary)
#calculate dimorphism indices and mutate them onto dataset
Ratio <- Summary %>% mutate(SVL_dimorphism=mean_SVL_Male/mean_SVL_Female, mean_lat=mean_lat_Female)
Ratio
summary(Ratio)
#get final form of data set by selecting columns of interest 
snake_dimorphism_data<- Ratio%>%dplyr::select(Tree_taxon, SVL_dimorphism:mean_lat)
#basic plots to get general idea
plot(SVL_dimorphism~mean_lat, data=snake_dimorphism_data)

#check normality 
hist(snake_dimorphism_data$SVL_dimorphism)#not very good but we will see the residuals 
hist(snake_dimorphism_data$mean_lat)

#run basic linear model 
base_model <- glm(SVL_dimorphism ~ mean_lat, data = snake_dimorphism_data)
# Inspect our linear model
summary(base_model)
#check residuals and plots
plot(base_model)#not good: can see two species with more leverage but expected because not many samples at some latitudes 
plot(density(base_model$residuals))#not great
#plot base model
plot(SVL_dimorphism ~ mean_lat, data = snake_dimorphism_data)
abline(base_model)

# now test effect of evolution using phylogenetically controlled model
#import tree
squamate_tree<- read.tree("Downloads/snakephylo/squam_shl_names.tre")
str(squamate_tree)
#select species names in primary dataset
drop.species <- snake_dimorphism_data$Tree_taxon
#get matching species in tree
selected.tips <- unlist(sapply(drop.species,grep,squamate_tree$tip.label))
drop.species <- squamate_tree$tip.label[selected.tips]
drop.species
#make pruned tree with species of interest
pruned_tree<-drop.tip(squamate_tree,
                      setdiff(squamate_tree$tip.label,drop.species))
#take out tree node labels for comparison of tips after
pruned_tree$node.label<-NULL
#plot tree simple tree to visualise 
plotTree(pruned_tree,ftype="i")

#get dimorphism data and tree together 
#match data and tree
snake_dimorphism_data<-as.data.frame(snake_dimorphism_data[match(pruned_tree$tip.label, snake_dimorphism_data$Tree_taxon),])
snake_comp <- comparative.data(phy = pruned_tree, data = snake_dimorphism_data, names.col = "Tree_taxon")

#pgls to see if info on phylogeny has any effect
#look for lambda
phylosig(snake_comp$phy, snake_comp$data$SVL_dimorphism, method="lambda", test=FALSE, nsim=5000, se=NULL, start=NULL,
         control = list())
#make pgls model
snake_pgls <- pgls(SVL_dimorphism ~ mean_lat, data = snake_comp, lambda = "ML", bounds = list(lambda = c(0.85, 0.90)))
summary(snake_pgls)
plot(density(snake_pgls$residuals))#not good because of lack of data = drags even more than when don't consider phylogeny 
#plot model 
plot(SVL_dimorphism ~ mean_lat, data = snake_dimorphism_data, ylab = "Sexual size dimorphism index", xlab= "Absolute Latitude (degrees)", ylim=c(0.7,1.20), xlim=c(9,60))
abline(a = coef(snake_pgls)[1], b = coef(snake_pgls)[2])
#ggplot 
snake.predict <- cbind(snake_dimorphism_data, predict(snake_pgls))
plot1 <- ggplot(snake.predict, aes(mean_lat,SVL_dimorphism))
plot1 <- plot1 + geom_point()+
  labs(
    x = "Latitude (degrees)",
    y = "Sexual size dimorphism index")+ 
  xlim(9,70)+
  ylim(0.7,1.20)
plot1 <- plot1 + geom_line(aes(mean_lat, `predict(snake_pgls)`))
plot1 <- plot1 + theme_classic()+
  theme(legend.position = "None")
plot1

#plot tree with colour
#get ggtree for this version of R
install.packages("remotes")
remotes::install_github("YuLab-SMU/ggtree")
library(ggtree)
citation(package = "ggtree")

#simple tree
ggtree(pruned_tree) + geom_treescale()+ geom_tiplab(size=3)

library(treeio)
citation(package = "treeio")
library(tidytree)
citation(package = "tidytree")
library(ggplot2)
citation(package = "ggplot2")
library(TDbook)
citation(package = "TDbook")

#match trait and tree names
snake_dimorphism_data$Tree_taxon %in% pruned_tree$tip.label
#check no species missing 
index <- !(snake_dimorphism_data$Tree_taxon %in% pruned_tree$tip.label)
snake_dimorphism_data[index,]
#make dataset with dilorphism index 
dimorph_data <- as.data.frame(snake_dimorphism_data[,2])
rownames(dimorph_data) <- snake_dimorphism_data$Tree_taxon
colnames(dimorph_data) <- "Sexual Dimorphism"
#plot tree 
p3 <- ggtree(pruned_tree) +
  geom_tiplab(size=3) 
#add heatmap next to it 
p4 <-  gheatmap(p3,
                dimorph_data ,
                offset=0.2, low="blue", high="orange", colnames_position = "top", font.size=4)

p4

#citation of R
RStudio.Version()

