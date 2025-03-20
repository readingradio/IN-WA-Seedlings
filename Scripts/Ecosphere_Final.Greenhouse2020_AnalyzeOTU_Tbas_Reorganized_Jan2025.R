library(vegan)
library(ape)
library(doBy)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
#library(tanagR)
library(car)
#library(ggpubr)
library(gridExtra)
library(grid)
library(gridBase)
library(lme4)
library(nlme)
library(ggrepel)
library(ggh4x)
library(patchwork)

# set working directory

	#setwd("/Users/will1809/OneDrive - purdue.edu/Dissertation/GMW.Dissertation.Analyses/CH5")
	
	isolate_table <- read.csv("Isolates/Greenhouse2020_IsolateData_Compiled.June2021.csv")
	
	isolate_table[!is.na(isolate_table$voucher),"Isolate"] %>% length
	isolate_table[!is.na(isolate_table$voucher),"Tissue"] %>% table
	isolate_table[!is.na(isolate_table$voucher) & isolate_table$Tissue=="Canker","Reps.numerator"] %>% sum
	sum(isolate_table$otu %>% unique %>% is.na %>% !.)
	sum(isolate_table[isolate_table$Tissue=="Root","otu"] %>% unique %>% is.na %>% !.)
	sum(isolate_table[isolate_table$Tissue=="Canker","otu"] %>% unique %>% is.na %>% !.)
	
	str(isolate_table)
	# 46 genera
	# 32 families
	# 15 orders
	# 6 classes
	# 2 phyla
	# 71 OTUs
	# 265 unique sequences (voucher -1 for NA)
	# 783 isolates WRONG
	# 286 morphospecies (among multiple batches of culturing)
	#source("Code/OTU.functions.R")
	
	# put all relative abundances of isolates on same scale (out of 10 for roots, out of 5 for canker tissue)
	table(isolate_table$Reps.denominator)
	isolate_table$weighted_abundance <- (isolate_table$Reps.numerator / isolate_table$Reps.denominator) * 5 * as.integer(as.factor(isolate_table$Tissue))
	
	isolate_table$otu <- as.character(isolate_table$otu)
	isolate_table$SampleID <- as.character(isolate_table$SampleID)


# read in experiment metadata

	design <- read.csv('Phenotypes/MorphoGH2020Expt1.csv')[,c(1:3,6)]
	str(design)
	design$SampleID <- paste(as.character(design$State), as.character(design$Block), as.character(design$Rep), sep="-")
	design$SampleID # 75 plants
	
	design$SampleID %>% unique() %>% setdiff(unique(isolate_table$SampleID))
	samplecodes_isolates_nodata <- isolate_table$SampleID %>% unique() %>% setdiff(unique(design$SampleID))
	
	# check for erroneously named plant IDs sample names in isolate data
	isolate_table[isolate_table$SampleID %in% samplecodes_isolates_nodata,]

	cankers <- isolate_table[which(isolate_table$Tissue == "Canker"),]
	roots   <- isolate_table[which(isolate_table$Tissue == "Root"),]
	
	# check design balance
	xtabs(~ Inoculation + Block, data=design[which(design$State == "WA"),])
	xtabs(~ State + Block, data=design[which(design$Inoculation == "Gm"),])
	
	# check coverage of isolates
	design$SampleID %>% as.factor %>% intersect(roots$SampleID) %>% length # 69 plants
	design$SampleID %>% as.factor %>% intersect(cankers$SampleID) %>% length # 68 plants

	
# make a releative abundance table

	#canker_otu_table <- convert_otu_table_tbas_meta(cankers, design)
	
	######## just inoculated
	#canker2 <- left_join(cankers, design, by="SampleID") %>% filter(Inoculation == "Gm")
	
	#canker_otu_table <- canker2[which(!is.na(canker2$otu) & !is.na(canker2$SampleID)), c("SampleID","otu","weighted_abundance")] %>% reshape2::dcast(SampleID ~ otu, fun.aggregate=function(x) min(c(sum(x),5)),value.var="weighted_abundance")
	##########
	
	######### 
	
	canker_otu_table <- cankers[which(!is.na(cankers$otu) & !is.na(cankers$SampleID)), c("SampleID","otu","weighted_abundance")] %>% reshape2::dcast(SampleID ~ otu, fun.aggregate=function(x) min(c(sum(x),5)),value.var="weighted_abundance")
##########
	
	
	rownames(canker_otu_table) <- canker_otu_table$SampleID
	
	
	c.empty.samples <- setdiff(design$SampleID, canker_otu_table$SampleID)
	c.empty.rows <- as.data.frame(matrix(0,nrow=length(c.empty.samples), ncol=dim(canker_otu_table)[2]))
	rownames(c.empty.rows) <- c.empty.samples
	colnames(c.empty.rows) <- colnames(canker_otu_table)
	
	canker_otu_table <- rbind(canker_otu_table[,-1], c.empty.rows[,-1])
	canker_otu_table <- canker_otu_table[design$SampleID,]
	
	treat.gm.cank <- rownames(canker_otu_table) %in% (design %>% filter(Inoculation=="Gm") %>% dplyr::select(SampleID) %>% unlist)
	
	#root_otu_table <- convert_otu_table_tbas_meta(roots, design)

	
		###just inoculated
	#root2 <- left_join(roots, design, by="SampleID") %>% filter(Inoculation == "Gm")
	
	#root_otu_table <- root2[which(!is.na(root2$otu) & !is.na(root2$SampleID)), c("SampleID","otu","weighted_abundance")] %>% reshape2::dcast(SampleID ~ otu, fun.aggregate=function(x) min(c(sum(x),50)),value.var="weighted_abundance")
	
	
	############
	
	
############ all
	
	
	root_otu_table <- roots[which(!is.na(roots$otu) & !is.na(roots$SampleID)), c("SampleID","otu","weighted_abundance")] %>% reshape2::dcast(SampleID ~ otu, fun.aggregate=function(x) min(c(sum(x),50)),value.var="weighted_abundance")

	###########
	
	rownames(root_otu_table) <- root_otu_table$SampleID
	r.empty.samples <- setdiff(design$SampleID, root_otu_table $SampleID)
	r.empty.rows <- as.data.frame(matrix(0,nrow=length(r.empty.samples), ncol=dim(root_otu_table)[2]))
	rownames(r.empty.rows) <- r.empty.samples
	colnames(r.empty.rows) <- colnames(root_otu_table)
	
	root_otu_table <- rbind(root_otu_table[,-1], r.empty.rows[,-1])
	root_otu_table <- root_otu_table[design$SampleID,]

	row.names(design) <- design$SampleID
	treat.gm.root <- rownames(root_otu_table) %in% (design %>% filter(Inoculation=="Gm") %>% dplyr::select(SampleID) %>% unlist)
	treat.gm.cank==	treat.gm.root
## diversity plot

design$State=="WA"
design$State=="IN"

r.diversity <- diversity(root_otu_table)
c.diversity <- diversity(canker_otu_table)

gmtag <- which(design$Inoculation=='Gm')

Anova(lm(r.diversity ~ design$State), type='III')
Anova(lm(c.diversity ~ design$State), type='III')

summary(lm(r.diversity ~ design$State), type='III')
summary(lm(c.diversity ~ design$State), type='III')

Anova(lmer(r.diversity ~ design$State + design$State:design$Inoculation + (1|as.factor(design$State):as.factor(design$Block))), type='III')
Anova(lmer(c.diversity ~ design$State + design$State:design$Inoculation + (1|as.factor(design$State):as.factor(design$Block))), type='III')

windows();DHARMa::simulateResiduals(lmer(r.diversity ~ design$State + design$State:design$Inoculation + (1|as.factor(design$State):as.factor(design$Block)))) %>% plot
windows();DHARMa::simulateResiduals(lmer(c.diversity ~ design$State + design$State:design$Inoculation + (1|as.factor(design$State):as.factor(design$Block)))) %>% plot

windows(); hist(r.diversity)
windows(); hist(c.diversity)

Anova(lm(r.diversity[gmtag] ~ design$State[gmtag]), type='III')
summary(lm(r.diversity[gmtag] ~ design$State[gmtag]), type='III')

Anova(lm(c.diversity[gmtag] ~ design$State[gmtag]), type='III')

#plot(lm(r.diversity ~ design$State))

design$treat_long <- paste(design$State, design$Inoculation)
design$treat_long %>% unique

Anova(lm(r.diversity ~ design$treat_long))
Anova(lm(c.diversity ~ design$treat_long))

design[24:75,]

Anova(lm(r.diversity[27:75] ~ design$Inoculation[27:75]))
Anova(lm(c.diversity[27:75] ~ design$Inoculation[27:75]))
summary(lm(c.diversity[27:75] ~ design$Inoculation[27:75]))

newdata <- data.frame(cbind(DR=r.diversity, CR=c.diversity, State=design$State, Treat=design$treat_long))

lm.d <- lm(DR ~ State + State:Treat, data=newdata)
anova(lm.d)

lm.c <- lm(CR ~ State + State:Treat, data=newdata)
anova(lm.c)
############### ACROSS ALL

r.shannon.wa <- root_otu_table[design$State=="WA",] %>% diversity
r.shannon.in <- root_otu_table[design$State=="IN",] %>% diversity

r.wa.mean<-	root_otu_table[design$State=="WA",] %>% diversity %>% mean
r.wa.sd<-	root_otu_table[design$State=="WA",] %>% diversity %>% sd
r.wa.se<-	root_otu_table[design$State=="WA",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))
	
r.in.mean<-		root_otu_table[design$State=="IN",] %>% diversity %>% mean
r.in.sd<-		root_otu_table[design$State=="IN",] %>% diversity %>% sd
r.in.se<-		root_otu_table[design$State=="IN",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))

c.shannon.wa <- canker_otu_table[design$State=="WA",] %>% diversity
c.shannon.in <- canker_otu_table[design$State=="IN",] %>% diversity

c.wa.mean<-	canker_otu_table[design$State=="WA",] %>% diversity %>% mean
c.wa.sd<-	canker_otu_table[design$State=="WA",] %>% diversity %>% sd
c.wa.se<-	canker_otu_table[design$State=="WA",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))
	
c.in.mean<-		canker_otu_table[design$State=="IN",] %>% diversity %>% mean
c.in.sd<-		canker_otu_table[design$State=="IN",] %>% diversity %>% sd
c.in.se<-		canker_otu_table[design$State=="IN",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))

rc.se <- c(r.in.se, r.wa.se, c.in.se, c.wa.se)

toplot <- c(r.in.mean, r.wa.mean, c.in.mean, c.wa.mean)
names(toplot) <- c("IN", "WA", "IN", "WA")

#####################
### FIG 3 PANEL A ###
#####################

windows(width=4, height=5); par(xpd=T, oma=c(1,1,1,1), mar=c(0,4,0,0))
plot(toplot, pch=15, srt=90, ylim=c(0,2), xlim=c(.5,4.5), ylab="Shannon diversity", xlab="", axes=F)
axis(1, labels=c("Roots","Cankers"), at=c(1.5,3.5), las=0, lty=0, line=-2)
axis(2, las=2)
#axis(1, labels=names(toplot), at=c(1:4), lty=0, las=0)
arrows(x0=c(1:4), y0=toplot-rc.se, y1=toplot+rc.se, length=.06, angle=90, code=3)
text(x=c(1:4), y=toplot+rc.se+.2, labels=names(toplot))
text(x=c(1.5,3.5), y=c(2,2), labels=c("p = 0.02", "p = 0.69"))



###################
####################
###################
##################

## X ALL INOCULATION-STATE COMBOS

r.shannon.wa.gm <- root_otu_table[design$treat_long=="WA Gm",] %>% diversity
r.shannon.control.gm <- root_otu_table[design$treat_long=="WA Sham",] %>% diversity
r.shannon.in <- root_otu_table[design$treat_long=="IN Gm",] %>% diversity

r.wa.gm.mean<-	root_otu_table[design$treat_long=="WA Gm",] %>% diversity %>% mean
r.wa.gm.sd<-	root_otu_table[design$treat_long=="WA Gm",] %>% diversity %>% sd
r.wa.gm.se<-	root_otu_table[design$treat_long=="WA Gm",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))

r.in.mean<-		root_otu_table[design$treat_long=="IN Gm",] %>% diversity %>% mean
r.in.sd<-		root_otu_table[design$treat_long=="IN Gm",] %>% diversity %>% sd
r.in.se<-		root_otu_table[design$treat_long=="IN Gm",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))

r.wa.control.mean<-	root_otu_table[design$treat_long=="WA Sham",] %>% diversity %>% mean
r.wa.control.sd<-	root_otu_table[design$treat_long=="WA Sham",] %>% diversity %>% sd
r.wa.control.se<-	root_otu_table[design$treat_long=="WA Sham",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))

c.shannon.wa.gm <- canker_otu_table[design$treat_long=="WA Gm",] %>% diversity
c.shannon.wa.control <- canker_otu_table[design$treat_long=="WA Sham",] %>% diversity
c.shannon.in <- canker_otu_table[design$treat_long=="IN Gm",] %>% diversity

c.wa.gm.mean<-	canker_otu_table[design$treat_long=="WA Gm",] %>% diversity %>% mean
c.wa.gm.sd<-	canker_otu_table[design$treat_long=="WA Gm",] %>% diversity %>% sd
c.wa.gm.se<-	canker_otu_table[design$treat_long=="WA Gm",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))

c.wa.control.mean<-	canker_otu_table[design$treat_long=="WA Sham",] %>% diversity %>% mean
c.wa.control.sd<-	canker_otu_table[design$treat_long=="WA Sham",] %>% diversity %>% sd
c.wa.control.se<-	canker_otu_table[design$treat_long=="WA Sham",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))

c.in.mean<-		canker_otu_table[design$treat_long=="IN Gm",] %>% diversity %>% mean
c.in.sd<-		canker_otu_table[design$treat_long=="IN Gm",] %>% diversity %>% sd
c.in.se<-		canker_otu_table[design$treat_long=="IN Gm",] %>% diversity %>% (function(x) sd(x)/sqrt(length(x)))

rc.se2 <- c(r.in.se, r.wa.gm.se, r.wa.control.se, c.in.se, c.wa.gm.se, c.wa.control.se)

toplot2 <- c(r.in.mean, r.wa.gm.mean, r.wa.control.mean, c.in.mean, c.wa.gm.mean, c.wa.control.mean)
toplot2.above <- c("IN", "WA", "IN", "WA")
toplot2.below <- c("Gm", "Gm", "Sham", "Gm", "Gm", "Sham")
### FIG 3 PANEL A

windows(width=5, height=5); par(xpd=T, oma=c(1,1,1,1), mar=c(0,4,0,0))
plot(toplot2,
     pch=15,
     srt=90,
     ylim=c(0,2),
     xlim=c(.5,6.5),
     ylab="Shannon diversity", xlab="", axes=F)

axis(1, labels=c("Roots","Cankers"), at=c(2,5), las=0, lty=0, line=-2)
axis(2, las=2)
#axis(1, labels=names(toplot), at=c(1:4), lty=0, las=0)
arrows(x0=c(1:6), y0=toplot2-rc.se2, y1=toplot2+rc.se2, length=.06, angle=90, code=3)
text(x=c(1,2.5,4,5.5), y=toplot2[c(1,2,4,5)]+rc.se2[c(1,2,4,5)]+.3, labels=toplot2.above, cex=.75)
text(x=c(1:6), y=toplot2+rc.se2+.1, labels=toplot2.below, cex=.75)

text(x=c(2,5), y=c(2,1.25), labels=c("p = 0.02", "p = 0.68"), cex=.75)
text(x= 5.5, y=toplot2[5]+rc.se2[5]+.2, labels="p = 0.06", cex=.75)


###################
#####################
###################

####################
# stacked barplots #
####################


	#set up color pallette

######
# from the internets
	get_legend <- function(p) {
		tmp <- ggplot_gtable(ggplot_build(p))
		leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
		legend <- tmp$grobs[[leg]]
		legend
	}

	library(scales)

#	colorss<-c(tanagr_palette("buthraupis_montana")[c(1,4,5)], # 3
#			tanagr_palette("tangara_velia")[c(1,3,4,5)],			#+4=7
#			tanagr_palette("stilpnia_preciosa")[c(1:5)],			#+5=12
#			tanagr_palette("chlorornis_riefferii")[c(1:2)],		#+2=14
#			tanagr_palette("ramphocelus_sanguinolentus")[c(4:5)],#16
#			tanagr_palette("dacnis_berlepschi")[3:4],				#+2=18
#			tanagr_palette("cyanerpes_cyaneus")[2:3],				#+2=20
#			tanagr_palette("bangsia_edwardsi")[c(1,3)],			#+2=22
#			tanagr_palette("tangara_chilensis")[c(1,5)])

#	tanagr_select <- scale_fill_manual(values=colorss)

	colorscheme1<-	 c("#ACD000",
							"#8B5B06",
							"#345896",
							"#4D820C",
							"#F1ED7F",
							"#3A170E",
							"#004D6B",
							"#302F35",
							"#E93924",
							"#73DBDA",
							"#D4940E",
							"#AEC7E0",
							"#020104")
	
	colorss<-colorscheme1
	tanagr_select <- scale_fill_manual(values=colorss)
	
# themes

stackplottheme<- theme(
				axis.text.x=element_text(size=14, angle = 45, hjust = 1, color="black"),
				axis.text.y=element_text(size = 14,color="black"),
				axis.line.y=element_line(color = "black"),
				axis.title=element_text(size=24),
				axis.line.x = element_blank(),
				axis.ticks.x = element_blank(),
				panel.border = element_blank(),
				panel.background = element_rect(fill = "transparent"),
    			plot.background = element_rect(fill = "transparent", color = NA),
    			panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
    			legend.background = element_rect(fill = "transparent"),
    			legend.box.background = element_blank(),
    			legend.text = element_text(size=9),
    			legend.spacing.x = unit(.75, 'lines'),
    			legend.spacing.y = unit(1, 'lines'),
    			legend.title = element_text(size=10),
    			legend.key.size = unit(.75, "lines"),
    			legend.key.height = unit(1, "lines"),
    			legend.key.width = unit(1, "lines"))

### Roots

	roots.t <- transpose(as.data.frame(root_otu_table)) %>% as.data.frame
	rownames(roots.t)<- colnames(root_otu_table)
	colnames(roots.t)<- rownames(root_otu_table)
	roots.t$otu = rownames(roots.t)
	class(roots.t$otu)

	roots$otu <- as.character(roots$otu)
	roots$Phylum <- as.character(roots$Phylum)
	roots$Class <-  as.character(roots$Class)
	roots$Order <-  as.character(roots$Order)
	roots$Family <- as.character(roots$Family)
	roots$Genus <-  as.character(roots$Genus)

	roots$Taxon.assignment <-  as.character(roots$Taxon.assignment)
	otu_translation_table <-   NULL
	for (i in rownames(roots.t))
		otu_translation_table <- rbind(
			otu_translation_table,
			roots[which(roots$otu==i)[1],
			c("otu","Phylum","Class","Order","Family","Genus","Taxon.assignment")])

	roots.joined <- left_join(otu_translation_table, roots.t, by = "otu")

	design$StateBlockInoc <- paste(design$State, design$Block, design$Inoc, sep="") 
	design$StateInoc      <- paste(design$State, design$Inoc, sep="") 
	design$StateBlock     <- paste(design$State, design$Block, sep="") 


### Cankers

	canks.t <- transpose(as.data.frame(canker_otu_table)) %>% as.data.frame
	rownames(canks.t)<- colnames(canker_otu_table)
	colnames(canks.t)<- rownames(canker_otu_table)
	canks.t$otu = rownames(canks.t)
	class(canks.t$otu)

	cankers$otu <- as.character(cankers$otu)
	cankers$Phylum <- as.character(cankers$Phylum)
	cankers$Class <-  as.character(cankers$Class)
	cankers$Order <-  as.character(cankers$Order)
	cankers$Family <- as.character(cankers$Family)
	cankers$Genus <-  as.character(cankers$Genus)

	cankers$Taxon.assignment <-  as.character(cankers$Taxon.assignment)
	otu_translation_table <-   NULL
	for (i in rownames(canks.t))
		otu_translation_table <- rbind(
			otu_translation_table,
			cankers[which(cankers$otu==i)[1],
			c("otu","Phylum","Class","Order","Family","Genus","Taxon.assignment")])

	canks.joined <- left_join(otu_translation_table, canks.t, by = "otu")

	design$StateBlockInoc <- paste(design$State, design$Block, design$Inoc, sep="") 
	design$StateInoc      <- paste(design$State, design$Inoc, sep="") 
	design$StateBlock     <- paste(design$State, design$Block, sep="") 


### Roots

	# organize by genus
	root.genus <- roots.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>%
		group_by(State, Block, Inoculation, Genus) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(State, Block, Inoculation, desc(count))

	# organize by family
	root.family <- roots.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>% group_by(State, Block, Inoculation, Family) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(State, Block, Inoculation, desc(count))

	# organize by order
	root.order <- roots.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>%
		group_by(State, Block, Inoculation, Order) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(State, Block, Inoculation, desc(count))

	# organize by class
	root.class <- roots.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>%
		group_by(State, Block, Inoculation, Class) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(State, Block, Inoculation, desc(count))
### Cankers

	# organize by genus
	cank.genus <- canks.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>% group_by(State, Block, Inoculation, Genus) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(State, Block, Inoculation, desc(count))

	# organize by family
	cank.family <- canks.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>% group_by(State, Block, Inoculation, Family) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(State, Block, Inoculation, desc(count))

	# organize by order
	cank.order <- canks.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>%
		group_by(State, Block, Inoculation, Order) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(State, Block, Inoculation, desc(count))

	# organize by class
	cank.class <- canks.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>%
		group_by(State, Block, Inoculation, Class) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(State, Block, Inoculation, desc(count))


###########################
# plot by state and order #
###########################

	r.plot.state.order.data <- root.order %>%
	  mutate(StateInoculation = paste(State, Inoculation)) %>%
		group_by (StateInoculation, Order) %>%
		dplyr::summarize (state.count = sum(count)) %>%
		left_join(., root.order %>%
		            mutate(StateInoculation = paste(State, Inoculation)) %>%
			group_by (StateInoculation) %>%
			dplyr::summarize (t = sum(count))) %>%
		dplyr::mutate(relabund= state.count/t)
	
	r.plot.state.order.data$Tissue <- "Root"

	c.plot.state.order.data <- cank.order %>%
	  mutate(StateInoculation = paste(State, Inoculation)) %>%
	  group_by (StateInoculation, Order) %>%
		dplyr::summarize (state.count = sum(count)) %>%
		left_join(., cank.order %>%
		            mutate(StateInoculation = paste(State, Inoculation)) %>%
		            group_by (StateInoculation) %>%
			dplyr::summarize (t = sum(count))) %>%
		dplyr::mutate(relabund= state.count/t)
	
	c.plot.state.order.data$Tissue <- "Canker"

	all.plot.state.order.data <- rbind(r.plot.state.order.data, c.plot.state.order.data)
	all.plot.state.order.data$Tissue <- factor(all.plot.state.order.data$Tissue, levels=c("Root","Canker"))

	n.order <- length(unique(all.plot.state.order.data$Order))

	all.plot.state.order <- all.plot.state.order.data %>%
	ggplot(aes(x=StateInoculation, y=relabund, fill=Order)) +
		geom_bar(aes(), stat="identity", position="fill", size=0) +
		facet_wrap(~Tissue, strip.position="bottom")+
		scale_x_discrete(position = "top")+
		xlab("") +
  		ylab("") +
  		ggtitle("B")+
  		stackplottheme+
  		theme(
  			axis.text.x=element_text(size=14, angle = 45, hjust = 0, color="black"),
  			strip.background = element_rect(color="white", size=0, fill="white"),
  			strip.text =element_text(size=13))

	
	all.plot.state.order2 <- all.plot.state.order.data %>%
	  ggplot(aes(x=StateInoculation, y=relabund, fill=Order)) +
	  geom_bar(aes(), stat="identity", position="fill", size=0) +
	  facet_wrap(~Tissue)+#, strip.position="bottom")+
	  #scale_x_discrete(position = "top")+
	  xlab("") +
	  ylab("") +
	#  ggtitle("B")+
	  stackplottheme+
	  theme(
	    axis.text.x=element_text(angle=45, hjust=1, size=14, vjust = 1),
	    strip.background = element_rect(color="white", size=0, fill="white"),
	    strip.text = element_text(vjust=1,size=15))
	
	#windows(width=6,height=4)
	#all.plot.state.order+ scale_fill_manual(values=colorscheme1)
	all.plot.state.order+ scale_fill_manual(values= colorss[1:n.order])
# composite

	windows(width=8.5, height=5); layout(mat=rbind(c(1,1,1,2,2,2,2,2),c(1,1,1,2,2,2,2,2))); par(xpd=T, oma=c(0,0,0,0), mar=c(3,4,5,0), cex=1.1)
	plot(toplot, pch=15, srt=90, ylim=c(0,2), xlim=c(.5,4.5), ylab="Shannon diversity", xlab="", axes=F)
	axis(1, labels=c("Root","Canker"), at=c(1.5,3.5), las=0, lty=0, line=-1.4)
	axis(2, las=2)
	arrows(x0=c(1:4), y0=toplot-rc.se, y1=toplot+rc.se, length=.06, angle=90, code=3)
	text(x=c(1:4), y=toplot+rc.se+.2, labels=names(toplot))
	text(x=c(1.5,3.5), y=c(2,2), labels=c("p = 0.02", "p = 0.68"), cex=1)
	text(0.49,2.69, labels="A")

	plot.new()
	vps <- baseViewports()
	pushViewport(vps$figure)
	vp1<-plotViewport(c(2,0,.25,0))
	print(all.plot.state.order+ scale_fill_manual(values= colorss[1:n.order]), vp=vp1)

	######################################################
	#### plot by state, block, inoculation and family ####
	######################################################
	
	newdata <- data.frame(cbind(DR=r.diversity, CR=c.diversity, State=as.factor(design$State), Treat=as.factor(design$treat_long), Block=as.factor(design$StateBlock)))
	
	#lm.d <- lm(DR ~ State + State:Treat, data=newdata)
	#s.d <- anova(lm.d)[4:5]
	#s.d
	
	lmer.d <- lmer(DR ~ State + State:Treat + (1|State:Block), data=newdata)
	s.d <- as.data.frame(Anova(lmer.d))
	s.d$p <- as.character(paste("p =", round(s.d$`Pr(>Chisq)`,2)))
	s.d

	#lm.c <- lm(CR ~ State + State:Treat, data=newdata)
	lmer.c <- lmer(CR ~ State + State:Treat + (1|State:Block), data=newdata)
	s.c <- as.data.frame(Anova(lmer.c))
	s.c$p <- as.character(paste("p =", round(s.c$`Pr(>Chisq)`,2)))
	s.c
	
	toplot2_df <-
	  data.frame(
	    mean = toplot2,
	    se   = rc.se2,
	    tissuelabs = c(rep("Root",3),rep("Canker",3)) %>% factor(levels = c("Roots","Cankers")),
	    Inoc    = rep(c(rep("Gm",2),"Sham"),2),
	    statelab   = rep(c("IN",rep("WA",2)),2)
	  ) %>%
	  
	  mutate(treatlong = paste(statelab, Inoc))
	
	newdata_gg <- data.frame(cbind('Root'=r.diversity,
	                               'Canker'=c.diversity,
	                               statelab=design$State,
	                               Inoc=design$Inoculation)) %>%
	  reshape2::melt(variable.name = 'tissuelabs',
	                 id.vars = c("statelab","Inoc")) %>%
	  mutate(value = as.numeric(value)) %>%
	  mutate(treatlong = paste(statelab, Inoc))
	
	windows(5,5)
	
#	ggplot()+
	library(ggpubr)
	
	levels(newdata_gg$tissuelabs) <- paste(levels(newdata_gg$tissuelabs), "\n", sep="")
	
	ggbp <-
	  
	  ggboxplot(
	    data = newdata_gg,
	    x = 'treatlong',
	    y = 'value',
	    color = 'Inoc')+
	  stat_compare_means(
	    mapping = aes(label = paste0("p = ", after_stat(p.format))),
	    method = 't.test',
	    comparisons = list(statelab=c("IN Gm", "WA Gm"), c("WA Gm", "WA Sham")),
	  ) +	
	  facet_wrap('tissuelabs') + theme_minimal() +
	    theme(axis.text.x = element_text(angle=45, hjust=1, size=14, vjust = 1),
	          strip.text = element_text(vjust=1,size=15),
	          axis.text.y = element_text(size=15),
	          legend.text = element_text(size=12),
	          legend.title =element_blank(),
	          axis.title   =element_text(size=15),
	          legend.position = 'right')+
	  labs(y = "\nShannon diversity\n", x = "\nLocation & Inoculation") +
	  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5,3), limits=c(0,2.75))#+
#	  geom_point(data= toplot2_df,
#	             aes(x = treatlong,
#	                 y = mean,
#	                 color = Inoculation), shape = 15, size = 1.5, inherit.aes=FALSE) +
#	  geom_linerange(data= toplot2_df,
#	                 aes(x = treatlong,
#	                     ymin = mean-se,
#	                     ymax = mean+se,
#	                     color = Inoculation), inherit.aes=FALSE)#, position=position_dodge(width = 1)) 
	
	
	
	########## COMPOSITE NEW
	windows(width=6.5, height=6.5)#; layout(mat=rbind(c(1,1,1,2,2,2,2,2),c(1,1,1,2,2,2,2,2))); par(xpd=T, oma=c(0,0,0,0), mar=c(3,4,5,0), cex=1.1)
	
	### LEFT PANEL
	library(devEMF)
	library(Cairo)
	setwd('Revision_Jan_2025')
	#emf(file = "Figures_ordered/Fig_2AB_new_correct_dimensions.emf",12,5,emfPlus=T)
	svg(file = "Fig_2AB_vertical_march2025.svg",6.5,6.5)
	(ggbp + theme(axis.text.x = element_blank(),
	              axis.line.y = element_line(),
	              axis.ticks.y = element_line()) +
	    coord_cartesian(clip = "off")) +
	  (all.plot.state.order2+
	         labs(y = "\nRelative abundance\n",x = "\nLocation & Inoculation")+
	         scale_fill_manual(values= colorss[1:n.order])+
	         theme(axis.title.y = element_text(size = 15),
	               axis.title   = element_text(size=15),
	               strip.text = element_blank())) +
	  plot_layout(ncol = 1, nrow = 2, axes = "collect_x")+#, widths = c(2,3))+
	  #egg::tag_facet()
	  plot_annotation(tag_levels = 'A')
	dev.off()
#	plot(toplot2,
#	     pch=15,
#	     srt=90,
#	     ylim=c(0,2),
#	     xlim=c(.5,6.5),
#	     ylab="Shannon diversity", xlab="", axes=F)
#	axis(1, labels=c("Roots","Cankers"), at=c(2,5), las=0, lty=0, line=-2)
#	axis(2, las=2)
#	#axis(1, labels=names(toplot), at=c(1:4), lty=0, las=0)
#	arrows(x0=c(1:6), y0=toplot2-rc.se2, y1=toplot2+rc.se2, length=.06, angle=90, code=3)
#	# pvals top
#	text(x=c(1.75,4.75), y=c(2,2), labels=c(s.d$p[1], s.c$p[1]), cex=.75)
#	# IN WA
#	text(x=c(1,2.5,4,5.5), y=
	       #toplot2[c(1,2,4,5)]+rc.se2[c(1,2,4,5)]+.3,
#	       rep(1.5+2*.5/3,2), labels=toplot2.above, cex=.75)
	# pvals inoc
#	text(x= c(2.5,5.5), y=
	       #toplot2[5]+rc.se2[5]+.2,
#	       rep(1.5+.5/3,2), labels=c(s.d$p[2], s.c$p[2]), cex=.75)
	#text(x=c(1:6), y=toplot2+rc.se2+.1, labels=c(), cex=.75)
	# Gm Sham
#	text(x=c(1:6), y=toplot2+rc.se2+.1, labels=toplot2.below, cex=.75)
	# A
#	text(0.49,2.69, labels="A")
	
	### RIGHT PANEL
#	plot.new()
#	vps <- baseViewports()
#	pushViewport(vps$figure)
#	vp1<-plotViewport(c(2,0,.25,0))
#	print(all.plot.state.order+ scale_fill_manual(values= colorss[1:n.order]), vp=vp1)

######################################################
#### FIG 5A													    	#
#### plot by state, block, inoculation and family ####
######################################################

	root.fam.2 <- roots.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		left_join(design, by="SampleID") %>%
		group_by(StateBlockInoc, Family) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(StateBlockInoc, desc(count))

	r.plot.state.fam.data <- root.fam.2 %>%
		left_join(., root.fam.2 %>%
				group_by (StateBlockInoc) %>%
				dplyr::summarize (t = sum(count))) %>%
				dplyr::mutate(relabund= count/t)
				
	r.plot.state.fam.data$StateBlockInoc <- factor(r.plot.state.fam.data$StateBlockInoc, levels=unique(r.plot.state.fam.data$StateBlockInoc)[c(1,2,4,3,6,5,8,7)])

	levels(r.plot.state.fam.data$StateBlockInoc) <-
	  c("IN P1 (Gm)", 		"IN P2 (Gm)",
	    "WA P1 (Sham)", 	"WA P1 (Gm)",
	    "WA P2 (Sham)", 	"WA P2 (Gm)",
	    "WA P3 (Sham)", 	"WA P3 (Gm)")

	n.fams <- length(unique(r.plot.state.fam.data$Family))

	r.plot.state.fam <- r.plot.state.fam.data %>%
		ggplot(aes(x=StateBlockInoc, y=relabund, fill=Family)) +
			geom_bar(aes(), stat="identity", position="fill", size=0)+
  			xlab("") +
  			ylab("") + stackplottheme

####################################################
# FIG 5B                                           #
# plot gm inoculated by necrotic area quantile X 8 #
####################################################


# canker necrosis area and composition

	necrosis <- read.csv ("Phenotypes/WA-IN-Transplant2019-20.Necrosis.csv")
	necrosis$SampleID <- paste(necrosis$State, necrosis$Block, necrosis$Plant, sep="-")
	necrosis <- left_join(necrosis, design[,4:5], by="SampleID")
	necrosis <- necrosis[necrosis$Inoculation == "Gm",]

# adonis based on quantiles

	necr.summary <- summaryBy(Area~SampleID, data=necrosis, FUN=c(mean,sd), id=c("State","Inoculation"))
	necr.summary <-necr.summary[order(necr.summary$Area.mean),]
	#dim(necr.summary)
	necr.summary$Area.8.bins <- c(rep(1,7),rep(2,6),rep(3,7),rep(4,6),rep(5,7),rep(6,6),rep(7,7),rep(8,6))

	root.quantile8.plot.fam <- roots.joined %>%
		reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>%
		right_join(necr.summary, by="SampleID") %>%
		group_by(Area.8.bins, Family) %>%
		dplyr::summarize (count = sum(Number)) %>%
		arrange(Area.8.bins, desc(count))

	root.quantile8.plot.data.fam <- root.quantile8.plot.fam %>%
		left_join(., root.quantile8.plot.fam %>%
			group_by (Area.8.bins) %>%
			dplyr::summarize (t = sum(count))) %>%
			dplyr::mutate(relabund= count/t)

	necr.summary$Area.mean <- round(necr.summary$Area.mean,0)

	root.quantile8.plot.data.fam$Area.8.bins <- factor(root.quantile8.plot.data.fam$Area.8.bins)
	levels(root.quantile8.plot.data.fam$Area.8.bins) <- c(
		with(necr.summary[necr.summary$Area.8.bins == 1, ], paste(min(Area.mean), "-", max(Area.mean), "mm²")),
		with(necr.summary[necr.summary$Area.8.bins == 2, ], paste(min(Area.mean), "-", max(Area.mean), "mm²")),
		with(necr.summary[necr.summary$Area.8.bins == 3, ], paste(min(Area.mean), "-", max(Area.mean), "mm²")),
		with(necr.summary[necr.summary$Area.8.bins == 4, ], paste(min(Area.mean), "-", max(Area.mean), "mm²")),
		with(necr.summary[necr.summary$Area.8.bins == 5, ], paste(min(Area.mean), "-", max(Area.mean), "mm²")),
		with(necr.summary[necr.summary$Area.8.bins == 6, ], paste(min(Area.mean), "-", max(Area.mean), "mm²")),
		with(necr.summary[necr.summary$Area.8.bins == 7, ], paste(min(Area.mean), "-", max(Area.mean), "mm²")),
		with(necr.summary[necr.summary$Area.8.bins == 8, ], paste(min(Area.mean), "-", max(Area.mean), "mm²")))
	n.fams <- length(unique(root.quantile8.plot.data.fam$Family))

	library(tanagR)
	
		colorss<-c(tanagr_palette("buthraupis_montana")[c(1,4,5)], # 3
				tanagr_palette("tangara_velia")[c(1,3,4,5)],			#+4=7
				tanagr_palette("stilpnia_preciosa")[c(1:5)],			#+5=12
				tanagr_palette("chlorornis_riefferii")[c(1:2)],		#+2=14
				tanagr_palette("ramphocelus_sanguinolentus")[c(4:5)],#16
				tanagr_palette("dacnis_berlepschi")[3:4],				#+2=18
				tanagr_palette("cyanerpes_cyaneus")[2:3],				#+2=20
				tanagr_palette("bangsia_edwardsi")[c(1,3)],			#+2=22
				tanagr_palette("tangara_chilensis")[c(1,5)])
	
	#	tanagr_select <- scale_fill_manual(values=colorss)
	
	colorscheme1<-	 c("#ACD000",
	                  "#8B5B06",
	                  "#345896",
	                  "#4D820C",
	                  "#F1ED7F",
	                  "#3A170E",
	                  "#004D6B",
	                  "#302F35",
	                  "#E93924",
	                  "#73DBDA",
	                  "#D4940E",
	                  "#AEC7E0",
	                  "#020104")
	
	#colorss<-colorscheme1
	tanagr_select <- scale_fill_manual(values=colorss)
	
r.plot.8.fam <- root.quantile8.plot.data.fam %>%
		ggplot(aes(x=Area.8.bins, y=relabund, fill=Family)) +
			geom_bar(aes(), stat="identity", position="fill", size=0)+
  			xlab("") +
  			ylab("") +	stackplottheme

####################################################
# COMPOSITE FIG 5A + 5B                            #
####################################################

p1 <- r.plot.state.fam + ggtitle("A")+ #labs(tag="A")+ 
  tanagr_select +
  theme(#legend.position = "none",
        axis.text.x=element_text(size = 10, angle = 45, hjust = 1),
        axis.title.y=element_text(size = 15),
        axis.text.y=element_text(size = 12),
        plot.margin = unit(c(0,.5,-1,0), "cm"))
p2 <- r.plot.8.fam + ggtitle("B")+ #labs(tag="B")+ 
  tanagr_select +
  theme(#legend.position = "none",
        axis.text.x=element_text(size = 10, angle = 45, hjust = 1),
        axis.title.y=element_text(size = 15),
        axis.text.y=element_text(size = 12),
        plot.margin = unit(c(-1,.5,0,0), "cm"))
#leg<- get_legend(r.plot.8.fam + tanagr_select + guides(fill=guide_legend(ncol=1)) )

#windows();ggarrange(ggarrange(p1, p2, ncol=1), leg, ncol=2)

#windows(width=6.5,height=7)

svg("Figure_5_March2025.svg", 6.5, 7)
#grid.arrange(p1 ,p2, leg, heights=c(1.1,1), widths=c(6,2), layout_matrix=rbind(c(1,3),c(2,3)),
#             padding = unit(0, "cm"))

(p1 + ylab("Relative abundance\n") + guides(fill=guide_legend(ncol=1)))+
  (p2 + ylab("Relative abundance\n") + guides(fill=guide_legend(ncol=1)))+
  plot_layout(ncol = 1, guides="collect", axis_titles="collect")

dev.off()

#########################################################
#                                                       #
# FIGURE 7                                              #
# plot cankers by state, block, inoculation, and genus  #
#                                                       #
#########################################################

cank.genus.2 <- canks.joined %>% reshape2::melt(id.vars=1:7, variable.name="SampleID", value.name="Number") %>% left_join(design, by="SampleID") %>% group_by(StateBlockInoc, Genus) %>% dplyr::summarize (count = sum(Number)) %>% arrange(StateBlockInoc, desc(count))

c.plot.state.block.inoc.genus.data <- cank.genus.2 %>%
	left_join(., cank.genus.2 %>%
			group_by (StateBlockInoc) %>%
			dplyr::summarize (t = sum(count))) %>%
			dplyr::mutate(relabund= count/t)

print(c.plot.state.block.inoc.genus.data, n = 104)

c.plot.state.block.inoc.genus.data$StateBlockInoc <- factor(c.plot.state.block.inoc.genus.data$StateBlockInoc, levels=unique(c.plot.state.block.inoc.genus.data$StateBlockInoc)[c(1,2,4,3,6,5,8,7)])

levels(c.plot.state.block.inoc.genus.data$StateBlockInoc) <-
#                        c("IN Plot 1 (Gm)", 		"IN Plot 2 (Gm)",
#													"WA Plot 1 (Control)", 	"WA Plot 1 (Gm)",
#													"WA Plot 2 (Control)", 	"WA Plot 2 (Gm)",
#													"WA Plot 3 (Control)", 	"WA Plot 3 (Gm)")
  
  c("IN P1 (Gm)", 		"IN P2 (Gm)",
    "WA P1 (Sham)", 	"WA P1 (Gm)",
    "WA P2 (Sham)", 	"WA P2 (Gm)",
    "WA P3 (Sham)", 	"WA P3 (Gm)")

n.genera <- length(unique(c.plot.state.block.inoc.genus.data$Genus))

c.plot.state.genus <- c.plot.state.block.inoc.genus.data %>%
	ggplot(aes(x=StateBlockInoc, y=relabund, fill=Genus)) +
		geom_bar(aes(), stat="identity", position="fill", size=0)+
  		xlab("") +
  		ylab("") + stackplottheme

windows(width=7, height=5);

svg("Figure_4_March2025.svg",7,5)
c.plot.state.genus+
  ylab("Relative abundance\n")+
  scale_fill_manual(values= colorss[1:n.genera+5])+
  theme(legend.text=element_text(size=9, face="italic"),
        axis.title.y=element_text(size = 15),
        axis.text.y=element_text(size = 12))
dev.off()
##############################
##                           #
## ADONIS , NMDS & HEATMAPS  #
##                           #
##############################



############################
### ADONIS ROOTS ###########
############################

# set up variables for ordinations and adonis




	r.empty.samples.logical <- rowSums(root_otu_table) > 1
	r.singletons.logical    <- colSums(root_otu_table) > 1
	r_otu.rarified <- root_otu_table[r.empty.samples.logical , r.singletons.logical]
rowSums(r_otu.rarified)>0

	r.gm    <- design[row.names(root_otu_table), "Inoculation"] %>% as.factor %>% .[r.empty.samples.logical]
	r.state <- design[row.names(root_otu_table), "State"] %>% as.factor %>% .[r.empty.samples.logical]
	r.block <- design[row.names(root_otu_table),] %>% with(paste(State,Block,sep="-")) %>% as.factor %>% .[r.empty.samples.logical]

	names_r_otu <- data.frame(otu=names(r_otu.rarified)) %>% left_join(roots.joined[,1:7])
	
	r_otu.rarified_named <- r_otu.rarified
	#names(r_otu.rarified_named) <- names_r_otu$Genus
	
	indicstate <- indicspecies::multipatt(r_otu.rarified, r.state, control = how(nperm=9999)) 
	summary(indicstate)
	
	indicstate$sign$otu <- rownames(indicstate$sign)
	indicstate$sign <- left_join(indicstate$sign, roots.joined[,1:7]) %>% filter(p.value<.05)
	arrange(indicstate$sign, p.value)%>%
	  write.csv(row.names=F, 'Revision_Jan_2025/indicspecies_tables/Table_S3_Indicspp_Roots_States.csv')
	
	indicblock <- indicspecies::multipatt(r_otu.rarified , r.block, control = how(nperm=9999)) 
	summary(indicblock)
	indicblock$sign$otu <- rownames(indicblock$sign)
	indicblock$sign <- left_join(indicblock$sign, roots.joined[,1:7]) %>% filter(p.value<.05)
	arrange(indicblock$sign, p.value)%>%
	  write.csv(row.names=F, 'Revision_Jan_2025/indicspecies_tables/Table_S4_Indicspp_Roots_Plots.csv')
	
	r_otu.gm.only <- r_otu.rarified[r.gm == "Gm",]
	r.state.gm.only <- r.state[r.gm == "Gm"]
	r.block.gm.only <- r.block[r.gm == "Gm"]
	
	
	#added to analyze necrosis in adonis
	r.necrosis.gm.only<-
	  design %>% left_join(necrosis) %>%
	  filter(SampleID %in% rownames(r_otu.gm.only)) %>%
	  filter(Inoculation == "Gm") %>%
	  summaryBy(Area ~ State + Block + Rep,., id='SampleID')# %>%
	  #na.omit
	rownames(r.necrosis.gm.only) <- r.necrosis.gm.only$SampleID 
	r.necrosis.gm.only <- r.necrosis.gm.only[rownames(r_otu.gm.only),]
	
	unique(r.necrosis.gm.only$SampleID) %>% setdiff(rownames(r_otu.rarified_named))
  rownames(r_otu.rarified_named)	%>% setdiff(unique(r.necrosis.gm.only$SampleID) )
	
  cankersizeclass <- as.factor(log(na.omit(r.necrosis.gm.only$Area.mean)) > mean(log(na.omit(r.necrosis.gm.only$Area.mean))))
  levels(cankersizeclass) <- c('small','big')
  
  indiccanks<- indicspecies::multipatt(r_otu.rarified[na.omit(r.necrosis.gm.only)$SampleID,], cankersizeclass, control = how(nperm=9999)) 
  summary(indiccanks)
  
  indiccanks$sign$otu <- rownames(indiccanks$sign)
  indiccanks$sign <- left_join(indiccanks$sign, roots.joined[,1:7]) %>% filter(p.value<=.05)
  arrange(indiccanks$sign, p.value)%>%
    write.csv(row.names=F, 'Revision_Jan_2025/indicspecies_tables/Table_S6_Indicspp_Roots_necrosis_all.csv')
  
  cankersizeclass_wa <- as.factor(log(na.omit((r.necrosis.gm.only%>%filter(State=='WA'))$Area.mean)) >= mean(log(na.omit((r.necrosis.gm.only%>%filter(State=='WA'))$Area.mean))))
  levels(cankersizeclass_wa) <- c('small','big')
  
  indiccanks_wa<- indicspecies::multipatt(
    r_otu.rarified[na.omit(r.necrosis.gm.only%>%filter(State=='WA'))$SampleID,],
    cankersizeclass_wa, control = how(nperm=99999)) 
  summary(indiccanks_wa)
  indiccanks_wa$sign$otu <- rownames(indiccanks_wa$sign)
  indiccanks_wa$sign <- left_join(indiccanks_wa$sign, roots.joined[,1:7]) #%>% filter(p.value<=.05)
  arrange(indiccanks_wa$sign, p.value) %>% head(10)
  
  cankersizeclass_in <- as.factor(log(na.omit((r.necrosis.gm.only%>%filter(State=='IN'))$Area.mean)) > mean(log(na.omit((r.necrosis.gm.only%>%filter(State=='IN'))$Area.mean))))
  levels(cankersizeclass_in) <- c('small','big')
  
  indiccanks_in<- indicspecies::multipatt(
    r_otu.rarified[na.omit(r.necrosis.gm.only%>%filter(State=='IN'))$SampleID,],
    cankersizeclass_in, control = how(nperm=99999)) 
  summary(indiccanks_in)
  indiccanks_in$sign$otu <- rownames(indiccanks_in$sign)
  indiccanks_in$sign <- left_join(indiccanks_in$sign, roots.joined[,1:7]) %>% filter(p.value<=.05)
  arrange(indiccanks_in$sign, p.value)%>%
    write.csv(row.names=F, 'Revision_Jan_2025/indicspecies_tables/Table_S7_Indicspp_Roots_necrosis_IN.csv')
    	################################ added Jan 2025
    	# overall adonis cankers roots
  
	leavout.gm <- which(is.na(r.necrosis.gm.only$Area.mean))
	
	#indicspecies::multipatt(r_otu.rarified_named[-leavout.gm,], r.necrosis.gm.only$Area.mean[-leavout.gm])
	
    	d.gm <- vegdist(
    	  #decostand(r_otu.gm.only[-leavout.gm,-leavout.gm],method="total"),
    	  r_otu.gm.only[-leavout.gm,-leavout.gm],
    	  method='jaccard')
    	  #method='robust.aitchison')
    	
    	
    	# adonis w necrosis
    	adonis2(d.gm ~ r.necrosis.gm.only$Area.mean[-leavout.gm], by='term', permutations=9999)
    	adonis2(d.gm ~ r.necrosis.gm.only$Area.mean[-leavout.gm]+
    	          r.state.gm.only[-leavout.gm], by='term', permutations=9999)
    	adonis2(d.gm ~ r.necrosis.gm.only$Area.mean[-leavout.gm]+
    	          r.state.gm.only[-leavout.gm]+
    	          r.state.gm.only[-leavout.gm]:
    	          r.block.gm.only[-leavout.gm], by='term', permutations=9999)

    	# without necrosis
    	adonis2(d.gm ~ #r.necrosis.gm.only$Area.mean[-leavout.gm]+
    	          r.state.gm.only[-leavout.gm]+
    	          r.state.gm.only[-leavout.gm]:
    	          r.block.gm.only[-leavout.gm], by='term', permutations=9999)
    	
	
	r.block.gm.only <- r.block[r.gm == "Gm"]

	r_otu.wa.only <- r_otu.rarified[r.state == "WA",]
	r.gm.wa.only <- r.gm[r.state == "WA"]
	r.block.wa.only <- r.block[r.state == "WA"]
	
	r_otu.wa.gm.only <- r_otu.rarified[r.state == "WA" & r.gm == "Gm",]

# distance matrix jaccard on gm vs control in WA

	d.gm.wa <- vegdist(decostand(r_otu.wa.only, method="total"), method='jaccard')
	
	#adonis
	adonis2(d.gm.wa ~ r.block.wa.only + r.gm.wa.only, permutations=9999, by='terms')
	#adonis2(d.gm.wa ~ r.block.wa.only : r.gm.wa.only+ r.gm.wa.only, permutations=9999)
	#adonis2(d.gm.wa ~ r.block.wa.only : r.gm.wa.only+ r.block.wa.only, permutations=9999)
		
		
		#adonis2(d.gm.wa ~ r.block.wa.only : r.gm.wa.only+ r.gm.wa.only + r.block.wa.only , permutations=9999)

	#adonis2(d.gm.wa ~ r.block.wa.only * r.gm.wa.only, permutations=9999)
	
	## results
	adonis2(d.gm.wa ~ r.block.wa.only + r.gm.wa.only, permutations=9999, by='terms')
	#adonis2(d.gm.wa ~ r.block.wa.only, permutations=9999)
	#adonis2(d.gm.wa ~ r.gm.wa.only, permutations=9999)

# distance matrix jaccard on WA vs IN (gm-only)

	d.in.wa <- vegdist(decostand(r_otu.gm.only, method="total"), method='jaccard')
	d.in.wa.n <- vegdist(decostand(r_otu.gm.only[-leavout.gm,], method="total"), method='jaccard')
	
	#adonis
	
	#results
	adonis2(d.in.wa ~ r.state.gm.only + r.block.gm.only, permutations=9999)
	adonis2(d.in.wa ~ r.state.gm.only * r.block.gm.only, permutations=9999)
	
	# WITH NECROSIS in WA
	str(d.in.wa.n)
	length(r.state.gm.only[-leavout.gm])
	length(r.block.gm.only[-leavout.gm])
	length(r.necrosis.gm.only$Area.mean[-leavout.gm])
	adonis2(d.in.wa.n ~ r.state.gm.only[-leavout.gm] + r.block.gm.only[-leavout.gm], permutations=9999, by='terms')
	adonis2(d.in.wa.n ~
	          r.state.gm.only[-leavout.gm] +
	          r.block.gm.only[-leavout.gm] + 
	          r.necrosis.gm.only$Area.mean[-leavout.gm], permutations=9999, by='terms')

	#adonis2(d.in.wa ~ r.state.gm.only, permutations=9999)

	
	sin<-r.state.gm.only[-leavout.gm] == "IN"
	swa<-r.state.gm.only[-leavout.gm] == "WA"
	d.in.wa.nin <- vegdist(decostand(r_otu.gm.only[-leavout.gm,][sin,], method="total"), method='jaccard')
	d.in.wa.nwa <- vegdist(decostand(r_otu.gm.only[-leavout.gm,][swa,], method="total"), method='jaccard')

	# IN
	adonis2(d.in.wa.nin ~
	          r.block.gm.only[-leavout.gm][sin] + 
	          r.necrosis.gm.only$Area.mean[-leavout.gm][sin], permutations=9999, by='terms')
  # WA
	adonis2(d.in.wa.nwa ~
	          r.block.gm.only[-leavout.gm][swa] + 
	          r.necrosis.gm.only$Area.mean[-leavout.gm][swa], permutations=9999, by='terms')
	
	state.gmi<-r.state.gm.only[-leavout.gm]
	  block.gmi<-r.block.gm.only[-leavout.gm]
	  area.gmi<-r.necrosis.gm.only$Area.mean[-leavout.gm]
	  
	# INTERACTIVE
	adonis2(d.in.wa.n ~
	          state.gmi *
	          block.gmi * 
	          area.gmi, permutations=9999)
	adonis2(d.in.wa.n ~
	          state.gmi +
	          block.gmi + 
	          area.gmi +
	          state.gmi*area.gmi, permutations=9999)
# distance matrix jaccard on everything

	d.all <- vegdist(decostand(r_otu.rarified, method="total"), method='jaccard')

	#ordination
	root.nmds <- metaMDS(decostand(r_otu.rarified, method="total"), distance="jaccard", k =3, trymax=250)
	root.nmds$points


############################
### ADOINIS CANKS ##########
############################

# set up variables for ordinations and adonis

	c.empty.samples.logical <- rowSums(canker_otu_table) > 1
	c.singletons.logical    <- colSums(canker_otu_table) > 1
	c_otu.rarified <- canker_otu_table[c.empty.samples.logical , c.singletons.logical]

	c.gm    <- design[row.names(canker_otu_table), "Inoculation"] %>% as.factor %>% .[c.empty.samples.logical]
	c.state <- design[row.names(canker_otu_table), "State"] %>% as.factor %>% .[c.empty.samples.logical]
	c.block <- design[row.names(canker_otu_table),] %>% with(paste(State,Block,sep="-")) %>% as.factor %>% .[c.empty.samples.logical]

	c_otu.gm.only <- c_otu.rarified[c.gm == "Gm",]
	c.state.gm.only <- c.state[c.gm == "Gm"]
	c.block.gm.only <- c.block[c.gm == "Gm"]

	c.necrosis.gm.only<-
	  design %>% left_join(necrosis) %>%
	  filter(SampleID %in% rownames(c_otu.gm.only)) %>%
	  filter(Inoculation == "Gm") %>%
	  summaryBy(Area ~ State + Block + Rep,., id='SampleID')# %>%
	#na.omit
	rownames(c.necrosis.gm.only) <- c.necrosis.gm.only$SampleID 
	c.necrosis.gm.only <- c.necrosis.gm.only[rownames(c_otu.gm.only),]
	leavout.gm.canker <- which(is.na(c.necrosis.gm.only$Area.mean))
	
	c.gm.comm <- vegdist(decostand(c_otu.gm.only, method="total"), method='jaccard')
	
	adonis2(c.gm.comm ~ c.necrosis.gm.only$Area.mean, permutations = 9999, by = 'terms')
	
	adonis2(c.gm.comm ~ c.necrosis.gm.only$Area.mean, permutations = 9999, by = 'terms')
	adonis2(c.gm.comm ~ c.necrosis.gm.only$Area.mean+
	          c.necrosis.gm.only$State, permutations = 9999, by = 'terms')
	adonis2(c.gm.comm ~
	          c.necrosis.gm.only$Area.mean+
	          c.necrosis.gm.only$State+
	          c.necrosis.gm.only$State:c.necrosis.gm.only$Block,
	        permutations = 9999, by = 'terms')
	adonis2(c.gm.comm ~
	    #      c.necrosis.gm.only$Area.mean+
	          c.necrosis.gm.only$State+
	          c.necrosis.gm.only$State:c.necrosis.gm.only$Block,
	        permutations = 9999, by = 'terms')
	
	c_otu.wa.only <- c_otu.rarified[c.state == "WA",]
	c.gm.wa.only <- c.gm[c.state == "WA"]
	c.block.wa.only <- c.block[c.state == "WA"]

	dim(c_otu.wa.only)
	length(c.gm.wa.only)

# distance matrix jaccard on gm vs control in WA

	cd.gm.wa <- vegdist(decostand(c_otu.wa.only, method="total"), method='jaccard')

	#adonis
	
	# results
	dim(cd.gm.wa)
	length(c.block.wa.only)
	dim(c.gm.wa.only)
	#adonis2(cd.gm.wa ~ c.block.wa.only * c.gm.wa.only, permutations=9999)
	adonis2(cd.gm.wa ~ c.block.wa.only + c.gm.wa.only, permutations=9999, by='terms')
	#adonis2(cd.gm.wa ~ c.gm.wa.only, permutations=9999)

# distance matrix jaccard on WA vs IN (gm-only)

	cd.in.wa <- vegdist(decostand(c_otu.gm.only, method="total"), method='jaccard')

	#adonis
	adonis2(cd.in.wa ~ c.state.gm.only + c.block.gm.only, permutations=9999, by='terms')
	#adonis2(cd.in.wa ~ c.state.gm.only, permutations=9999)
	str(cd.in.wa)
	length(c.necrosis.gm.only$Area.mean)

	adonis2(cd.in.wa ~
	          c.state.gm.only +
	          c.block.gm.only + 
	          c.necrosis.gm.only$Area.mean , permutations=9999)
	
	####### canker indicators
	
	
	c.indicstate <- indicspecies::multipatt(c_otu.rarified, c.state, control = how(nperm=9999)) 
	summary(c.indicstate)
	
	c.indicstate$sign$otu <- rownames(c.indicstate$sign)
	c.indicstate$sign <- left_join(c.indicstate$sign, canks.joined[,1:7]) %>% filter(p.value<.05)
	arrange(c.indicstate$sign, p.value) %>% write.csv(row.names=F, 'Revision_Jan_2025/indicspecies_tables/Table_S1_Indicspp_Cankers_States.csv')
	
	c.indicblock <- indicspecies::multipatt(c_otu.rarified , c.block, control = how(nperm=9999)) 
	summary(c.indicblock)
	c.indicblock$sign$otu <- rownames(c.indicblock$sign)
	c.indicblock$sign <- left_join(c.indicblock$sign, canks.joined[,1:7]) %>% filter(p.value<.05)
	arrange(c.indicblock$sign, p.value)%>%
	  write.csv(row.names=F, 'Revision_Jan_2025/indicspecies_tables/Table_S2_Indicspp_Cankers_Plots.csv')
	
	c.cankersizeclass <- as.factor(log(na.omit(c.necrosis.gm.only$Area.mean)) > mean(log(na.omit(c.necrosis.gm.only$Area.mean))))
	levels(c.cankersizeclass) <- c('small','big')
	
	c.indiccanker_all <- indicspecies::multipatt(c_otu.gm.only , c.cankersizeclass, control = how(nperm=9999)) 
	summary(c.indiccanker_all)
	c.indiccanker_all$sign$otu <- rownames(c.indiccanker_all$sign)
	c.indiccanker_all$sign <- left_join(c.indiccanker_all$sign, canks.joined[,1:7]) #%>% filter(p.value<.1)
	arrange(c.indiccanker_all$sign, p.value)%>%
	  write.csv(row.names=F, 'Revision_Jan_2025/indicspecies_tables/Table_S5_Indicspp_Cankers_necrosis_all.csv')
	
	c.cankersizeclass_wa <- as.factor(log(na.omit((c.necrosis.gm.only%>%filter(State=='WA'))$Area.mean)) >= mean(log(na.omit((c.necrosis.gm.only%>%filter(State=='WA'))$Area.mean))))
	levels(c.cankersizeclass_wa) <- c('small','big')
	
	c.indiccanks_wa<- indicspecies::multipatt(
	  c_otu.gm.only[na.omit(c.necrosis.gm.only%>%filter(State=='WA'))$SampleID,],
	  c.cankersizeclass_wa, control = how(nperm=99999)) 
	summary(c.indiccanks_wa)
	
	c.cankersizeclass_in <- as.factor(log(na.omit((c.necrosis.gm.only%>%filter(State=='IN'))$Area.mean)) >= mean(log(na.omit((c.necrosis.gm.only%>%filter(State=='IN'))$Area.mean))))
	levels(c.cankersizeclass_in) <- c('small','big')
	
	c.indiccanks_in<- indicspecies::multipatt(
	  c_otu.gm.only[na.omit(c.necrosis.gm.only%>%filter(State=='IN'))$SampleID,],
	  c.cankersizeclass_in, control = how(nperm=99999)) 
	summary(c.indiccanks_in)
	
#######################
##### NMDS. ###########
#######################

# distance matrix jaccard on everything

	cd.all <- vegdist(decostand(c_otu.rarified, method="total"), method='jaccard')

	#ordination
	cank.nmds <- metaMDS(decostand(c_otu.rarified, method="total"), distance="jaccard", k =3, trymax=1500)
	cank.nmds$points

	
#############################	
###### ECOSPHERE VERSION ####
#############################
	
## FIG 3	
	
	# redo nmds
	root.nmds2 <- metaMDS(decostand(r_otu.rarified, method="total"), distance="jaccard", k =2, trymax=500)
	root.nmds2 <- metaMDS(decostand(r_otu.rarified, method="total"), distance="jaccard", k =2, trymax=500)
	cank.nmds2 <- metaMDS(decostand(c_otu.rarified, method="total"), distance="jaccard", k =2, trymax=500)
	cank.nmds2 <- metaMDS(decostand(c_otu.rarified, method="total"), distance="jaccard", k =2, trymax=500)
	
	# pcoa
	root.pcoa <- pcoa(vegdist(decostand(r_otu.rarified, method="total"), "jaccard"))
  cank.pcoa <- pcoa(vegdist(decostand(c_otu.rarified, method="total"), "jaccard"))
	
  # cca
  
  design %>% head
  
  cankers<-read.csv("Phenotypes/WA-IN-Transplant2019-20.Necrosis.csv" )  %>%
    summaryBy(Area ~ State + Block + Plant, .) %>%
    rename(Rep=Plant) %>%
    left_join(read.csv("Phenotypes/MorphoGH2020Expt1.csv")[,1:6]) %>%
    left_join(design)
  
  rownames(cankers)<-cankers$SampleID
  
### ROOTS  
  rowSums(r_otu.rarified)
  rows.shared <-
    rownames(filter(na.omit(cankers), Inoculation=="Gm")) %>%
    intersect(rownames(r_otu.rarified))# %>% filter(rowSums(.) > 0)))
  
  #selection <- rownames(r_otu.rarified) %in% rows.shared

  community.data.roots.gm <-
    decostand(r_otu.rarified[rows.shared,], method="hellinger", MARGIN=1)
#    decostand(r_otu.rarified[rows.shared,], method="standardize", MARGIN=1)
  rowSums(community.data.roots.gm)
  Xr <- #r_otu.rarified[rows.shared,]
    community.data.roots.gm
  Yr <- cankers[rows.shared,] %>% dplyr::select(Area.mean, Caliper, Displacement) %>%
    decostand("log", MARGIN=2) %>% decostand("standardize")
  Zr <- cankers[rows.shared,] %>% dplyr::select(State,Block)
  dim(Xr);dim(Yr);dim(Zr)

  roots.cca <- cca(Xr,Yr,Zr)
  summary(roots.cca)
  RsquareAdj(roots.cca)
  roots.rda <- rda(Xr,Yr,Zr,scale=T)
  roots.pcoa.gm <- pcoa(vegdist(Xr, method="jaccard"))
  
  #windows()
  p.roots.sites<-plot(roots.cca, scaling=3)$sites %>%
    as.data.frame %>%
    dplyr::rename(Axis.1 = CCA1)%>% dplyr::rename(Axis.2=CCA2)
  p.roots.sites<-p.roots.sites%>%
    mutate(SampleID=rownames(p.roots.sites)) %>%
    separate(SampleID, into = c("State","Block","Rep"))
  
  
  p.roots.sp<-plot(roots.cca, scaling=2)$species %>%
    as.data.frame %>%
    dplyr::rename(Axis.1 = CCA1)%>% dplyr::rename(Axis.2=CCA2)
  p.roots.sp<-p.roots.sp%>%
    mutate(otu=rownames(p.roots.sp)) %>%
    left_join(roots.joined, by="otu") %>%
    separate(Taxon.assignment, into=c('genus','epithet'), sep="[\\_\\,]", remove=F)%>%
    mutate(sp_abbr = paste(substr(genus,1,1), ".", epithet, sep=""))
  p.roots.sp<-p.roots.sp%>%
    mutate(hj =
             1-(Axis.1/Axis.2-min(Axis.1/Axis.2))/(max(Axis.1/Axis.2)-min(Axis.1/Axis.2))
#             -(sign(Axis.1)-1)/2
           )%>%
    mutate(vj =
             1-(Axis.2/Axis.1-min(Axis.2/Axis.1))/(max(Axis.2/Axis.1)-min(Axis.2/Axis.1))
           )
  
  p.roots.cca <- plot(roots.cca, scaling=0)$biplot %>%
    as.data.frame %>%
    dplyr::rename(Axis.1 = CCA1)%>% dplyr::rename(Axis.2=CCA2)
  p.roots.cca <- p.roots.cca%>%
    mutate(name=rownames(p.roots.cca))
  p.roots.cca <- p.roots.cca%>%
    mutate(hj =
             1-(Axis.1/Axis.2-min(Axis.1/Axis.2))/(max(Axis.1/Axis.2)-min(Axis.1/Axis.2))
    )%>%
    mutate(vj =
             1-(Axis.2/Axis.1-min(Axis.2/Axis.1))/(max(Axis.2/Axis.1)-min(Axis.2/Axis.1))
    )
  
  ggplot() +
    geom_point(data=p.roots.sites, aes(x=Axis.1, y=Axis.2, colour=State), size=3) +
    ggplot2::stat_ellipse(data=p.roots.sites, aes(x=Axis.1, y=Axis.2, colour=State), geom="polygon", fill=NA, alpha=0.2, type="t", level=.48)
    
  # try a facet
  ggplot() +
    geom_point(data=p.roots.sites, aes(x=Axis.1, y=Axis.2, colour=Block), size=3) +
    ggplot2::stat_ellipse(data=p.roots.sites, aes(x=Axis.1, y=Axis.2, colour=Block), geom="polygon", fill=NA, alpha=0.2, type="t", level=.48)+
    facet_wrap('State')
  
  
  ggca1 <- ggplot() +
    geom_segment(data=p.roots.sp, aes(xend=Axis.1, yend=Axis.2),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0, col="black")+
    geom_text_repel(
      data=p.roots.sp %>%
        filter(abs(Axis.1)>=.5 | abs(Axis.2)>=.5),
      mapping=
        aes(x=Axis.1,
            y=Axis.2,
            label=Genus,
       #     hjust=hj,
      #      vjust=vj
            ),
      vjust="outward",
      hjust="outward",
      #nudge_x=.1,
      #nudge_y=.1,
      cex=3,
      direction="y",
      box.padding=0,
      min.segment.length=0.1
      )+
    geom_segment(data=p.roots.cca, aes(xend=Axis.1, yend=Axis.2),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0, col="blue")+
    geom_text(
      data=p.roots.cca,
      mapping=
        aes(x=Axis.1,
            y=Axis.2,
            label=name#,
            #hjust=hj,
            #vjust=vj
        ),
      vjust="outward",
      hjust="outward",
     # nudge_x=.5,
    #  nudge_y=.5,
      cex=5, col="blue")+
    xlim(-1.75,1.75)+ylim(-1.25,2.5)+
    labs(x = "\nCCA1", y = "CCA2\n")+
    theme_bw()
  
  ggca1
  
  p.roots.biplot.all <-
    rbind(
      cbind(p.roots.sp%>%dplyr::rename(name=Genus)%>%dplyr::select(Axis.1,Axis.2,name,hj,vj), level='sp', lcol="black", sz=3),
      cbind(p.roots.cca,level='cov',lcol="blue", sz=5)
    )
  
  ggcca <- ggplot() +
    geom_segment(data=p.roots.biplot.all, aes(xend=Axis.1, yend=Axis.2, col=lcol),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0)+
    geom_text(
      data=p.roots.biplot.all %>%
        filter((abs(Axis.1)>=.5 | abs(Axis.2)>=.5) | level=="cov"),
      mapping=
        aes(x=Axis.1*1.1,
            y=Axis.2*1.1,
            label=name,
            hjust=hj,
            vjust=vj,
            col=lcol
        ),
      #vjust="outward",
      #hjust="outward",
      #direction="y",
      #box.padding=0.1
    )+
    xlim(-1.75,1.75)+ylim(-1.25,2.5)+
    labs(x = "\nCCA1", y = "CCA2\n")+
    theme_bw()+
    theme(legend.position="none")
  ggcca
  
## CANKERS

  rows.shared.c <-
    rownames(filter(na.omit(cankers), Inoculation=="Gm")) %>%
    intersect(rownames(c_otu.rarified))
  
  community.data.canks.gm <-
    decostand(c_otu.rarified[rows.shared.c,], method="total", MARGIN=1)
  
  Xc <- #r_otu.rarified[rows.shared,]
    community.data.canks.gm
  Yc <- cankers[rows.shared.c,] %>% dplyr::select(Area.mean, Caliper, Displacement) %>%
    decostand("log", MARGIN=2) %>% decostand("standardize")
  Zc <- cankers[rows.shared.c,] %>% dplyr::select(State,Block)
  dim(Xc);dim(Yc);dim(Zc)
  
  canks.cca <- cca(Xc,Yc,Zc)
  RsquareAdj(canks.cca)
  canks.rda <- rda(Xc,Yc,Zc,scale=T)
  #canks.pcoa.gm <- pcoa(vegdist(X, method="jaccard"))
  
  #windows()
  p.canks.sites<-plot(canks.cca, scaling=3)$sites %>%
    as.data.frame %>%
    dplyr::rename(Axis.1 = CCA1)%>% dplyr::rename(Axis.2=CCA2)
  p.canks.sites<-p.canks.sites%>%
    mutate(SampleID=rownames(p.canks.sites)) %>%
    separate(SampleID, into = c("State","Block","Rep"))
  
  
  p.canks.sp<-plot(canks.cca, scaling=2)$species %>%
    as.data.frame %>%
    dplyr::rename(Axis.1 = CCA1)%>% dplyr::rename(Axis.2=CCA2)
  p.canks.sp<-p.canks.sp%>%
    mutate(otu=rownames(p.canks.sp)) %>%
    left_join(canks.joined, by="otu") %>%
    separate(Taxon.assignment, into=c('genus','epithet'), sep="[\\_\\,]", remove=F)%>%
    mutate(sp_abbr = paste(substr(genus,1,1), ".", epithet, sep=""))
  p.canks.sp<-p.canks.sp%>%
    mutate(hj =
             1-(Axis.1/Axis.2-min(Axis.1/Axis.2))/(max(Axis.1/Axis.2)-min(Axis.1/Axis.2))
           #             -(sign(Axis.1)-1)/2
    )%>%
    mutate(vj =
             1-(Axis.2/Axis.1-min(Axis.2/Axis.1))/(max(Axis.2/Axis.1)-min(Axis.2/Axis.1))
    )
  
  p.canks.cca <- plot(canks.cca, scaling=0)$biplot %>%
    as.data.frame %>%
    dplyr::rename(Axis.1 = CCA1)%>% dplyr::rename(Axis.2=CCA2)
  p.canks.cca <- p.canks.cca%>%
    mutate(name=rownames(p.canks.cca))
  p.canks.cca <- p.canks.cca%>%
    mutate(hj =
             1-(Axis.1-min(Axis.1))/(max(Axis.1)-min(Axis.1))
    )%>%
    mutate(vj =
             1-(Axis.2-min(Axis.2))/(max(Axis.2)-min(Axis.2))
    )

  ggca2 <- ggplot() +
    geom_segment(data=p.canks.sp, aes(xend=Axis.1, yend=Axis.2),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0, col="black")+
    geom_text_repel(
      data=p.canks.sp ,
      mapping=
        aes(x=Axis.1,
            y=Axis.2,
            label=sp_abbr,
            #     hjust=hj,
            #      vjust=vj
        ),
      vjust="outward",
      hjust="outward",
      #nudge_x=.1,
      #nudge_y=.1,
      cex=3,
      direction="y",
      box.padding=0,
      min.segment.length=0.1
    )+
    geom_segment(data=p.canks.cca, aes(xend=Axis.1, yend=Axis.2),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0, col="blue")+
    geom_text(
      data=p.canks.cca,
      mapping=
        aes(x=Axis.1*1.25,
            y=Axis.2*1.25,
            label=name#,
           # hjust=hj,
          #  vjust=vj
        ),
      vjust="outward",
      hjust="outward",
   #   nudge_x=.5,
   #   nudge_y=.5,
      cex=5, col="blue")+
    xlim(-2.25,1.5)+ylim(-1.25,1)+
    labs(x = "\nCCA1", y = "CCA2\n")+
    theme_bw()
  
  p.roots.sp %>% dplyr::select(Axis.1,Axis.2,Genus)
  
#### put it together
  allrows<-
    bind_rows(
      cbind(p.roots.sp%>%dplyr::select(Axis.1,Axis.2,Genus,hj,vj), panel="Roots", part="Spp Scores"),
      cbind(p.roots.cca, panel="Roots", part="Pheno Scores"),
      cbind(p.canks.sp%>%dplyr::select(Axis.1,Axis.2,sp_abbr,hj,vj), panel="Cankers", part="Spp Scores"),
      cbind(p.canks.cca, panel="Cankers", part="Pheno Scores"))

  scaleFUN <- function(x) sprintf("%.2f", x)
  scaleFUN(.5)
  
  ggca_all <-
    ggplot(data=allrows %>% filter(Axis.2<2)) +
    geom_segment(data= ~ filter(.x, panel=="Roots", part=="Spp Scores"), aes(xend=Axis.1, yend=Axis.2, col=part),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0)+#, col="black")+
    geom_text_repel(
      data= ~ filter(.x, panel=="Roots", part=="Spp Scores", abs(Axis.1)>=.5 | abs(Axis.2)>=.5),
      mapping=aes(x=Axis.1,y=Axis.2,label=Genus, col=part),
      vjust="outward",hjust="outward",cex=3,direction="y",box.padding=0,min.segment.length=0.1)+
    geom_segment(data= ~ filter(.x, panel=="Roots", part=="Pheno Scores"), aes(xend=Axis.1, yend=Axis.2, col=part),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0)+#, col="blue")+
    geom_text(
      data= ~ filter(.x, panel=="Roots", part=="Pheno Scores"),
      mapping=aes(x=Axis.1,y=Axis.2,label=name, col=part, size=2),
      vjust="outward",hjust="outward",cex=5)+#, col="blue") +
    geom_segment(data=~ filter(.x, panel=="Cankers", part=="Spp Scores"), aes(xend=Axis.1, yend=Axis.2, col=part),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0)+#, col="black")+
    geom_text_repel(
      data=~ filter(.x, panel=="Cankers", part=="Spp Scores"),
      mapping=aes(x=Axis.1,y=Axis.2,label=sp_abbr, col=part),
      vjust="outward", hjust="outward",cex=3,direction="y",box.padding=0)+#,min.segment.length=1000)+
    geom_segment(data=~ filter(.x, panel=="Cankers", part=="Pheno Scores"), aes(xend=Axis.1, yend=Axis.2, col=part),
                 arrow=arrow(length=unit(.05, "inches")), x=0, y=0)+#, col="blue")+
    geom_text(
      data=~ filter(.x, panel=="Cankers", part=="Pheno Scores"),
      mapping=aes(x=Axis.1*1.25,y=Axis.2*1.25,label=name, col=part, size=2),
      vjust="outward",hjust="outward",cex=5)+#, col="blue")+
    xlim(-2,2)+
    guides(colour=guide_legend(title="Eigenvalue"))+labs(x="\nCCA1",y="CCA2\n")+
    facet_wrap(~panel, scales="free_y")+
    theme_bw()#+
    #theme(strip.background = element_blank())

  windows(10,5)
  ggca_all
## NMDS
  
	# roots
	# r.data.scores <- vegan::scores(root.nmds2)  #Using the scores function from vegan to extract the site scores and convert to a data.frame

	# optionally, pcoa
	r.data.scores <- root.pcoa
	r.data.scores$sites <- r.data.scores$vectors[,1:2]

	names(r.data.scores$sites)[1:2] <- c("Axis.1","Axis.2")
	r.data.scores$sites <- as.data.frame(r.data.scores$sites)
	r.data.scores$sites$site <- rownames(r.data.scores$sites) # create a column of site names, from the rownames of data.score
	r.data.scores$species <- as.data.frame(r.data.scores$species)
	r.data.scores$species$sp <- rownames(r.data.scores$species) # create a column of site names, from the rownames of data.scores
	r.data.scores$sites <- r.data.scores$sites %>% separate(site, into = c("State","Block","Rep"), sep="-")
	r.data.scores$sites$Rep <- as.integer(r.data.scores$sites$Rep)
	str(r.data.scores$sites)
	str(design)
	r.data.scores$sites <- r.data.scores$sites %>% left_join(design)
	r.data.scores$sites <- r.data.scores$sites %>% mutate(grp = paste(State, Block, Inoculation))
	#data.scores$grp <- bioregions_vector$biogeographic_region  %>%
	#  gsub("\\.", " ", .) %>% gsub("([A-Z]{1}[a-z]+)([A-Z]{1}[a-z]+)","\\1 \\2", .)#  add the grp variable created earlier
	#head(data.scores)  #look at the data
	
	# cankers
	# c.data.scores <- vegan::scores(cank.nmds2)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
	
	# optionally, pcoa
	c.data.scores <- cank.pcoa
	c.data.scores$sites <- c.data.scores$vectors[,1:2]
	
	names(c.data.scores$sites)[1:2] <- c("Axis.1","Axis.2")
	c.data.scores$sites <- as.data.frame(c.data.scores$sites)
	c.data.scores$sites$site <- rownames(c.data.scores$sites) # create a column of site names, from the rownames of data.score
	c.data.scores$species <- as.data.frame(c.data.scores$species)
	c.data.scores$species$sp <- rownames(c.data.scores$species) # create a column of site names, from the rownames of data.scores
	c.data.scores$sites <- c.data.scores$sites %>% separate(site, into = c("State","Block","Rep"), sep="-")
	c.data.scores$sites$Rep <- as.integer(c.data.scores$sites$Rep)
	str(c.data.scores$sites)
	str(design)
	c.data.scores$sites <- c.data.scores$sites %>% left_join(design)
	c.data.scores$sites <- c.data.scores$sites %>% mutate(grp = paste(State, Block, Inoculation))
	
	# for nmds
	all.datscores <- cbind(r.data.scores$sites[-27,], tissue="Roots") %>% rbind(cbind(c.data.scores$sites[-1,], tissue="Cankers"))
	
	# optionall, for pcoa
	all.datscores <- cbind(r.data.scores$sites, tissue="Roots") %>% rbind(cbind(c.data.scores$sites, tissue="Cankers"))
	
	levels(factor(all.datscores$StateBlock))
	levels(factor(all.datscores$StateBlock))
	
	all.datscores$StateBlock <- factor(all.datscores$StateBlock)
	levels(all.datscores$StateBlock)<-c(
	                                     "IN P1",
	                                     "IN P2",
	                                     "WA P1",
	                                     "WA P2",
	                                     "WA P3"
	                                   )
	
	#levels(as.factor(data.scores$grp))
	#r.data.scores$sites
	#gg1<- ggplot() +
	#  geom_point(data=r.data.scores$sites, aes(x=NMDS1, y=NMDS2, colour=treat_long, pch=StateBlock), size=3) +
	#  ggplot2::stat_ellipse(data=r.data.scores$sites, aes(x=NMDS1, y=NMDS2, colour=treat_long, linetype=StateBlock), geom="polygon", fill=NA, alpha=0.2, type="t", level=.48) +
	  #  ggplot2::stat_ellipse(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=grp), linetype=3, type="t", level=.95) +
	  #scale_fill_manual("Unifrac Floristic Cluster", values=c('black','green','red')) +
	  #scale_color_manual("Unifrac Floristic Cluster", values=c('black','green','red')) +
	#  coord_cartesian()+#xlim=c(-.5,.6), ylim=c(-.4,.5)) +
	#  theme_bw()
	
	gg2<- ggplot()+
	  geom_point(data=all.datscores,aes(x=Axis.1, y=Axis.2, colour=treat_long, pch=StateBlock, fill=treat_long), size=3) +
	  ggplot2::stat_ellipse(data=all.datscores, aes(x=Axis.1, y=Axis.2, colour=treat_long, linetype=StateBlock), geom="polygon", fill=NA, alpha=0.2, type="t", level=.48) +
	  #  ggplot2::stat_ellipse(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=grp), linetype=3, type="t", level=.95) +
	  #scale_fill_manual("Unifrac Floristic Cluster", values=c('black','green','red')) +
	  #scale_color_manual("Unifrac Floristic Cluster", values=c('black','green','red')) +
	  coord_cartesian()+#xlim=c(-.5,.6), ylim=c(-.4,.5)) +
	    scale_shape_manual(values=21:25)+
	  xlab("\nPCOA1\n")+ylab("PCOA2\n")+
	  guides(colour=guide_legend(title="State & Inoc"))+
	  guides(pch=guide_legend(title="State & Block"))+
	  guides(linetype=guide_legend(title="State & Block"))+
	  guides(fill=guide_legend(title="State & Inoc"))+
	  theme_bw()+
	  facet_wrap(~ tissue, scales="free")
	
	
	# try another facet Jan 2025
	
	gg3 <-
	  ggplot()+
	  geom_point(data=all.datscores,aes(x=Axis.1, y=Axis.2, colour=treat_long, pch=StateBlock, fill=treat_long), size=3) +
	  ggplot2::stat_ellipse(data=all.datscores, aes(x=Axis.1, y=Axis.2, colour=treat_long, linetype=StateBlock), geom="polygon", fill=NA, alpha=0.2, type="t", level=.48) +
	  #  ggplot2::stat_ellipse(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=grp), linetype=3, type="t", level=.95) +
	  #scale_fill_manual("Unifrac Floristic Cluster", values=c('black','green','red')) +
	  #scale_color_manual("Unifrac Floristic Cluster", values=c('black','green','red')) +
	  coord_cartesian()+#xlim=c(-.5,.6), ylim=c(-.4,.5)) +
	  scale_shape_manual(values=21:25)+
	  xlab("\nPCOA1\n")+ylab("PCOA2\n")+
	  guides(colour=guide_legend(title="State & Inoc"))+
	  guides(pch=guide_legend(title="State & Block"))+
	  guides(linetype=guide_legend(title="State & Block"))+
	  guides(fill=guide_legend(title="State & Inoc"))+
	  theme_bw()+
	  facet_grid(StateBlock~tissue);windows(8,5);gg3
	
	#mean(necrosis$Area)
	# ~ 49.5
	
	# Jan 2025
	one.more.thing.to.try <-
	  all.datscores %>% 
	  mutate(SampleID = gsub("-","",SampleID)) %>% inner_join(necr.summary) %>%
	  #mutate(Area.8.bins = as.factor(Area.8.bins)) %>%
	  #mutate(Area.2.bins = (Area.mean > 50) %>% factor(ordered=T))
	  (function (x)
	    x %>% left_join(
	      (group_by(x, tissue, StateBlock, Inoculation) %>%
	         summarise(group_median = median(Area.mean),
	                   transformed_group_mean = exp(mean(log(Area.mean))))%>%ungroup())
	    )) %>%  mutate(Area.2.bins.median = (Area.mean > group_median) %>% factor(ordered=T)) %>%
	  mutate(Area.2.bins.log.mean = (Area.mean > transformed_group_mean) %>% factor(ordered=T))
	
	levels(one.more.thing.to.try$Area.2.bins.median) <- c('lower 50%','bottom 50%')
	levels(one.more.thing.to.try$Area.2.bins.log.mean) <- c('Small (< plot logmean)','Large (\u2265 plot logmean)')
	
	r_otu.rarified_named2<-r_otu.rarified_named
	
	rownames(r_otu.rarified_named2) <- gsub('-','',rownames(r_otu.rarified_named2) )
	
	one.more.thing.to.try_in1 <- one.more.thing.to.try %>%
	  dplyr::select(StateBlock, SampleID, Area.2.bins.log.mean) %>%
	  filter(SampleID %in% rownames(r_otu.rarified_named2)) %>%
	  filter(StateBlock == "IN P1") %>% distinct
	indicnecr_in1 <- indicspecies::multipatt(r_otu.rarified_named2[one.more.thing.to.try_in1$SampleID,], one.more.thing.to.try_in1$Area.2.bins.log.mean)#$sign
	summary(indicnecr_in1)
	#indicnecr_in1$otu <- rownames(indicnecr_in1)
	#indicnecr_in1 <- left_join(indicnecr_in1, roots.joined[,1:7]) %>% filter(p.value<.05)
	#arrange(indicnecr_in1, p.value)
	
	one.more.thing.to.try_in2 <- one.more.thing.to.try %>%
	  dplyr::select(StateBlock, SampleID, Area.2.bins.log.mean) %>%
	  filter(SampleID %in% rownames(r_otu.rarified_named2)) %>%
	  filter(StateBlock == "IN P2") %>% distinct
	indicnecr_in2 <- indicspecies::multipatt(r_otu.rarified_named2[one.more.thing.to.try_in2$SampleID,], one.more.thing.to.try_in2$Area.2.bins.log.mean)#$sign
  summary(indicnecr_in2)
#	indicnecr_in2$otu <- rownames(indicnecr_in2)
#	indicnecr_in2 <- left_join(indicnecr_in2, roots.joined[,1:7]) %>% filter(p.value<.05)
#	arrange(indicnecr_in2, p.value)
	
	one.more.thing.to.try_wa1 <- one.more.thing.to.try %>%
	  dplyr::select(StateBlock, SampleID, Area.2.bins.log.mean) %>%
	  filter(SampleID %in% rownames(r_otu.rarified_named2)) %>%
	  filter(StateBlock == "WA P1") %>% distinct
	indicnecr_wa1 <- indicspecies::multipatt(r_otu.rarified_named2[one.more.thing.to.try_wa1$SampleID,], one.more.thing.to.try_wa1$Area.2.bins.log.mean)$sign
	indicnecr_wa1$otu <- rownames(indicnecr_wa1)
	indicnecr_wa1 <- left_join(indicnecr_wa1, roots.joined[,1:7]) %>% filter(p.value<.05)
	arrange(indicnecr_wa1, p.value)
	
	one.more.thing.to.try_wa2 <- one.more.thing.to.try %>%
	  dplyr::select(StateBlock, SampleID, Area.2.bins.log.mean) %>%
	  filter(SampleID %in% rownames(r_otu.rarified_named2)) %>%
	  filter(StateBlock == "WA P2") %>% distinct
	indicnecr_wa2 <- indicspecies::multipatt(r_otu.rarified_named2[one.more.thing.to.try_wa2$SampleID,], one.more.thing.to.try_wa2$Area.2.bins.log.mean)$sign
	indicnecr_wa2$otu <- rownames(indicnecr_wa2)
	indicnecr_wa2 <- left_join(indicnecr_wa2, roots.joined[,1:7]) %>% filter(p.value<.05)
	arrange(indicnecr_wa2, p.value)
	
	one.more.thing.to.try_wa3 <- one.more.thing.to.try %>%
	  dplyr::select(StateBlock, SampleID, Area.2.bins.log.mean) %>%
	  filter(SampleID %in% rownames(r_otu.rarified_named2)) %>%
	  filter(StateBlock == "WA P3") %>% distinct
	indicnecr_wa3 <- indicspecies::multipatt(r_otu.rarified_named2[one.more.thing.to.try_wa3$SampleID,], one.more.thing.to.try_wa3$Area.2.bins.log.mean)$sign
	indicnecr_wa3$otu <- rownames(indicnecr_wa3)
	indicnecr_wa3 <- left_join(indicnecr_wa3, roots.joined[,1:7]) %>% filter(p.value<.05)
	arrange(indicnecr_wa3, p.value)
	
	one.more.thing.to.try_wa <- one.more.thing.to.try %>%
	  dplyr::select(State, SampleID, Area.2.bins.log.mean) %>%
	  filter(SampleID %in% rownames(r_otu.rarified_named2)) %>%
	  filter(State == "WA") %>% distinct
	indicnecr_wa <- indicspecies::multipatt(r_otu.rarified_named2[one.more.thing.to.try_wa$SampleID,], one.more.thing.to.try_wa$Area.2.bins.log.mean)$sign
	indicnecr_wa$otu <- rownames(indicnecr_wa)
	indicnecr_wa <- left_join(indicnecr_wa, roots.joined[,1:7]) #%>% filter(p.value<.05)
	arrange(indicnecr_wa, p.value)
	
	one.more.thing.to.try_in <- one.more.thing.to.try %>%
	  dplyr::select(State, SampleID, Area.2.bins.log.mean) %>%
	  filter(SampleID %in% rownames(r_otu.rarified_named2)) %>%
	  filter(State == "IN") %>% distinct
	indicnecr_in <- indicspecies::multipatt(r_otu.rarified_named2[one.more.thing.to.try_in$SampleID,], one.more.thing.to.try_in$Area.2.bins.log.mean)$sign
	indicnecr_in$otu <- rownames(indicnecr_in)
	indicnecr_in <- left_join(indicnecr_in, roots.joined[,1:7]) #%>% filter(p.value<.05)
	arrange(indicnecr_in, p.value)
	
	#windows(12,5);
	
	gg4 <- ggplot()+
	  geom_point(data=one.more.thing.to.try %>% filter(Inoculation == 'Gm'),
	             aes(x=Axis.1, y=Axis.2,
	                 colour=Area.2.bins.log.mean,
	  #               pch=StateBlock,
	                 fill=Area.2.bins.log.mean),
	             size=3) +
	  ggplot2::stat_ellipse(data=one.more.thing.to.try, aes(x=Axis.1, y=Axis.2,
	                                                colour=Area.2.bins.log.mean),
	     #                                           linetype=StateBlock),
	                        geom="polygon", fill=NA, alpha=0.2, type="t", level=.48) +

	  coord_cartesian()+#xlim=c(-.5,.6), ylim=c(-.4,.5)) +
#	  scale_shape_manual(values=21:25)+
	  xlab("\nPCOA1\n")+ylab("PCOA2\n")+
	  guides(colour=guide_legend(title=paste0("Seedling canker size<br><br>",
	                                          "<span style='font-size: 8pt'>",
	                                          "Plant avg. necrotic area<br>compared to seedlings in plot<br>")))+
	                                          
	                                  #        "Mean necrotic area compared<br>to other seedlings in the plot<br>")))+
#	  guides(pch=guide_legend(title="State & Block"))+
#	  guides(linetype=guide_legend(title="State & Block"))+
	  guides(fill=guide_legend(title=paste0("Seedling canker size<br><br>",
	                                        "<span style='font-size: 8pt'>",
	                                        "Plant avg. necrotic area<br>compared to seedlings in plot<br>")))+
	  theme_bw() + theme(legend.title = ggtext::element_markdown()) +
	  facet_grid(tissue~StateBlock)
	  
	windows();gg4
	
	#ggplot_build(gg3) %>% str
	
	windows();ggplot()+
	  geom_point(data = one.more.thing.to.try,
	             aes(x = one.more.thing.to.try$Axis.1,
	                 y = one.more.thing.to.try$Area.sd))+
	  facet_grid(tissue~StateBlock)
	
	windows(12,5)
	gg3
	
	windows(10,5)
	gg2

	windows(10,5);ggca_all+guides(colour=guide_legend(title="Eigenvalue"))+labs(x="\nCCA1",y="CCA2\n")
##########

##### COMPOSITE
	
#	windows(10,8)
	#grid.arrange(gg2, ggca2, ggca1 + theme(axis.title.y=element_blank()),
	#             layout_matrix=rbind(c(1,1),c(2,3)))
	#ggarrange(gg2, ggca_all+labs(x="\nCCA1",y="CCA2\n"), nrow=2)#, common.legend=T, legend="right")

#	ggarrange(
#	  gg2,
#	  ggca_all,
#	  nrow=2, common.legend=T, legend="right",
#	  legend.grob = grid.arrange(get_legend(gg2), get_legend(ggca_all)))
	windows(10,8)
	library(devEMF)
	library(Cairo)
	#setwd('Revision_Jan_2025')
	#emf(file = "Fig_composite_new.emf",10,8,emfPlus=T)
	svg("Figure_3_March2025_3.svg", 10, 8)
	grid.arrange(
	  egg::tag_facet(gg2+scale_y_continuous(labels=scaleFUN)+theme(legend.position="none"), tag_pool=c("A","B"), open="", close="", size = 6),
	  egg::tag_facet(ggca_all+scale_y_continuous(labels=scaleFUN)+theme(legend.position="none"), tag_pool=c("C","D"), open="", close="", size = 6),
	  get_legend(gg2),
	  get_legend(ggca_all),
	  ggplot()+ggtitle(label="Cankers")+theme_void()+theme(plot.title.position="panel", plot.title = element_text(hjust=.6, face="bold", size=16)),
	  ggplot()+ggtitle(label="Roots")+theme_void()+theme(plot.title.position="panel", plot.title = element_text(hjust=.6, face="bold", size=16)),
	  NULL,
	  widths=c(4.25,4.25,1.5),
	  heights=c(.5,4.75,4.75),
	  layout_matrix=rbind(c(5,6,7),c(1,1,3),c(2,2,4)))
	dev.off()
	 # legend.grob = grid.arrange(get_legend(gg2), get_legend(ggca_all)))
			
###################################

	
	###################################
	###################################
	
	
	###################################
	
	#root.pcoa <- pcoa(d.all)
	#root.pcoa$vectors
	#with(root.pcoa, plot(vectors[,1], vectors[,2], col=r.gm, pch=as.numeric(r.block)))
	#with(root.pcoa, plot(vectors[,2], vectors[,3], col=r.block))

	c.ordi.summary.data <- cbind(cank.nmds$points[,1:3], state = as.character(c.state), block = as.character(c.block), gm = as.character(c.gm)) %>% as.data.frame
	c.ordi.summary.data$state.col <- as.numeric(c.state)
	c.ordi.summary.data$gm.symbol <- as.numeric(c.gm)
	c.ordi.summary.data$block2    <- c.ordi.summary.data$block#; levels(c.ordi.summary.data$block2) <- c("Plot 1", "Plot 2", "Plot 1", "Plot 2", "Plot 3")
	c.ordi.summary.data$state.block  <- with(c.ordi.summary.data, as.numeric(as.factor(paste(state,block2))))+20
	c.ordi.summary.data$state.block.gm <- with(c.ordi.summary.data, paste(state, block2, gm))
	c.ordi.summary.data$state.gm <- with(c.ordi.summary.data, paste(state, gm)) %>% as.factor %>% as.numeric

	c.toplot <- summaryBy(MDS1 + MDS2 + MDS3 ~ block + gm, data = c.ordi.summary.data, FUN = c(mean, function(x) mean(x)+ sd(x)/ sqrt(length(x)), function(x) mean(x)-sd(x)/sqrt(length(x))), fun.names=c("mean","plus","minus"), id=c("state.block", "state.block.gm", "gm.symbol", "state.gm"))

	r.ordi.summary.data <- cbind(root.nmds$points[,1:3], state = as.character(r.state), block = as.character(r.block), gm = as.character(r.gm)) %>% as.data.frame
	r.ordi.summary.data$state.col <- as.numeric(r.state)
	r.ordi.summary.data$gm.symbol <- as.numeric(r.gm)
	r.ordi.summary.data$block2    <- r.ordi.summary.data$block; #levels(r.ordi.summary.data$block2) <- c("Plot 1", "Plot 2", "Plot 1", "Plot 2", "Plot 3")
	r.ordi.summary.data$state.block  <- with(r.ordi.summary.data, as.numeric(as.factor(paste(state,block2))))+20
	r.ordi.summary.data$state.block.gm <- with(r.ordi.summary.data, paste(state, block2, gm))
	r.ordi.summary.data$state.gm <- with(r.ordi.summary.data, paste(state, gm)) %>% as.factor %>% as.numeric
str(r.ordi.summary.data)
r.ordi.summary.data$MDS1 <- as.numeric(r.ordi.summary.data$MDS1)
r.ordi.summary.data$MDS2 <- as.numeric(r.ordi.summary.data$MDS2)
r.ordi.summary.data$MDS3 <- as.numeric(r.ordi.summary.data$MDS3)

	r.toplot <- summaryBy(MDS1 + MDS2 + MDS3 ~ block + gm, data = r.ordi.summary.data, FUN = c(mean, function(x) mean(x)+ sd(x)/ sqrt(length(x)), function(x) mean(x)-sd(x)/sqrt(length(x))), fun.names=c("mean","plus","minus"), id=c("state.block", "state.block.gm", "gm.symbol", "state.gm"))

#######################
##### FIG 4 ###########
#######################

windows(width=8, height=6); layout(mat=rbind(c(1,1,2,2,5),c(1,1,2,2,5),c(3,3,4,4,5),c(3,3,4,4,5))); par(mar=c(2,5,4,1), xpd=T)

	with(r.toplot, plot(MDS1.mean, MDS2.mean, col = state.gm, bg= state.gm, pch = state.block, xlim=c(min(MDS1.minus),max(MDS1.plus)), ylim=c(min(MDS2.minus),max(MDS2.plus)), xlab = "", ylab = "NMDS 2", cex=2))
	with(r.toplot, arrows(x0 = MDS1.mean, y0 = MDS2.minus, y1 = MDS2.plus, length=.1, angle=90, code=3, col = state.gm))
	with(r.toplot, arrows(x0 = MDS1.minus, x1 = MDS1.plus, y0 = MDS2.mean, length=.1, angle=90, code=3, col = state.gm))
	with(r.toplot,text(min(MDS1.minus)+1,max(MDS2.plus)-1,"A"))
	par(mar=c(2,2,4,4), xpd=T)
	with(r.toplot, plot(MDS3.mean, MDS2.mean, col = state.gm, bg= state.gm, pch = state.block, xlim=c(min(MDS3.minus),max(MDS3.plus)), ylim=c(min(MDS2.minus),max(MDS2.plus)), xlab = "", ylab = "", cex=2, axes=F)); box(); axis(1)
	with(r.toplot, arrows(x0 = MDS3.mean, y0 = MDS2.minus, y1 = MDS2.plus, length=.1, angle=90, code=3, col = state.gm))
	with(r.toplot, arrows(x0 = MDS3.minus, x1 = MDS3.plus, y0 = MDS2.mean, length=.1, angle=90, code=3, col = state.gm))
	with(r.toplot,text(min(MDS3.minus)+1,max(MDS2.plus)-1,"A"))
# canks
# windows(width=12, height=4); par(mfrow=c(1,3), mar=c(5,5,1,1), oma=c(1,1,1,1), xpd=T)
par(mar=c(5,5,1,1), xpd=T)

	with(c.toplot, plot(MDS1.mean, MDS2.mean, col = state.gm, bg= state.gm, pch = state.block, xlim=c(min(MDS1.minus),max(MDS1.plus)), ylim=c(min(MDS2.minus),max(MDS2.plus)), xlab = "NMDS 1", ylab = "NMDS 2", cex=2))
	with(c.toplot, arrows(x0 = MDS1.mean, y0 = MDS2.minus, y1 = MDS2.plus, length=.1, angle=90, code=3, col = state.gm))
	with(c.toplot, arrows(x0 = MDS1.minus, x1 = MDS1.plus, y0 = MDS2.mean, length=.1, angle=90, code=3, col = state.gm))
	with(c.toplot,text(min(MDS1.minus)+1,max(MDS2.plus)-1,"B"))
	par(mar=c(5,2,1,4), xpd=T)
	with(c.toplot, plot(MDS3.mean, MDS2.mean, col = state.gm, bg= state.gm, pch = state.block, xlim=c(min(MDS3.minus),max(MDS3.plus)), ylim=c(min(MDS2.minus),max(MDS2.plus)), xlab = "NMDS 3", ylab = "", axes=F, cex=2)); box(); axis(1)
	with(c.toplot, arrows(x0 = MDS3.mean, y0 = MDS2.minus, y1 = MDS2.plus, length=.1, angle=90, code=3, col = state.gm))
	with(c.toplot, arrows(x0 = MDS3.minus, x1 = MDS3.plus, y0 = MDS2.mean, length=.1, angle=90, code=3, col = state.gm))
	with(c.toplot,text(min(MDS3.minus)+1,max(MDS2.plus)-1,"B"))
	par(mar=c(0,0,0,0), xpd=T);plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));with(c.toplot, legend("left", legend=state.block.gm, col= state.gm[c(1,2,4,3,6,5,8,7)], pch=state.block[c(1,2,4,3,6,5,8,7)], pt.bg= state.gm[c(1,2,4,3,6,5,8,7)], bty="n", cex=1.2))



##############################
## Figs. 6 & 8              ##
## Heatmaps.                ##
##############################


########## Roots

	names(roots.joined)<-sub("\\-","",names(roots.joined),perl=T)%>% sub("\\-","",.,perl=T)
	roots.joined.fam<-summaryBy(.~Family, data=roots.joined[,c(5,8:82)], FUN=sum, keep.names=T)
	dim(roots.joined.fam)
	roots.joined.t <- roots.joined.fam[,2:76] %>% t %>% data.frame
	roots.joined.t$SampleID <- rownames(roots.joined.t)
	necr.summary$SampleID <- sub("\\-","",necr.summary$SampleID,perl=T)%>% sub("\\-","",.,perl=T)
	matrix.roots<-right_join(roots.joined.t, necr.summary[,c("SampleID","Area.mean")], by ="SampleID")[,-25]
	names(matrix.roots)[-25] <-roots.joined.fam[,1]
	matrix.roots$Area.mean <- log(matrix.roots$Area.mean)
	
	#matrix.roots[-25]<-matrix.roots[-25]+.0000001%>%log
	
	library(Hmisc)
	corr.mat <- rcorr(as.matrix(matrix.roots), type="spearman")
	
#	corr.mat$r%>%as.matrix%>%as.data.frame%>%write.csv("roots.correlations.csv", row.names=T)
	#nas<-which(is.na(corr.mat$P)[,1])
	
	#corr.mat$r <- corr.mat$r[-nas,-nas]
	#corr.mat$P <- corr.mat$P[-nas,-nas]
	#corr.mat$n <- corr.mat$n[-nas,-nas]
	
	#dim(corr.mat$r)
	ticks<-(0:(dim(corr.mat$r)[1]-1))/(dim(corr.mat$r)[1]-1)
	length(ticks)
	tick.names<-dimnames(corr.mat$r)[[1]]
	length(tick.names)
	tick.names[25]<-"Mean necrotic area"
	
	corr.mat$r
	
	image.plot<-corr.mat$r
	image.plot[upper.tri(corr.mat$r)]<- NA
	image.plot2<-(corr.mat$P < 0.05)
	image.plot2[lower.tri(corr.mat$r)]<-NA

# Fig 6

	windows(height=6,width=8);par(mar=c(8,2,2,16))



	image(z=image.plot, zlim= c(min(corr.mat$r),max(corr.mat$r)), axes=F, col=topo.colors(12) );axis(1,at=ticks,labels= tick.names,las=2,cex.axis=.75, font=1);axis(4,at=ticks,labels=tick.names,las=2,cex.axis=.75, font=1)
		image(image.plot2,axes=F, add=T, col=c("white","red"))
		l<-seq(min(corr.mat$r),max(corr.mat$r),length.out=12)%>% round(2) %>% c("p < 0.05")
		legend(1.4,.75, fill=c(topo.colors(12),"red"), legend=l, bty='n', xpd=T)
		grid(length(tick.names))
	#image(image.plot,axes=F);axis(1,at=ticks,labels= tick.names,las=2,cex.axis=.5);axis(4,at=ticks,labels=tick.names,las=2,cex.axis=.75)

### Cankers

	names(canks.joined)<-sub("\\-","",names(canks.joined),perl=T)%>% sub("\\-","",.,perl=T)
	dim(canks.joined)
	canks.joined.gen<-summaryBy(.~Genus, data= canks.joined[,c(6,8:82)], FUN=sum, keep.names=T)
	canks.joined.t <- canks.joined.gen[,2:76] %>% t %>% data.frame
	canks.joined.t$SampleID <- rownames(canks.joined.t)

	matrix.canks<-right_join(canks.joined.t, necr.summary[,c("SampleID","Area.mean")], by ="SampleID")[,-14]
	names(matrix.canks)[-14] <-canks.joined.gen[,1]
	matrix.canks $Area.mean <- log(matrix.canks $Area.mean)
	matrix.canks<-matrix.canks[,-which(colSums(matrix.canks)==0)]
	
	matrix.canks[,-14]<- matrix.canks[,-14]+.0000001 %>% log
	
	cank.corr.mat <- rcorr(as.matrix(matrix.canks), type="pearson")
#	cank.corr.mat$r%>%as.matrix%>%as.data.frame%>%write.csv("cankers.correlations.csv", row.names=T)
	
	#dim(cank.corr.mat $r)
	ticks<-(0:(dim(cank.corr.mat $r)[1]-1))/(dim(cank.corr.mat $r)[1]-1)
	length(ticks)
	tick.names<-dimnames(cank.corr.mat $r)[[1]]
	length(tick.names)
		tick.names[12]<-"Mean necrotic area"

	cank.corr.mat $P
	
	image.plot<-cank.corr.mat $r
	image.plot[upper.tri(cank.corr.mat $r)]<- NA
	image.plot2<-(cank.corr.mat $P < 0.1)
	image.plot2[lower.tri(cank.corr.mat $r)]<-NA

	windows(height=6,width=8);par(mar=c(8,2,2,16))

# Fig 7

	image(z=image.plot, zlim= c(min(cank.corr.mat $r),max(cank.corr.mat $r)), axes=F, col=topo.colors(12) );axis(1,at=ticks,labels= tick.names,las=2,cex.axis=.75,font=3);axis(4,at=ticks,labels=tick.names,las=2,cex.axis=.75,font=3)
		image(image.plot2,axes=F, add=T, col=c("white","red"))
		l<-seq(min(cank.corr.mat $r),max(cank.corr.mat $r),length.out=12)%>% round(2) %>% c("p < 0.1")
		legend(1.4,.75, fill=c(topo.colors(12),"red"), legend=l, bty='n', xpd=T)
		grid(length(tick.names))

# fig 1
		
		windows();
		
		fig1a <-
		  
		ggplot(data = necrosis %>%
		         left_join(distinct(dplyr::select(one.more.thing.to.try, State, Block, StateBlock)))%>%
		         
		         #mutate(StateBlock = paste(State, Block)) %>%
		         filter(Area < 150))+
		  geom_boxplot(
		    mapping = aes(
		      y = StateBlock,
		      x = Area,
		      fill = State
		    ), outliers =F) +
		  geom_jitter(mapping = aes(
		    y = StateBlock,
		    x = Area
		  ), height = .2) +
		  #guides("State & Block")
		  labs(
		    x = "\nNecrotic area (mm\u00b2)",
		    y = "Seedling germination and\ngrowth plot (year 1)\n"
		  )+
		  theme_bw() +
		  theme(panel.grid = element_blank())
		
		letterlayer <- ggplot_build(fig1a)$data[[1]] %>% cbind(let = c('a','ab','b','ab','ab'))
		
		windows(4.5,3.5); fig1a + geom_text(data = letterlayer,
		                                    mapping = aes(y=y, x=xmax+10, label=let),
		                                    vjust=0.25, size=5)
		
		fig1a <- fig1a + geom_text(data = letterlayer,
		                  mapping = aes(y=y, x=xmax+10, label=let),
		                  vjust=0.25, size=5)

		
		str(one.more.thing.to.try)
		
		# testing
		one.more.thing.to.try %>% filter(Inoculation == "Gm") %>%
		  filter(tissue=='Roots', StateBlock=="IN P1") %>%
		  (function(w)
		    dist(select(w, Axis.1, Axis.2)) %>% as.matrix %>%
		     .[which(as.numeric(w$Area.2.bins.log.mean)==1),
		       which(as.numeric(w$Area.2.bins.log.mean)==2)] %>% mean) 
		
		abc <- function (a,b,c) dist(cbind(a,b)) %>% as.matrix %>%
		  .[which(as.numeric(c)==1),
		    which(as.numeric(c)==2)] %>% mean
		
		one.more.thing.to.try %>% filter(Inoculation == "Gm") %>%
		  filter(tissue=='Roots', StateBlock=="IN P1") %>%
		  with(., abc(Axis.1, Axis.2, Area.2.bins.log.mean))
		
		
		sum(necrosis$Area>=150)
		sum(necrosis$Area>=200)
		necrosis[which(necrosis$Area>=200),]
		
		# calculate the distance between the big cankers and the little cankers within each group
		
		one.more.thing.to.try2 <-
		one.more.thing.to.try %>% filter(Inoculation == "Gm") %>%
		  #filter(Area.mean < 200)%>%#150) %>%
		  filter(Area.mean < 150) %>%
		  dplyr::group_by(tissue, State, Block, StateBlock) %>%
		  dplyr::summarise(
		    microbial_distance_necrosis = abc(Axis.1, Axis.2, Area.2.bins.log.mean)) %>%
		  
		  left_join(necrosis %>% filter(Inoculation == "Gm") %>%
		              #filter(Area < 200)%>%#150) %>%
		              filter(Area < 150) %>%
		            group_by(State, Block, SampleID) %>% dplyr::summarise(necrosis_plant_mean = mean(Area)) %>%
		              group_by(State, Block) %>% dplyr::summarise (necrosis_var = var(necrosis_plant_mean)))
		
		             #group_by(State, Block) %>% summarise (necrosis_var = var(Area)))
		
		

		  
		  #,
		  #  necrosis_var = (function (m, v) var(m)#+mean(v^2)
		  #                  )(Area.mean, Area.sd)) #%>%
		  
		  #windows();
		  
		  rootshootdata<-read.csv("Revision_Jan_2025/downloaded-purdue-sharepoint/cank_rootshoot_int_params.csv")%>%
		    with(
		      data.frame(
		        Param = Parameter[7:10] %>% gsub('^treat','',.) %>% sub('\\:cov$','',.) %>% c('IN P1',.),
		        Slope = Coefficient[6]+c(0,Coefficient[7:10])
		      )
		    ) %>% left_join(one.more.thing.to.try2, by = join_by('Param'=='StateBlock')) 
		    
		  
		  r1 <- range(rootshootdata$Slope) %>% (function (x) x[2]-x[1])
		  r2 <- range(rootshootdata$necrosis_var) %>% (function (x) x[2]-x[1])

		  a <- range(rootshootdata$Slope)[1]
		  b <- range(rootshootdata$necrosis_var)[1]
		  
		  center_necr <- (1-range(rootshootdata$Slope)[2]/r1)*r2+b
		  
		  m <- r2/r1
		  
		  rootshootdata_melted <- rootshootdata %>%
		    mutate(slope_transformed = center_necr + Slope*m)%>%
		    reshape2::melt(measure.vars=c('slope_transformed', 'necrosis_var'))
		  
		  levels(rootshootdata_melted$variable) <-c("Slope coefficient","Necrotic area variance")
		    
		  
		  position_scales_x <-
		    list(
		      scale_x_continuous(limits = c(0.34,.5), breaks = c(.36,.42,.48)),
		      scale_x_continuous(limits = c(.15,.45),breaks = c(.2,.3,.4)))
		  position_scales_y <-
		    list(
		      scale_y_continuous(limits = c(300,500), breaks = c(320,380,440,500)),
		      scale_y_continuous(limits = c(300,1200), breaks= c(300,600,900,1200)))
		  
		  position_scales_y2 <-
		    list(
		      scale_y_continuous(name = "Among-plant variance in\navg. necrotic area (mm\u2074)\n", limits = c(300,700), breaks = (3:7)*100,
		                         sec.axis = sec_axis(~ (.-center_necr)/m, breaks = c(-4:0)/25+.02,
		                                             name = "Slope coefficient\nroot-shoot ratio to necrotic area")),
		      scale_y_continuous(name = "Among-plant variance in\navg. necrotic area (mm\u2074)\n", limits = c(175,1200), breaks= c(200,400,600,800,1000,1200),
		                         sec.axis = sec_axis(~ (.-center_necr)/m, breaks = c(-2:2)/10,
		                                             name = "Slope coefficient\nroot-shoot ratio to necrotic area")))
		  
		  fig1b_alt <-
		  ggplot(rootshootdata_melted) +
		    geom_point(aes(x = microbial_distance_necrosis ,
		                   y = value,
		                   color = variable,
		                   pch = Param
		    ), size=2.5)+
		      scale_y_continuous(#name = "Variance average necrotic area among plants\n",
		                         sec.axis = sec_axis(~ (.-center_necr)/m))+#, name = "Slope coefficient of root-shoot ratio and necrotic area\n"))+
		    scale_color_manual(name = 'Y Axis', limits = c("Slope coefficient","Necrotic area variance"), values=c('darkgreen','darkviolet'))+
		    scale_shape_discrete(name = "State & Block")+
		    labs(x = "\nCommunity difference between seedlings\nwith large and small cankers") +
		    facet_grid(State ~ tissue, scales="free")+
		    facetted_pos_scales(x = position_scales_x,
		                        y = position_scales_y2)+#,
		    #y= position_scales_y2)+
		    theme_bw()+theme(panel.grid = element_blank(),
		                     axis.text.x = element_text(angle=90,vjust=.5),
		                     strip.background = element_blank(),
		                     strip.placement = 'outside',
		                     axis.title.y.left = element_text(colour='darkviolet'),
		                     axis.text.y.left  = element_text(colour='darkviolet'),
		                     axis.title.y.right = element_text(colour='darkgreen'),
		                     axis.text.y.right  = element_text(colour='darkgreen'))
		  
		  fig1b <-
		  ggplot(one.more.thing.to.try2) +
		    geom_point(aes(x = microbial_distance_necrosis ,
		                   y = necrosis_var,
		                 #  color = Block,
		                   pch = StateBlock
		                 ), size=2.5) +
		    #geom_point(data = rootshootdata, mapping = aes(x = microbial_distance_necrosis, y = Slope, pch=StateBlock))+
		    labs(
		      x = "\nCommunity difference between seedlings\nwith large and small cankers",
		      y = "Variance in individual\nplant necrotic area means\n") +
		    facet_grid(State ~ tissue, scales="free")+
		    facetted_pos_scales(x = position_scales_x, y= position_scales_y)+
		    
		    theme_bw()+theme(panel.grid = element_blank(),
		                     axis.text.x = element_text(angle=90,vjust=.5),
		                     strip.background = element_blank(),
		                     strip.placement = 'inside')
		
		  windows(4.5,3.5); fig1b
		  
		  
		  # and one analysis

		with(one.more.thing.to.try2,
		     lmer(necrosis_var ~ microbial_distance_necrosis*(1|tissue)*(1|State)) %>% 
		       #summary())
		       car::Anova(type=2))
		with(one.more.thing.to.try2,
		     lmer(necrosis_var ~ microbial_distance_necrosis*(1|tissue)*(1|State)) %>% 
		       summary())
		       #car::Anova(type=3))
		with(one.more.thing.to.try2,
		     lmer(necrosis_var ~ microbial_distance_necrosis*(1|tissue)*(1|State))) %>% 
		       #summary())
		       #car::Anova(type=2))
		       simulateResiduals %>% plot
		
		#windows(9.25,7)
		devEMF::emf("Fig1_new.emf", emfPlus = T, 10,6.5)
		#svg('Revision_Jan_2025/Figures_ordered/Fig1ABC.svg', 10,6.5)
		#windows(10,6.5)
		((free(fig1a, side='b') + fig1b_alt)/
		    #free(
		      free(gg4+theme_bw()+
		        theme(axis.text = element_blank(),
		              axis.ticks = element_blank(),
		              strip.background = element_blank(),
		              legend.title = ggtext::element_markdown())))+
		  plot_layout(heights=c(1,1.6), tag_level='new')+
		  plot_annotation(tag_levels = 'A') & theme(strip.text = element_text(face='bold'))
		dev.off()
		
		read.csv ("../WA-IN-Transplant2019-20.Necrosis.csv") %>% 
		  mutate(SampleID = paste(State, Block, Plant, sep="-")) %>%
		  left_join(design[,4:5], by="SampleID") %>%
		  with(lm(log1p(Area) ~ Inoculation)%>% with(., list(summary(.),Anova(., type=3))))
		
		read.csv ("../WA-IN-Transplant2019-20.Necrosis.csv") %>% 
		  mutate(SampleID = paste(State, Block, Plant, sep="-")) %>%
		           left_join(design[,4:5], by="SampleID") %>%
		  
		  ggplot()+geom_density(aes(x=Area, fill=Inoculation), bw='nrd', alpha=.5)+
		  scale_x_continuous(limits = c(0,150))+
		  labs(x =  "\nNecrotic area (cm\u00b2)", y = "Density\n")+theme_bw()
		  
		  filter(Inoculation == "Sham")
		