#!/usr/bin/env Rscript

###
# Plot Ancestry Lineages
###

#takes in an output from CW's ancestry plotting 
#as of right now, only takes a single chromosome for 1-30 samples

#to do:
	#Generalize a bit more....
	#add in errors for missing packages
	#colour by scaffold?
	
#usage: Rscript ancestryPlots.R <FILE1 .ahmm.maxpost> <chromosome>
#e.g: Rscript ancestryPlots.R nofilter.16.Ne5k.ahmm.maxpost  16

#Load packages --------------------
require(ggplot2)
require(ggthemes)
require(plyr)
require(wesanderson)
require(reshape)
require(scales)
require(grid)

#load data.frames --------------
args = commandArgs(trailingOnly=TRUE)
anc =  read.table(file= args[1] ,header=F, skip=1)
chr = args[2]
samps = ncol(anc) -1
names(anc) = c("POS", paste("SAMPLE",c(1:samps),sep=""))
scaffs = read.table(file="../resources/scaffold_chr_map_4.5.txt",header=T)

#anc$relPOS = seq(1,nrow(anc),1)
#anc$POS = NULL
#anc.lineage = 1-anc
#sudo Rscript ancestryPlots.R /home/charlie/PooledAfz/pooledahmm_data/nofilter.16.Ne5kHiRes.ahmm.maxpost 16

#Calculate MN and SD for lineage -------------------------
mean.anc = apply(anc[-1],1 , mean,na.rm=T)
sd.anc = apply(anc[-1],1 , sd,na.rm=T)

#create dataframes -------------------------------------
mean.anc.df = data.frame(cbind(mn = mean.anc,sd = sd.anc,pos = anc$POS))
scaffs = scaffs[scaffs$chromosome==chr,]



# plot --------------------------------


mean.anc.plot = ggplot(mean.anc.df,aes(x=pos,y=mn)) +
				geom_rect(data = scaffs, 
							aes(xmin = start, xmax = stop, ymin = 0, ymax = 0.65) , #add in rectangles for scaffolds
							fill="grey", alpha=0.5, inherit.aes = FALSE) +
				geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), color = wes_palette("GrandBudapest")[1]) +				
				geom_point(aes(x=pos,y=mn), color = wes_palette("GrandBudapest")[3],  size = 2.3) +
				
				theme_bw() +
				scale_x_continuous(labels = comma) +
				scale_y_continuous(limits = c(0, .65)) +
					theme(strip.background = element_blank(),
							axis.line.x = element_line(size = 1, colour = "black"),
							axis.line.y = element_line(size = 1, colour = "black"),
							axis.line = element_line(size = 1, colour = "black"),
							text = element_text(size=18),
							axis.ticks = element_line(size = 1), 
							panel.grid.major = element_blank(),
							panel.grid.minor = element_blank(),
							panel.background = element_blank(), 
							panel.border = element_blank()
							) +						
				labs(x = "Position (bp)", y = "Proportion M Ancestry") 

				
				
				
#mean.anc.plot



ggsave(filename = paste(args[1], ".png",sep=""), mean.anc.plot, width = 12, height = 6, units = "in" )	











#
#
##rotate and assemble df for plotting --------------
#M.anc = melt(anc, id=c("relPOS"))
#M.anc$Ancestry = rep("M", nrow(M.anc))
#A.anc = melt(anc.lineage , id=c("relPOS"))
#A.anc$Ancestry = rep("A", nrow(A.anc))
#A.anc$relPOS = M.anc$relPOS
#AM.anc = data.frame(rbind(M.anc, A.anc))
#
#
## plot --------------------------------
#pal <- wes_palette("Chevalier",2)
#AM.anc = AM.anc[AM.anc$variable=="SAMPLE1" | AM.anc$variable=="SAMPLE2" | AM.anc$variable=="SAMPLE10" | AM.anc$variable=="SAMPLE16",]
#
#anc.plot = ggplot(AM.anc,aes(x=relPOS,y=value,fill=Ancestry)) +
#				facet_wrap( ~ variable, ncol = 1) +
#				geom_bar(stat="identity") +
#				theme_bw() +
#				scale_fill_manual(values=pal) +
#				scale_x_continuous(labels = comma) +
#					theme(strip.background = element_blank(),
#							strip.text = element_text(size = rel(3.0), vjust = -4.0), 
#							panel.margin.y = unit(-0.01, "lines"),
#							axis.text.y = element_blank(),
#							axis.ticks.y=element_blank(),
#							axis.ticks.y=element_blank(),
#							strip.text.x = element_blank(),
#							panel.grid.major = element_blank(),
#							panel.grid.minor = element_blank(),
#							panel.background = element_blank(), 
#							panel.border = element_blank(),
#							legend.position="bottom",
#							legend.title = element_blank()
#							) +
#				labs(x = "Position (bp)", y = "")
#anc.plot
#
#
##rotate and assemble df for plotting --------------
#M.anc = melt(anc, id=c("POS"))
#M.anc$Ancestry = rep("M", nrow(M.anc))
#A.anc = melt(anc.lineage , id=c("POS"))
#A.anc$Ancestry = rep("A", nrow(A.anc))
#A.anc$POS = M.anc$POS
#AM.anc = data.frame(rbind(M.anc, A.anc))
#
#
#
#
## plot --------------------------------
#pal <- wes_palette("Chevalier",2)
#AM.anc = AM.anc[AM.anc$variable=="SAMPLE1" | AM.anc$variable=="SAMPLE2" | AM.anc$variable=="SAMPLE10" | AM.anc$variable=="SAMPLE16",]
#
#anc.plot = ggplot(AM.anc,aes(x=POS,y=value,fill=Ancestry)) +
#				facet_wrap( ~ variable, ncol = 1) +
#				geom_bar(stat="identity") +
#				theme_bw() +
#				scale_fill_manual(values=pal) +
#				scale_x_continuous(labels = comma) +
#					theme(strip.background = element_blank(),
#							strip.text = element_text(size = rel(3.0), vjust = -4.0), 
#							panel.margin.y = unit(-1, "lines"),
#							axis.text.y = element_blank(),
#							axis.ticks.y=element_blank(),
#							axis.ticks.y=element_blank(),
#							strip.text.x = element_blank(),
#							panel.grid.major = element_blank(),
#							panel.grid.minor = element_blank(),
#							panel.background = element_blank(), 
#							panel.border = element_blank(),
#							legend.position="bottom",
#							legend.title = element_blank()
#							) +
#				labs(x = "Position (bp)", y = "")
#anc.plot
#
#



