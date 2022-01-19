#!/usr/bin/env Rscript
#
# Use this filr along with depth file
# samtools depth -d 0 -aa $file.bam > $file.depth 
# awk 'BEGIN{print "Location \t Depth";}{print $(NF-1)"\t"$NF}' $file.depth > depth.$$
# Rscript  /home4/nCov/Sreenu/Scripts/getCoverage.R depth.$$ $file
#
#Developed by: Dr. Sreenu Vattipally

library(ggplot2)

args <- (commandArgs(TRUE))


ttlName<-paste0("Coverage of ",args[2]);
plotName<-paste0(args[2], "-coverage.pdf");

pdf(file=plotName, height=7, width=15)

plotTheme<-theme( 
		 plot.title = element_text(color="#705AA5", size=14, face="bold.italic",hjust = 0.5), 
		 axis.title.x = element_text(color="#5AA570", size=14, face="bold"), 
		 axis.title.y = element_text(color="#A5705A", size=14, face="bold"), 
		 axis.text.x = element_text(face="bold",color="#5AA570",size=10), 
		 axis.text.y = element_text(face="bold",color="#A5705A", size=10))

data<-read.table(args[1],header = FALSE);

gplot<-ggplot(data,aes(x=V2,y=V3))+geom_line(color="#688197")
gplot<-gplot+geom_hline(yintercept = mean(data$V3),color="#A5705A")
gplot<-gplot+ggtitle(ttlName) + labs(x = "Location", y="Depth")
gplot+theme_light() + plotTheme

# Fush output to PDF
dev.off()
