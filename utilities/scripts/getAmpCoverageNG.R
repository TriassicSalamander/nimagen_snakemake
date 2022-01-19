#!/usr/bin/env Rscript
#
# Use this filr along with depth file
# samtools depth -d 0 -aa $file.bam > $file.depth 
# /home4/nCov/Sreenu/Scripts/ampliconDepth.awk /home4/nCov/Sreenu/Primers/nCoV-2019-V3.bed  $file.depth > $file.amplicon.cov
#
#Developed by: Dr. Sreenu Vattipally

library(ggplot2)

args <- (commandArgs(TRUE))

ttlName<-paste0("Amplicon depth of ",args[2]);
plotName<-paste0(args[2], "-Amplicon-Depth.pdf");

pdf(file=plotName, height=5, width=15)

plotTheme<-theme( 
		 plot.title = element_text(color="#705AA5", size=14, face="bold.italic",hjust = 0.5), 
		 axis.title.x = element_text(color="#5AA570", size=14, face="bold"), 
		 axis.title.y = element_text(color="#A5705A", size=14, face="bold"), 
		 axis.text.x = element_text(face="bold",color="#5AA570",size=7), 
		 axis.text.y = element_text(face="bold",color="#A5705A", size=10))

data<-read.table(args[1],header = FALSE);
gplot<-ggplot(data,aes(x=V1,y=V4))+geom_col()
gplot<-gplot+ggtitle(ttlName) + labs(x = "Amplicon number", y="Average Depth")
gplot+theme_light() + plotTheme + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
#gplot+theme_light() + plotTheme + scale_x_continuous(breaks = seq(1,154, by=5))


# Fush output to PDF
dev.off()
