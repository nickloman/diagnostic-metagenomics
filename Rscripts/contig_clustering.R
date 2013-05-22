
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(ggplot2)


blob.plot<-function(dfilt) {
  theme_set(theme_bw())
  
  colorbrewer=c("#DDDDDD", "#4575B4", "#ABD9E9", "#777777", "#313695", "#E0F3F8", "#FEE090", "#74ADD1", "#FFFFBF", "#F46D43", "#D73027", "#FDAE61", "#A50026")
  paultol    =c("#DDDDDD", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#44AA99", "#999933", "#AA4499", "#882255", "#777777")
  
  g <- ggplot() + scale_colour_manual(values=paultol, name="Taxonomic\nClassification", limits=levels(dfilt$taxon) )
  for (t in levels(dfilt$taxon)) {
    g <- g + geom_point(data=dfilt[dfilt$taxon==t,],aes(gc, cov, colour=taxon), size=2, alpha=I(1/3))
  }
  axis_breaks = c(1,2,5,10,20,50,100,200,500,1000,10000);
  g +      facet_wrap(~NSamples) +
    xlim(0,1) + ylim(0,max(dfilt$cov, na.rm=TRUE)) + scale_y_log10() +
    #breaks = axis_breaks, labels = axis_breaks) +
    labs(x="GC content", y="Read coverage") + 
    guides(colour = guide_legend(override.aes = list(alpha = 1,size=10))) + 
    opts (
      strip.text.x = theme_text(colour = "black", size = 25, vjust = 0.5),
      axis.text.x  = theme_text(colour = "black", size = 25, vjust = 1),
      axis.text.y  = theme_text(colour = "black", size = 25, vjust = 0.5),
      axis.title.x = theme_text(colour = "black", size = 25, vjust = 0),
      axis.title.y = theme_text(colour = "black", size = 25, hjust = 0.5, vjust = 0.5, angle=90),
      legend.text  = theme_text(colour = "black", size = 25, vjust = 0),
      legend.title = theme_text(colour = "black", size = 25, vjust = 0, hjust = 0, lineheight=1),
      legend.key.height = unit(2.2,"line"),
      legend.justification=c(1,0), legend.position=c(1,0)
    )
}

clean.blobs<-function(d,threshold) {
  annotated<-d[d$taxon!="Not annotated",]
  total<-dim(annotated)[1]
  levels(d$taxon)[which(table(d$taxon)<threshold*total)]<-"Not annotated"
  return(d)
}

orig <- read.table("10pc_contigstaxon.txt", header=TRUE, sep='\t')
dfilt <- clean.blobs(orig,0.005)
taxa  <- rownames(rev(sort(table(dfilt$taxon))))
taxa  <- taxa[taxa != "Not annotated"]
dfilt$taxon<-ordered(dfilt$taxon,levels=c("Not annotated",taxa)) #put the taxon category with the most hits first so it is plotted first (background)
efilt<-reshape(subset(dfilt, NSamples == 45), direction="wide", idvar=c("Contig", "NSamples"), timevar="Group")

genericimage <- function(data) {
  theme_set(theme_bw())
  
  paultol    =c("#DDDDDD", "#88CCEE",  "#DDCC77", "#CC6677", "#117733", "#332288", "#44AA99", "#999933", "#AA4499", "#882255", "#777777")
  g <- ggplot() + facet_wrap(~label) + scale_colour_manual(values=paultol, name="Taxonomic\nClassification", limits=levels(data$taxon) )
  for (t in levels(data$taxon)) {
    if(nrow(data[data$taxon==t,]) > 0) {
      g <- g + geom_point(data=data[data$taxon==t,],aes(gc, cov, colour=taxon), size=2, alpha=I(1/3))
    }
  }
  axis_breaks = c(1,2,5,10,20,50,100,200,500,1000,10000);
  g + xlim(0,1) +
    ylim(0,10000)  + scale_y_log10(breaks=axis_breaks) + scale_x_continuous(labels=percent) +
    labs(x="GC content", y="Read coverage")
}

stecs1<-subset(efilt, NSamples == 45 & PresentIn.stecs > 1)
juststecs1<-with(stecs1, data.frame(cov=cov.stecs, taxon=taxon.stecs, gc=gc.stecs,label="A" ))
#b<-genericimage(juststecs1)

stecs20<-subset(efilt, NSamples == 45 & PresentIn.stecs >= 20)

write.table(stecs20, file="stec20contigs.txt", sep="\t", quote=FALSE)

juststecs20<-with(stecs20, data.frame(cov=cov.stecs, taxon=taxon.stecs, gc=gc.stecs,label="B" ))

justecoli<-subset(stecs20, taxon.stecs == 'Enterobacteriales')
justecoliplot<-with(justecoli, data.frame(cov=cov.stecs, taxon=taxon.stecs, gc=gc.stecs,label="X"))
genericimage(justecoliplot)
c<-Mclust(with(justecoliplot, data.frame(log(cov), gc)))
plot(c)


#juststecsbloody<-with(stecsbloody,  data.frame(cov=cov.stecs, taxon=taxon.stecs, gc=gc.stecs,label="Bloody" ))
#stecsbloody<-subset(efilt, NSamples == 45 & PresentIn.stecs >= 20 & PresentIn.stecsbloody>=9)
b<-genericimage(juststecs)

subtract<-subset(efilt, PresentIn.stecs >= 20 & PresentIn.metahitsamps < 1 & NSamples == 45)

#subtract<-subset(efilt,
#                 PresentIn.stecs >= 20
#                 & PresentIn.metahitsamps < 1
#                 & PresentIn.stecsmiseq >= 1
#                 & NSamples == 45)

subtracted<-with(subtract, data.frame(cov=cov.stecs, taxon=taxon.stecs, gc=gc.stecs, label="C" ))

newdf<-rbind(juststecs1, juststecs20, subtracted)
genericimage(newdf)

# fun animation
library(animation)





saveMovie({
i<-1
while(i<=43) {
 stecs1<-subset(efilt, NSamples == 45 & PresentIn.stecs >= i)
juststecs1<-with(stecs1, data.frame(cov=cov.stecs, taxon=taxon.stecs, gc=gc.stecs,label=paste("In >=", i, " samples")))
 print(genericimage(juststecs1))
 i<-i+1
}
}, interval=0.5, movie.name="r-animation2.gif", width=800, height=600)

library(animation)

pdf("pathogenenrichment.pdf", width=10, height=6)
ggsave("pathogenenrichment.pdf", genericimage(newdf), width=10, height=6, units="cm")
dev.off()
grid.arrange(a, b, c, ncol=3, nrow=1)
dev.off()

simpleimage <- function(data) {
  theme_set(theme_bw())
  
  paultol    =c("#DDDDDD", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#44AA99", "#999933", "#AA4499", "#882255", "#777777")
  g <- ggplot() 
    #scale_colour_manual(values=paultol, name="Taxonomic\nClassification", limits=levels(data$taxon.stecs) )
  for (t in levels(data$taxon.stecs)) {
     if(nrow(data[data$taxon.stecs==t,]) > 0) {
       g <- g + geom_point(data=data[data$taxon.stecs==t,],aes(gc.stecs, cov.stecs, colour=taxon.stecs), size=2, alpha=I(1/3))
    }
  }
  axis_breaks = c(1,2,5,10,20,50,100,200,500,1000,10000);
  g + xlim(0,1) +
    ylim(0,max(dfilt$cov, na.rm=TRUE))  + scale_y_log10(breaks=axis_breaks) + 
    labs(x="GC content", y="Read coverage") + opts(legend.position = "none")
}

g1<-simpleimage(ffilt)
grid.arrange(g1, g1, g1, ncol=3, nrow=1)
dev.off()

pdf("pathogenenrichment.pdf", width=8, height=6)
simpleimage(ffilt)
dev.off()

dev.off()
png(paste("out2.png",sep="."), 1000*3, 1000*3, units="px",res=100)
blob.plot(subset(dfilt, Group == "metahitsamps"))

sum(subset(ffilt, taxon.stecs == 'Enterobacteriales')$ContigLen.stecs)
sum(subset(ffilt, taxon.stecs == 'Not annotated')$ContigLen.stecs)
sum(subset(ffilt, taxon.stecs != 'Not annotated' & taxon.stecs != 'Enterobacteriales')$ContigLen.stecs)
median(subset(ffilt, taxon.stecs == 'Enterobacteriales')$gc.stecs)
median(subset(ffilt, taxon.stecs == 'Not annotated')$gc.stecs)


blob.plot(ffilt)
dev.off()
subset(dfilt, V58 > 41)$gc



#model-based clustering

library(mclust)
formodel1<-cbind(subtracted, subtract$Contig)
formodel2<-with(subtracted, data.frame(log(cov), gc))
mdl<-Mclust(formodel2)

mdlgroups<-cbind(formodel1, class=mdl$classification)
mdlgroups<-cbind(mdlgroups, uncert=mdl$uncertainty)
write.table(mdlgroups, file="contiggroups.txt", sep="\t", quote=FALSE)

