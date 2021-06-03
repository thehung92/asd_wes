#!/usr/bin/env Rscript
#### library ####
library("Gviz")
library("tidyverse")
library("data.table")
library("lemon")
library("biomaRt")
library("gridtext")
library("ggpubr")
source("src/function_plot-gene-region.R")
# debug for ggplot v3.3.3
pdf(file="Rplots.pdf")
dev.off()
#### import data ####
tdt <- fread(file="./output/tdt/asd.290_ibd-clean.tdt", header=TRUE)
as_tibble(tdt) %>%
  arrange(P) %>%
  mutate(CHR=sub_chr_ensembl(CHR))%>%
  mutate(LOG=-log10(P)) %>%
  mutate_at("CHR", ~gsub("^","chr",.)) -> tdt
recomb.rate <- fread(file="./data/genetic_map_GRCh38_merged.tab.gz", header=TRUE)
# hg38; ensembl database.
gene.ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        dataset="hsapiens_gene_ensembl")
#### create annotation ####
# y axis on the right side
col.axis="blue"
p2 <- ggplot(head(recomb.rate), aes(pos, recomb_rate)) +
  geom_blank() +
  labs(y="Recombination Rate\n(cM/Mb)") +
  scale_y_continuous(position="right", limits=c(0,100)) +
  # capped the axis
  coord_capped_cart(right="both") +
  theme(axis.line.y=element_line(color=col.axis),
        axis.ticks.y=element_line(color=col.axis),
        axis.text.y=element_text(color=col.axis),
        text=element_text(family="Helvetica", face="bold",color=col.axis),
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank())
# create bottom legend
set.seed(12345)
df <- tibble(x=1:20,
             y=runif(20,0,8),
             z=runif(20,0,1))
df$range <- cut(df$z, breaks=seq(0,1,0.2))
# annotate lead dot
df$range <- as.character(df$range)
df$range[1] <- "Lead.variant"
# create breaks factor
break.factor <- as.factor(df$range) %>% levels()
# create color palette
col.palette <- c("darkblue","lightblue","green","orange","red","darkviolet")
# create label
break.label <- c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1","Lead variant")
legend.title <- parse(text="r^2~`:`")
# function for small legend
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 6, spaceLegend = 0.2) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize),
                                nrow=1, ncol=6)#,
           #color = guide_legend(override.aes = list(size = pointSize),
            #                    nrow=1, ncol=6)
           ) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
# plot scatter with horizontal legend
p <- ggplot(df, aes(x,y,color=range,shape=range, group=range)) +
  geom_point() +
  scale_color_manual(legend.title,
                     breaks=break.factor,
                     values=col.palette,
                     #name=parse(text="r^2~`:`"),
                     labels=break.label
  ) +
  scale_shape_manual(legend.title,
                     breaks=break.factor,
                     values=c(rep(20,length(break.factor)-1),18),
                     #name=parse(text="r^2~`:`"),
                     labels=break.label
  ) +
  theme_light() +
  theme(text=element_text(family="Helvetica", color="grey50"),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.background=element_rect(size=0.2,linetype="solid",
                                       fill="transparent")
  ) +
  guides(color=guide_legend(nrow=1,ncol=6))
p <- addSmallLegend(p,2.5,8,0.5)
# extract legend
leg <- get_legend(p) %>% as_ggplot()
#### create dashed line ####
text <- "Threshold:  <span style='color:red'>**- - -**</span>  Bonferroni correction"
#
g <- richtext_grob(text, hjust=0, vjust=0.5,
                   gp=gpar(fontfamily="Helvetica", fontsize=8, col="#808080"))
#### loop plot ####

for (i in 1:4) {
  # i=1
  print(paste("Plotting variant",i,"using gviz"))
  #
  tdt.sel <- tdt %>% dplyr::slice(i)
  sel.chr <- tdt.sel$CHR
  sel.pos <- tdt.sel$BP
  range <- 5e5
  gen <- "hg38"
  # create chromosome ideogram and genome axis track
  cat("creating ideogram & genomic axis...\n")
  gtrack <- GenomeAxisTrack(labelPos="below")
  itrack <- IdeogramTrack(genome=gen, chromosome= sel.chr)
  # create adjacent variants Pvalue track
  cat("creating pvalue data track...\n")
  tdt %>%
    filter(CHR==sel.chr, between(BP, sel.pos-range, sel.pos+range)) %>%
    dplyr::select(c(1,3,11)) -> tdt.sel.region
  # create annotation tracks of variants
  cat("creating annotation track...\n")
  tdt %>%
    filter(CHR==sel.chr, between(BP, sel.pos-range, sel.pos+range)) %>%
    dplyr::select(c(1,2,3)) -> sel.variant
  sel.variant[-1,2] <- " "
  sel.variant %>%
    dplyr::rename(chromosome=1, id=2, start=3) %>%
    mutate(end=start,group=c("sel",rep("adj",nrow(sel.variant)-1))) -> sel.variant
  # make transparent blue color
  col0 <- rgb(0, 0, 255, max = 255, alpha = 100, names = "blue40")
  # atrack with feature colored
  atrack <- AnnotationTrack(name="Variants", genome=gen, chromosome=sel.chr,
                            start=sel.variant$start, end=sel.variant$end,
                            id=sel.variant$id,featureAnnotation="id",
                            fontcolor.feature="black",
                            shape="box",stacking="dense", #below is new code
                            feature=rep(c("selected","adjacent"),c(1,nrow(sel.variant)-1)),
                            col="transparent", selected="red", adjacent=col0)
  displayPars(atrack) <- list(background.title="transparent",
                              fontcolor.title="black",
                              rotation.title=0,showTitle=TRUE,
                              cex.title=0.8)
  # create gene region track
  cat("creating gene region track...\n")
  # query needed biomart field
  out.bm.genes.region <- getBM(
    attributes = c('chromosome_name','exon_chrom_start','exon_chrom_end','strand',
                   'gene_biotype',
                   'ensembl_gene_id','ensembl_exon_id','ensembl_transcript_id',
                   'external_gene_name'), 
    filters = c('chromosome_name','start','end'), 
    values = list(gsub("chr","",sel.chr),sel.pos - range, sel.pos + range), 
    mart = gene.ensembl)
  # reformat dataframe for plotting gene region track
  out.bm.genes.region %>%
    dplyr::rename(chromosome=1,
                  start=2,
                  end=3,
                  strand=4,
                  feature=5,
                  gene=6,
                  exon=7,
                  transcript=8,
                  symbol=9) %>%
    mutate(strand=sub_strand_gviz(strand)) %>%
    filter(feature=="protein_coding") %>%
    mutate(symbol=paste0(symbol,strand)) %>%
    mutate_at("symbol",~gsub("\\+$","\u2192",.)) %>%
    mutate_at("symbol",~gsub("\\-$","\u2190",.)) -> genes.region
  grtrack <- GeneRegionTrack(genes.region,name="Known\nGenes\n(Ensembl)",
                             genome=gen,chromosome=sel.chr,
                             transcriptAnnotation="symbol",
                             collapseTranscripts="longest")
  displayPars(grtrack) <- list(rotation.title=0,
                               background.title="transparent",
                               fontcolor.title="black",
                               fontfamily.group="Helvetica",
                               just.group="above")
  # create recombination rate track
  cat("creating recombination rate track...\n")
  recomb.rate %>%
    filter(chrom==sel.chr, between(pos, sel.pos-range, sel.pos+range)) -> recomb.rate.sel
  recomb.rate.sel %>%
    arrange(pos) %>%
    mutate(pos_end=c(pos[-1],pos[nrow(recomb.rate.sel)])) %>%
    dplyr::select(1,2,5,3) -> recomb.df
  grange.recomb <- makeGRangesFromDataFrame(recomb.df, keep.extra.columns=TRUE,
                                            ignore.strand=TRUE,
                                            seqnames.field="chrom",
                                            start.field="pos",
                                            end.field="pos_end")
  rrtrack <- DataTrack(grange.recomb, name="Recomb.\nrate\n(cM/Mb)", genome=gen,
                       type="l")
  # create transparent blue
  col1 <- rgb(0, 0, 255, max = 255, alpha=round(0.7*255), names = "blue70")
  displayPars(rrtrack) <- list(ylim=c(0,100),
                               col=col1)
  # plot pvalue data track with grouping based on r2 correlation
  cat("creating pvalue data track with r2 grouping...\n")
  file=paste0("temp/ld-r2/matrix.r2.unrelated.variant",i,".txt")
  df.ldr2 <- read.csv(file, sep="")
  # select r2 value with lead variant
  df.ldr2 %>%
    dplyr::select(tdt.sel$SNP) %>%
    mutate(SNP=rownames(.)) %>%
    dplyr::rename(r2=1) -> df.ldr2
  df.ldr2$class=NA
  df.ldr2$class[1]="Lead.variant"
  for (j in 2:nrow(df.ldr2)) {
    x <- as.numeric(df.ldr2$r2[j])
    if (x<=0.2) {df.ldr2$class[j]="0-0.2"}
    if (x>0.2 & x<=0.4) {df.ldr2$class[j]="0.2-0.4"}
    if (x>0.4 & x<=0.6){df.ldr2$class[j]="0.4-0.6"}
    if (x>0.6 & x<=0.8){df.ldr2$class[j]="0.6-0.8"}
    if (x>0.8 & x<=1){df.ldr2$class[j]="0.8-1"}
  }
  # create data frame for dtrack
  df.ldr2 %>% 
    merge(x=., y=tdt[,c(1:3,11)],by="SNP") %>%
    pivot_wider(names_from = class, values_from=LOG) -> df.dtrack
  df.dtrack %>%
    dplyr::select(-1,-2) %>%
    rename(chromosome=1,start=2) %>%
    mutate(end=start) %>%
    relocate(end, .after=start) -> df.dtrack
  # create dtrack with grouping
  dtrack.group <- DataTrack(df.dtrack,genome="hg38",name="-log10(P-value)",
                            groups=colnames(df.dtrack)[-1:-3],
                            baseline=-log10(0.05/nrow(tdt)),col.baseline="red",
                            lty.baseline="dashed",type="p",ylim=c(0,9))
  # edit color palette
  col.palette0 <- c("darkblue","lightblue","green","orange","red","darkviolet")
  names(col.palette0) <- c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1","Lead.variant") # expected class for color palette
  r2.class <- levels(as.factor(df.ldr2$class)) # available class for r2
  col.palette1 <- col.palette0[names(col.palette0) %in% r2.class] # filter color palette
  # barplot(rep(2,6), col=col.palette)
  displayPars(dtrack.group) <- list(background.title = "transparent",
                                    col.axis="black",
                                    col.title="black",
                                    box.legend=FALSE, legend=FALSE,
                                    col=col.palette1,
                                    pch=c(rep(20,length(r2.class)-1),18), cex=1, # shape and size of dot
                                    cex.axis=0.8, cex.title=0.8)
  #
  cat("creating overlay and highlight lead variant...\n")
  # overlay rrtrack on top of dtrack.group
  ovltrack <- OverlayTrack(trackList=list(dtrack.group, rrtrack))
  displayPars(ovltrack) <- list(background.title="transparent")
  # hightlight track
  ht <- HighlightTrack(trackList=list(atrack, grtrack,ovltrack,gtrack),
                       start=sel.pos-2000, width=4000,
                       chromosome=sel.chr,
                       col="#FFE3E6")
  # start plotting with grid viewport on png device
  cat("plotting with grid viewport on png device...\n")
  output=paste0("./output/plot-table/variant_",i,".",sel.chr,"_",sel.pos,".png")
  png(file=output, bg="white",
      units="mm", width=175, height=100, res=300)
  grid.newpage()
  vp_1 <- viewport(x=0,y=0,width=0.9,height=1,just=c("left","bottom"))
  pushViewport(vp_1)
  Gviz::plotTracks(list(itrack, ht), sizes=c(1,1,3,5,2),
                   from=sel.pos-range, to=sel.pos+range,
                   title.width=0.8,
                   add=TRUE)
  popViewport()
  vp_2 <- viewport(x=0,y=1.85/12,width=0.99,height=5.15/12,just=c("left","bottom")) # manual edit to fit
  pushViewport(vp_2)
  grid.draw(ggplotGrob(p2))
  popViewport()
  vp_3 <- viewport(x=0.35,y=1.8/12, width=0.8, height=2/12, just="centre")
  pushViewport(vp_3)
  grid.draw(ggplotGrob(leg))
  popViewport()
  vp_4 <- viewport(x=0.65,y=1.8/12, width=0.8, height=2/12, just="centre")
  pushViewport(vp_4)
  grid.draw(g)
  popViewport()
  dev.off()
  # report output
  print(output)
}
# rename figure
DIR="output/plot-table/"
FILES <- list.files(DIR, pattern="variant.*png") %>% paste0(DIR, .)
FILES2 <- paste0(DIR, c("figure3", "figureS2", "figureS3", "figureS4"), ".png")
CMD=paste("mv", FILES, FILES2)
for (i in 1:length(CMD)) {
  system(CMD[i])
}
#### reset before exit ####
# rm(list=ls())
