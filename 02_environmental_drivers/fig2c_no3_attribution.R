# Fig. 2c: NO3-N variation partitioning and ARG Type/Class association.
# Run from repository root: Rscript 02_environmental_drivers/fig2c_no3_attribution.R

required <- c("openxlsx", "vegan", "ggplot2", "ggrepel", "svglite", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]; if(length(missing)) stop("Install required R packages: ",paste(missing,collapse=", "))
suppressPackageStartupMessages(invisible(lapply(required, library, character.only = TRUE)))
script <- commandArgs(trailingOnly=FALSE); script <- sub("^--file=","",script[grep("^--file=",script)[1]]); ROOT <- normalizePath(file.path(dirname(script),".."),mustWork=TRUE); DATA <- file.path(ROOT,"data","Supplementary_Tables.xlsx"); OUT <- file.path(ROOT,"outputs"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE); if(!file.exists(DATA)) stop("Missing input workbook: ",DATA); PERMUTATIONS <- 9999L
site_id <- function(x) sub("_?[0-9]+$","",sub("_memb$","",as.character(x))); read_table <- function(s) openxlsx::read.xlsx(DATA,sheet=s,startRow=2,check.names=FALSE,na.strings=c("","NA","NaN")); site_mean <- function(x,ids){x<-as.data.frame(lapply(x,function(v)as.numeric(as.character(v))));z<-aggregate(x,list(site=ids),mean,na.rm=FALSE);rownames(z)<-z$site;z$site<-NULL;z}; zscore <- function(x)as.data.frame(lapply(x,function(v)as.numeric(scale(v))))
ilr_land <- function(env){land<-c("Tree.cover","Shrubland","Grassland","Cropland","Built.up","Bare.vegetation","Herbaceous.wetland");z<-env[,land,drop=FALSE]/rowSums(env[,land,drop=FALSE]);positive<-z[z>0];z[z==0]<-min(positive)/2;z<-z/rowSums(z);out<-matrix(NA_real_,nrow(z),6,dimnames=list(rownames(z),paste0("land_ILR",1:6)));for(j in 1:6)out[,j]<-sqrt(j/(j+1))*(rowMeans(log(z[,land[seq_len(j)],drop=FALSE]))-log(z[,land[j+1]]));out}
t4<-read_table("Table 4");t5<-read_table("Table 5");t7<-read_table("Table 7");samples<-as.character(t7$Sample);sites<-site_id(samples);if(length(unique(sites))!=89L||any(table(sites)!=3L))stop("Expected 89 sites with 3 replicates each");env<-site_mean(t7[,setdiff(names(t7),"Sample"),drop=FALSE],sites);if(any(!complete.cases(env)))stop("Incomplete Table 7 site means")
a<-as.matrix(t4[,setdiff(names(t4),"Sample"),drop=FALSE]);storage.mode(a)<-"numeric";an<-sub("__.*$","",colnames(a));a<-sapply(unique(an),function(n)rowSums(a[,an==n,drop=FALSE]));colnames(a)<-unique(an);arg<-site_mean(as.data.frame(a),sites)[rownames(env),,drop=FALSE];arg[is.na(arg)]<-0;arg<-arg[,colSums(arg)>0,drop=FALSE]
g<-as.character(t5[[7]]);cm<-as.matrix(t5[,samples,drop=FALSE]);storage.mode(cm)<-"numeric";cm<-sapply(unique(g),function(n)colSums(cm[g==n,,drop=FALSE]));colnames(cm)<-unique(g);comm<-site_mean(as.data.frame(cm),sites)[rownames(env),,drop=FALSE];comm[is.na(comm)]<-0;comm<-comm[,colSums(comm)>0,drop=FALSE]
response<-function(x){h<-vegan::decostand(x,"hellinger");p<-prcomp(h,center=FALSE,scale.=FALSE);p$x[,which(p$sdev>max(p$sdev)*1e-10),drop=FALSE]}; arg_response<-response(arg);comm_response<-response(comm);land<-ilr_land(env);non_nitrogen<-c("pH","Temperature","EC","Turbidity","COD","TOC","TP")
groups<-list(NO3_N=zscore(env[,"NO3_N",drop=FALSE]),Other_WQ=zscore(env[,non_nitrogen,drop=FALSE]),Land_Urban=cbind(zscore(env[,"NTL",drop=FALSE]),land));
partition_one<-function(y,label){vp<-vegan::varpart(y,groups$NO3_N,groups$Other_WQ,groups$Land_Urban);ind<-vp$part$indfract[,"Adj.R.square"];data.frame(response=label,component=c("NO3-N unique","Other water quality","Land cover / urbanization","Shared / overlap","Unexplained"),fraction=c(ind[1],ind[2],ind[3],sum(ind[4:7],na.rm=TRUE),ind[8]))}
partition<-rbind(partition_one(arg_response,"ARG profile"),partition_one(comm_response,"Microbial community"))
assoc<-do.call(rbind,lapply(colnames(arg),function(v){ok<-is.finite(arg[,v])&is.finite(env$NO3_N);ct<-suppressWarnings(cor.test(env$NO3_N[ok],arg[ok,v],method="spearman",exact=FALSE));data.frame(arg_type=v,rho=unname(ct$estimate),P=ct$p.value,n_sites=sum(ok))}));assoc$q<-p.adjust(assoc$P,"BH");assoc$symbol<-ifelse(assoc$q<.001,"***",ifelse(assoc$q<.01,"**",ifelse(assoc$q<.05,"*","")))
partition$component <- factor(partition$component, levels = c(
  "NO3-N unique", "Other water quality", "Land cover / urbanization",
  "Shared / overlap", "Unexplained"))
partition$fraction <- pmax(partition$fraction, 0)
partition$ring <- ifelse(partition$response == "Microbial community", 1, 2)
partition$start <- ave(partition$fraction, partition$response,
  FUN = function(x) c(0, head(cumsum(x), -1)))
partition$end <- partition$start + partition$fraction

sector <- function(row) {
  angle <- seq(2 * pi * row$start + pi / 2,
    2 * pi * row$end + pi / 2, length.out = 32)
  inner <- if (row$ring == 1) .62 else 1.16
  outer <- if (row$ring == 1) .86 else 1.43
  data.frame(
    x = c(inner * cos(angle), rev(outer * cos(angle))),
    y = c(inner * sin(angle), rev(outer * sin(angle))),
    component = row$component,
    group = paste(row$response, row$component))
}
poly <- do.call(rbind, lapply(seq_len(nrow(partition)), function(i) sector(partition[i, ])))

assoc <- assoc[order(assoc$rho), , drop = FALSE]
assoc$angle <- 2 * pi * (seq_len(nrow(assoc)) - 1) / nrow(assoc) + pi / 2
assoc$x <- 2.20 * cos(assoc$angle)
assoc$y <- 2.20 * sin(assoc$angle)
assoc$label_x <- 2.62 * cos(assoc$angle)
assoc$label_y <- 2.62 * sin(assoc$angle)
assoc$direction <- ifelse(assoc$rho >= 0, "Positive", "Negative")

arg_colors <- c(
  "macrolide.lincosamide.streptogramin"="#E64B35", "aminoglycoside"="#4DBBD5",
  "glycopeptide"="#00A087", "multidrug"="#3C5488", "polymyxin"="#F39B7F",
  "mupirocin"="#8491B4", "pleuromutilin_tiamulin"="#7E6148", "beta_lactam"="#91D1C2",
  "quinolone"="#B09C85", "sulfonamide"="#DC0000", "tetracycline"="#C77CFF",
  "chloramphenicol"="#F4A261", "florfenicol"="#2A9D8F", "trimethoprim"="#0072B2",
  "antibacterial_fatty_acid"="#D4A72C", "rifamycin"="#6A3D9A", "vancomycin"="#1B9E77",
  "fosfomycin"="#D95F02", "tunicamycin"="#E7298A", "bicyclomycin"="#66A61E",
  "defensin"="#A6761D", "streptothricin"="#E6AB02", "bacitracin"="#8C510A",
  "unclassified"="#7D8790", "other_peptide_antibiotics"="#666666", "factumycin"="#8C6BB1",
  "fusidic_acid"="#1F78B4", "tetracenomycin_C"="#33A02C", "novobiocin"="#FF7F00",
  "puromycin"="#984EA3", "edeine"="#B15928", "bleomycin"="#A6CEE3")
missing_colors <- setdiff(assoc$arg_type, names(arg_colors))
if (length(missing_colors)) stop("Missing ARG class colors: ", paste(missing_colors, collapse = ", "))
assoc$color <- unname(arg_colors[assoc$arg_type])

partition_colors <- c(
  "NO3-N unique"="#1F78B4", "Other water quality"="#86B7DB",
  "Land cover / urbanization"="#AEB8C2", "Shared / overlap"="#D7DDE2",
  "Unexplained"="#F1F3F5")
ring_labels <- data.frame(
  x = c(-1.02, -1.62), y = c(.25, .02),
  label = c("Microbial community", "ARG profile"))

p <- ggplot() +
  geom_polygon(data = poly, aes(x, y, group = group, fill = component),
    colour = "white", linewidth = .35) +
  geom_point(data = assoc, aes(x, y, colour = color, size = abs(rho), shape = direction),
    fill = "white", stroke = .7) +
  ggrepel::geom_text_repel(data = assoc,
    aes(x, y, label = gsub("_", " ", arg_type)),
    nudge_x = assoc$label_x - assoc$x, nudge_y = assoc$label_y - assoc$y,
    colour = "#20252A", size = 1.65, family = "Arial", lineheight = .82,
    box.padding = .08, point.padding = .08, segment.colour = NA,
    force = .7, force_pull = 1.5, max.overlaps = Inf, seed = 20260908) +
  geom_text(data = assoc[assoc$symbol != "", ], aes(x, y, label = symbol),
    colour = "black", size = 1.15, family = "Arial", fontface = "bold", vjust = .52) +
  geom_text(data = ring_labels, aes(x, y, label = label),
    colour = "#20252A", size = 1.55, family = "Arial", fontface = "bold", hjust = 0) +
  scale_colour_identity() +
  scale_fill_manual(name = "Variation partitioning", values = partition_colors,
    guide = guide_legend(order = 1, nrow = 2, byrow = TRUE)) +
  scale_shape_manual(name = "rho direction", values = c(Positive = 16, Negative = 21),
    guide = guide_legend(order = 2, override.aes = list(colour = "#5D6870", size = 2.7))) +
  scale_size_continuous(name = "|Spearman rho|", range = c(1.8, 4.6),
    breaks = c(.1, .3, .5), guide = guide_legend(order = 3,
      override.aes = list(colour = "#5D6870"))) +
  coord_fixed(xlim = c(-3.55, 3.55), ylim = c(-3.35, 3.35), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(4, 18, 4, 18, unit = "mm"),
    legend.position = "bottom", legend.box = "vertical", legend.box.just = "center",
    legend.title = element_text(size = 5.4), legend.text = element_text(size = 5.2))

for (ext in c("pdf", "svg", "png")) {
  file <- file.path(OUT, paste0("02c_no3_attribution.", ext))
  if (ext == "pdf") ggsave(file, p, width = 145, height = 135, units = "mm", device = cairo_pdf)
  if (ext == "svg") ggsave(file, p, width = 145, height = 135, units = "mm", device = svglite::svglite)
  if (ext == "png") ggsave(file, p, width = 145, height = 135, units = "mm", dpi = 600, device = ragg::agg_png)
}
cat("Fig. 2c images written to ", OUT, "\n", sep = "")
