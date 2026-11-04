# Fig. 2b: unique adjusted R2 for four nitrogen forms in ARG/community profiles.
# Run from repository root: Rscript 02_environmental_drivers/fig2b_nitrogen_unique_r2.R

required <- c("openxlsx", "vegan", "ggplot2", "svglite", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(invisible(lapply(required, library, character.only = TRUE)))
script <- commandArgs(trailingOnly = FALSE); script <- sub("^--file=", "", script[grep("^--file=", script)[1]])
ROOT <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE); DATA <- file.path(ROOT, "data", "Supplementary_Tables.xlsx")
OUT <- file.path(ROOT, "outputs"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(DATA)) stop("Missing input workbook: ", DATA)
PERMUTATIONS <- 9999L
site_id <- function(x) sub("_?[0-9]+$", "", sub("_memb$", "", as.character(x)))
read_table <- function(s) openxlsx::read.xlsx(DATA, sheet = s, startRow = 2, check.names = FALSE, na.strings = c("", "NA", "NaN"))
site_mean <- function(x, ids) { x <- as.data.frame(lapply(x, function(v) as.numeric(as.character(v)))); z <- aggregate(x, list(site = ids), mean, na.rm = FALSE); rownames(z) <- z$site; z$site <- NULL; z }
ilr_land <- function(env) { land <- c("Tree.cover", "Shrubland", "Grassland", "Cropland", "Built.up", "Bare.vegetation", "Herbaceous.wetland"); z <- env[, land, drop=FALSE] / rowSums(env[, land, drop=FALSE]); positive <- z[z > 0]; z[z == 0] <- min(positive)/2; z <- z/rowSums(z); out <- matrix(NA_real_, nrow(z), 6, dimnames=list(rownames(z), paste0("land_ILR",1:6))); for (j in 1:6) out[,j] <- sqrt(j/(j+1)) * (rowMeans(log(z[,land[seq_len(j)],drop=FALSE])) - log(z[,land[j+1]])); out }
zscore <- function(x) as.data.frame(lapply(x, function(v) as.numeric(scale(v))))
t4 <- read_table("Table 4"); t5 <- read_table("Table 5"); t7 <- read_table("Table 7"); samples <- as.character(t7$Sample); sites <- site_id(samples)
if (length(unique(sites)) != 89L || any(table(sites) != 3L)) stop("Expected 89 sites with 3 replicates each")
env <- site_mean(t7[,setdiff(names(t7),"Sample"),drop=FALSE], sites); if (any(!complete.cases(env))) stop("Incomplete Table 7 site means")
a <- as.matrix(t4[,setdiff(names(t4),"Sample"),drop=FALSE]); storage.mode(a) <- "numeric"; an <- sub("__.*$", "", colnames(a)); a <- sapply(unique(an), function(n) rowSums(a[,an==n,drop=FALSE])); colnames(a) <- unique(an); arg <- site_mean(as.data.frame(a),sites)[rownames(env),,drop=FALSE]; arg[is.na(arg)] <- 0; arg <- arg[,colSums(arg)>0,drop=FALSE]
g <- as.character(t5[[7]]); cmat <- as.matrix(t5[,samples,drop=FALSE]); storage.mode(cmat) <- "numeric"; cmat <- sapply(unique(g), function(n) colSums(cmat[g==n,,drop=FALSE])); colnames(cmat) <- unique(g); comm <- site_mean(as.data.frame(cmat),sites)[rownames(env),,drop=FALSE]; comm[is.na(comm)] <- 0; comm <- comm[,colSums(comm)>0,drop=FALSE]
response <- function(x) { h <- vegan::decostand(x,"hellinger"); p <- prcomp(h,center=FALSE,scale.=FALSE); p$x[,which(p$sdev > max(p$sdev)*1e-10),drop=FALSE] }
arg_response <- response(arg); comm_response <- response(comm); land <- ilr_land(env); non_nitrogen <- c("pH","Temperature","EC","Turbidity","COD","TOC","TP"); forms <- c("NO3_N","TN","NO2_N","NH3_N")
make_x <- function(focal) { x <- cbind(zscore(env[,c(focal,non_nitrogen,"NTL"),drop=FALSE]), land); x[,c(focal,non_nitrogen,colnames(land),"NTL"),drop=FALSE] }
run_focals <- function(y, label) {
  do.call(rbind, lapply(forms, function(focal) {
    x <- make_x(focal)
    full <- vegan::rda(y ~ ., x)
    reduced <- vegan::rda(y ~ ., x[, setdiff(names(x), focal), drop = FALSE])
    test <- vegan::anova.cca(full, by = "margin", permutations = PERMUTATIONS)
    data.frame(response = label, focal = focal,
      unique_R2 = vegan::RsquareAdj(full)$adj.r.squared - vegan::RsquareAdj(reduced)$adj.r.squared,
      adjusted_R2 = vegan::RsquareAdj(full)$adj.r.squared, P = test[focal, "Pr(>F)"])
  }))
}
z <- rbind(run_focals(arg_response,"ARG profile"), run_focals(comm_response,"Microbial community")); z$q <- ave(z$P,z$response,FUN=function(x)p.adjust(x,"BH")); z$significant <- z$q < .05; z$focal <- factor(z$focal,forms,labels=c("NO3-N","TN","NO2-N","NH3-N")); z$q_label <- ifelse(z$q < .001,"q < 0.001",sprintf("q = %.3f",z$q)); z$y_value <- 100*pmax(z$unique_R2,0)
p <- ggplot(z,aes(focal,y_value,colour=focal)) + geom_segment(aes(xend=focal,y=0,yend=y_value),linewidth=.45,show.legend=FALSE) + geom_point(size=2.3,shape=16) + geom_text(aes(label=q_label),vjust=-.8,colour="#20252A",size=1.8) + facet_grid(response~.,scales="free_y",switch="y") + scale_colour_manual(values=c("NO3-N"="#FDAC2C","TN"="#4052A5","NO2-N"="#98145A","NH3-N"="#4B98B7"),guide="none") + scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.18)),labels=function(x)sprintf("%g",x)) + labs(x=NULL,y="Unique adjusted R² (%)") + theme_classic(base_family="Arial",base_size=7) + theme(panel.border=element_rect(fill=NA,colour="black",linewidth=.3),strip.background=element_blank(),strip.placement="outside",strip.text.y.left=element_text(face="bold",angle=90),panel.spacing.y=grid::unit(5,"mm"),axis.line=element_blank())
for (ext in c("pdf","svg","png")) { f <- file.path(OUT,paste0("02b_nitrogen_unique_r2.",ext)); if(ext=="pdf")ggsave(f,p,width=58,height=132,units="mm",device=cairo_pdf); if(ext=="svg")ggsave(f,p,width=58,height=132,units="mm",device=svglite::svglite); if(ext=="png")ggsave(f,p,width=58,height=132,units="mm",dpi=600,device=ragg::agg_png) }
cat("Fig. 2b images written to ",OUT,"\n",sep="")
