library(ape)
library(phytools)
library(viridis)
setwd("~/Dropbox/Data/p50_global/phylo")

# read in traits and tree
trait <- read.csv("trait_data_sPlot.csv", header=TRUE)
tree <- read.tree("p50_rd_phylo.tre")

# to ensure correct order of tips and traits
rownames(trait) <- trait$species %>% str_replace_all(" ", "_")       ###add row names
p.dist.mat <- cophenetic(tree)         ###to ensure correct order
trait <- trait[row.names(p.dist.mat),] ###to ensure correct order
tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")
rownames(trait) <- tree$tip.label

# radial phylo with two traits
p50 <- setNames(trait[,2],trait$species)
rd <- setNames(log(trait[,3]+1,base=10),trait$species)

tiff("Fig1.tiff",width=3000,height=3000,res=300)
obj <- contMap(tree, p50, plot=F, ftype="i", fsize = 2) %>% 
  setMap(viridis::viridis_pal(option = "inferno")(10))
plotTree.wBars(obj$tree, rd, method = "plotSimmap", colors = obj$cols,
               fsize=0.5,tip.labels=TRUE,lwd=3,
               type = "fan", scale = 30, col=1)
text(110,-110,"Gymnosperms",srt=45)
text(-50,-150,"Basal
eudicots", srt=-20)
text(-140,80,"Superrosids", srt=60)
text(160,80,"Superasterids", srt=-50)
add.color.bar(180, obj$cols, title = expression("P"[50]*" (MPa)"), lims = obj$lims, digits=1, prompt = FALSE,
              x = -80, y = 35)
text(-150,-412,"Rooting depth (m)", srt=-20)
text(-53,-520,"60")
text(-285,-430,"30")
text(-203,-455,"8")
text(-108,-465,"1")
dev.off()
