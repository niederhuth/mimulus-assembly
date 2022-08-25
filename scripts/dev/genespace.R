#Set variables
gids <- sids <- vids <- c("L1","S1")
ploidy <- c(1,1)
outgroup <- NULL
runwd <- getwd()
nCores <- 4
path2mcscanx <- "MCScanX"

#Check for BiocManager & install
if (!require("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
#Check for devtools & install
if (!requireNamespace("devtools", quietly = TRUE))
	install.packages("devtools")
#Check for Biostrings & install
if (!require("Biostrings", quietly = FALSE))
	BiocManager::install("Biostrings")
#Check for rtracklayer & install
if (!require("rtracklayer", quietly = FALSE))
	BiocManager::install("rtracklayer")
#Check for BiocManager & install 
if (!require("GENESPACE", quietly = FALSE))
	devtools::install_github("jtlovell/GENESPACE@v0.9.3")

#Load libraries
library(GENESPACE)

#Initialize genespace and set paramaters
gpar <- init_genespace(
	genomeIDs = gids,
	speciesIDs = sids,
	versionIDs = vids,
	outgroup = outgroup,
	ploidy = ploidy,
	wd = runwd,
	orthofinderInBlk = TRUE, 
	overwrite = FALSE, 
	verbose = TRUE,
	nCores = nCores,
	minPepLen = 50,
	gffString = ".gff",
	pepString = ".fa",
	path2orthofinder = "orthofinder",
	path2diamond = "diamond",
	path2mcscanx = path2mcscanx,
	rawGenomeDir = file.path(runwd, "rawGenomes")
)

#Parse and format the annotations
parse_annotations(
	gsParam = gpar,
	gffEntryType = "mRNA",
	gffIdColumn = "Name",
	gffStripText = "Name=",
	headerEntryIndex = 1,
	headerSep = " ",
	headerStripText = ""
)

#Run orthofinder
gpar <- run_orthofinder(gsParam = gpar)

#Set synteny parameters
gpar <- set_syntenyParams(gsParam = gpar)

#Run synteny
gpar <- synteny(gsParam = gpar)

#Call the pangenome
pg <- pangenome(gpar, genomeIDs = c("L1","S1"))

#Save the workspace
save.image("genespace.RData")




