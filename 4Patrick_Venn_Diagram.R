#  Lung COPD project - normal lung tissue comparisons to TP

#  Import Venn diagram library - Vennerable
library(Vennerable)

#  Set working directory
#setwd("Y://pflores/COPD/111713/COM_dds")

setwd("Y:///pflores/COPD/111713/2v8")

#  Import and assign data (from pep.xls file, which is tab-delimited, not excel)
redo <- read.delim("8ulredo.pep.xls", sep="\t", header=T)
redo <- unique(redo$peptide) # want unique peptides

#  TP Lung 2ul
#setwd("R:/tsd/Production/Lung/LongGradient/pat/cometOutput")
tp <- read.delim("2ul_comet.pep.iproph.pep.xls", sep="\t", header=T)
tp <- unique(tp$peptide)

#  Start a png file for saving to
png("Venn_Diagram_RedoLungVSTPLung_protein.png", height=8, width=8, units="in", res=200)

vd <- Venn(Sets = list("Redo Lung (peptides)" = redo,
		"TP Lung (peptides)" = tp))

plot(vd, type="circles", doWeights=T)

# Turn the device off
dev.off()

