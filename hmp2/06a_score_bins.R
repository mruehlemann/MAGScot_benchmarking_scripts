library(tidyverse)

sample=rev(strsplit(getwd(), split="/")[[1]])[1]

binners_outputs = read.table(paste0(sample,".contigs_to_bin.tsv"), head=F, stringsAsFactors=F)
colnames(binners_outputs) = c("bin","contig","set")

bisco_ms80 = read.table(paste0(sample,".MAGScoT_sm80.refined.contig_to_bin.out"), head=T, stringsAsFactors=F) %>% mutate(set="MAGScoT_ms80")
bisco_ms50 = read.table(paste0(sample,".MAGScoT_sm50.refined.contig_to_bin.out"), head=T, stringsAsFactors=F) %>% mutate(set="MAGScoT_ms50")
bisco_nomerge = read.table(paste0(sample,".MAGScoT_nomerge.refined.contig_to_bin.out"), head=T, stringsAsFactors=F) %>% mutate(set="MAGScoT_nomerge")
allbisco=bind_rows(bisco_ms80, bisco_ms50, bisco_nomerge)
colnames(allbisco) = c("bin","contig","set")

metawrap = read.table(paste0(sample,".metawrap.contigs_to_bin.tsv"), head=F, stringsAsFactors=F)
colnames(metawrap) = c("bin","contig","set")

das_tool = read.table(paste0(sample,".dastool.contig_to_bin.out"), head=F, stringsAsFactors=F)
colnames(das_tool) = c("bin","contig","set")

allbins = bind_rows(binners_outputs, allbisco, metawrap, das_tool) %>% mutate(bin=gsub("[.]fasta$","",bin))

contigs_lengths=read.table(paste0(sample,".contig_lengths.tsv"), head=F, stringsAsFactors=F)
colnames(contigs_lengths) = c("contig","length")

keepbins = left_join(allbins, contigs_lengths) %>% group_by(set, bin) %>% summarize(totlen=sum(length)) %>% filter(totlen>=250000) %>% pull(bin)

allbins_final = allbins %>% filter(bin %in% keepbins)

write.table(allbins_final, paste0(sample,".final_ctb.tsv"), col.names=F, row.names=F, quote=F, sep="\t")
