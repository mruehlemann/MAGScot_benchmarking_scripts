library(tidyverse)

samples=paste0("sample_",0:9)

all_scores_binners = sapply(samples, function(s) read.table(paste0(s,"/",s,".biscoreto.scores.out"), head=T, stringsAsFactors=F) %>% mutate(sample=s) %>% filter(grepl("Bisco", set)==F), simplify=F)
all_scores_biscoreto = sapply(samples, function(s) read.table(paste0(s,"/",s,".biscoreto.refined.out"), head=T, stringsAsFactors=F) %>% mutate(sample=s, set="BiScoReTo"), simplify=F)
all_scores_dastool = sapply(samples, function(s) read.table(paste0(s,"/",s,".dastool.scores.out"), head=T, stringsAsFactors=F) %>% mutate(sample=s), simplify=F)
all_scores_gsa = sapply(samples, function(s) read.table(paste0(s,"/",s,".gsa.scores.out"), head=T, stringsAsFactors=F) %>% mutate(sample=s), simplify=F)
all_scores_gsa2k = sapply(samples, function(s) read.table(paste0(s,"/",s,".gsa_2k.scores.out"), head=T, stringsAsFactors=F) %>% mutate(sample=s), simplify=F)

#all_scores_biscoreto = sapply(samples, function(s) read.table(paste0(s,"/",s,".biscoreto.refined.out"), head=T, stringsAsFactors=F) %>% mutate(sample=s), simplify=F)

all_scores = bind_rows(do.call("rbind", all_scores_binners), do.call("rbind", all_scores_biscoreto), do.call("rbind", all_scores_dastool), do.call("rbind", all_scores_gsa), do.call("rbind", all_scores_gsa2k))

saveRDS(all_scores, "all_scores.Rds")

all_scores %>% filter(max>=0.5) %>% arrange(set, -max) %>% group_by(set) %>% mutate(bin_rank = row_number()) %>% ggplot(aes(x=bin_rank, y=max, group=set)) + geom_line(aes(col=set))




runtimes_biscoreto = sapply(samples, function(s) read.table(paste0(s,"/",s,".biscoreto.log"), head=F, stringsAsFactors=F)  %>% reshape2::dcast(V2 ~ V3, value.var="V4") %>% mutate(set="BiScoReTo", sample=s), simplify=F)

runtimes_dastool = sapply(samples, function(s) read.table(paste0(s,"/",s,".dastool.log"), head=F, stringsAsFactors=F)  %>% reshape2::dcast(V2 ~ V3, value.var="V4") %>% mutate(set="DASTool", sample=s), simplify=F)

#runtimes_metawrap = sapply(samples, function(s) read.table(paste0(s,"/",s,".metawrap.log"), head=F, stringsAsFactors=F)  %>% reshape2::dcast(V2 ~ V3, value.var="V4") %>% mutate(set="metaWRAP", time = end - start, sample=s), simplify=F)
cc=sapply(as.character(c$V4), function(x) strptime(x, format="%Y%m%d%H%M%S"), simplify=F)

all_times = bind_rows(do.call("rbind", runtimes_biscoreto), do.call("rbind", runtimes_dastool))


#### amber output

subsets=c("Gold standard", "gsa2k", "biscoreto", "biscoreto_nomerge","concoct","dastool","maxbin2","metabat2","vamb")
all_amber = sapply(subsets, function(s) read.table(paste0("amber_results/genome/",s,"/metrics_per_bin.tsv"), head=T, stringsAsFactors=F, sep="\t") %>% mutate(set=s), simplify=F) %>% do.call("rbind",.)
all_amber = all_amber %>% mutate(score = Completeness..bp. - (1-Purity..bp.))
all_amber %>% filter(Bin.size..bp.>=1000000) %>% group_by(set) %>% summarize(n=n(), Purity = mean(Purity..bp.), Completeness=mean(Completeness..bp.), Score=mean(score))
all_amber %>% filter(Bin.size..bp.>=1000000) %>% group_by(set) %>% top_n(n=500, wt=score)  %>% summarize(n=n(), Purity = mean(Purity..bp.), Completeness=mean(Completeness..bp.), Score=mean(score))

saveRDS(all_amber, "all_amber.Rds")
