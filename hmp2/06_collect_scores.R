library(dplyr)
library(lubridate)

pipes=list()
pipes[["MAGScoT_ms50"]]=c("preproc_prodigal", "hmm","MAGScoT")
pipes[["MAGScoT_ms80"]]=c("preproc_prodigal", "hmm","MAGScoT_ms80")
pipes[["MAGScoT_nomerge"]]=c("preproc_prodigal", "hmm","MAGScoT_nomerge")
pipes[["dastool"]]=c("preproc_prodigal", "dastool")
pipes[["metawrap"]]=c("metawrap_AB", "metawrap_ABCD")

#### get runtimes
samples=list.files()
samples=samples[grepl("hmp|raw|meta|out", samples)==F]

timelogs=sapply(samples, function(s){
  print(s)
  logfiles=list.files(path=paste0("",s), pattern="*.log$")
  logs=sapply(paste0("", s,"/",logfiles), read.table, stringsAsFactors=F, head=F, simplify=F) %>% do.call("rbind",.)
  steps = sapply(unique(logs$V2), function(x) (logs %>% filter(V2 == x, V3=="end") %>% pull(V4) %>% as.character %>% parse_date_time(., orders="%Y%m%d%H%M%S") %>% tail(1) %>% seconds) - (logs %>% filter(V2 == x, V3=="start") %>% pull(V4) %>% as.character %>% parse_date_time(., orders="%Y%m%d%H%M%S") %>% tail(1) %>% seconds) ) %>% data.frame(time_in_s=.) %>% tibble::rownames_to_column("step")
  sapply(names(pipes), function(x) steps %>% mutate(step_this=match(step, pipes[[x]]), pipe=x) %>% arrange(step_this) %>% filter(!is.na(step_this)), simplify=F) %>% do.call("rbind",.) %>% mutate(sample=paste0("",s))
}, simplify=F) %>% do.call("rbind", .)

timelogs %>% group_by(sample, pipe)  %>% summarize(total=sum(time_in_s)) %>% data.frame %>% reshape2::dcast(pipe ~ sample)

##### get all binscores

all_bins = sapply(samples, function(s) read.table(paste0("", s,"/",s,".final_ctb.tsv"), stringsAsFactors=F, head=F) %>% mutate(sample=paste0("",s)),  simplify=F) %>%   do.call("rbind",.) %>% select(V1, V3, sample) %>% unique()

all_counts = read.table("hmpgut.counts.tsv", head=F, stringsAsFactors=F)
colnames(all_counts) = c("feature","sample","count")

all_bisco = sapply(samples, function(s) read.table(paste0("", s,"/",s,".final.scores.out"), stringsAsFactors=F, head=T) %>% mutate(sample=paste0("",s)),  simplify=F) %>%   do.call("rbind",.) %>%
mutate(completeness = ifelse(score.gtdb_rel207_ar53==max & !is.na(score.gtdb_rel207_ar53), uniqueSCGs.gtdb_rel207_ar53/52, uniqueSCGs.gtdb_rel207_bac120/120),
      contamination = ifelse(score.gtdb_rel207_ar53==max & !is.na(score.gtdb_rel207_ar53), multipleSCGs.gtdb_rel207_ar53/52, multipleSCGs.gtdb_rel207_bac120/120),
      score= completeness - 0.5*(contamination))

all_checkm = sapply(samples, function(s) read.table(paste0("", s,"/",s,".final.checkm.out"), head=T, comment.char="", stringsAsFactors=F, sep="\t"), simplify=F) %>%  do.call("rbind",.) %>%
        mutate(completeness= Completeness/100, contamination=Contamination/100, score =completeness - 0.5*(contamination)) %>% mutate(set=all_bins[match(Bin.Id, all_bins$V1),"V3"])


all_out = list(timelogs=timelogs, all_bisco=all_bisco, all_checkm=all_checkm, all_bins=all_bins, all_counts=all_counts)

saveRDS(all_out, "hmpgut.all_out.Rds")
