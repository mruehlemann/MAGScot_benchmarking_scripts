  library(dplyr)
  library(lubridate)

  pipes=list()
  pipes[["MAGScoT_ms50"]]=c("preproc_prodigal", "hmm","MAGScoT")
  pipes[["MAGScoT_ms80"]]=c("preproc_prodigal", "hmm","MAGScoT_ms80")
  pipes[["MAGScoT_nomerge"]]=c("preproc_prodigal", "hmm","MAGScoT_nomerge")
  pipes[["dastool"]]=c("preproc_prodigal", "dastool")
  pipes[["metawrap"]]=c("metawrap_AB", "metawrap_CD")

  #### get runtimes
  timelogs=sapply(0:9, function(s){
    print(s)
    logfiles=list.files(path=paste0("sample_",s), pattern="*.log$")
    logs=sapply(paste0("sample_", s,"/",logfiles), read.table, stringsAsFactors=F, head=F, simplify=F) %>% do.call("rbind",.)
    steps = sapply(unique(logs$V2), function(x) (logs %>% filter(V2 == x, V3=="end") %>% pull(V4) %>% as.character %>% parse_date_time(., orders="%Y%m%d%H%M%S") %>% tail(1) %>% seconds) - (logs %>% filter(V2 == x, V3=="start") %>% pull(V4) %>% as.character %>% parse_date_time(., orders="%Y%m%d%H%M%S") %>% tail(1) %>% seconds) ) %>% data.frame(time_in_s=.) %>% tibble::rownames_to_column("step")
    sapply(names(pipes), function(x) steps %>% mutate(step_this=match(step, pipes[[x]]), pipe=x) %>% arrange(step_this) %>% filter(!is.na(step_this)), simplify=F) %>% do.call("rbind",.) %>% mutate(sample=paste0("sample_",s))
  }, simplify=F) %>% do.call("rbind", .)

  timelogs %>% group_by(sample, pipe)  %>% summarize(total=sum(time_in_s)) %>% data.frame %>% reshape2::dcast(pipe ~ sample)

  ##### get all binscores

  sets=c("gsa2k","Gold standard","bisocreto_ms50","bisocreto_ms80","MAGScoT_nomerge","dastool","metawrap","concoct","maxbin2","metabat2","vamb")

  all_bins = sapply(0:9, function(s) read.table(paste0("sample_", s,"/sample_",s,".final_ctb.tsv"), stringsAsFactors=F, head=F) %>% mutate(sample=paste0("sample_",s)),  simplify=F) %>%   do.call("rbind",.) %>% select(V1, V3, sample) %>% unique()

  all_counts = read.table("marine.counts.tsv", head=F, stringsAsFactors=F)
  colnames(all_counts) = c("feature","sample","count")

  all_amber=sapply(sets, function(s) read.table(paste0("amber_final_out/genome/", s,"/metrics_per_bin.tsv"), head=T, sep="\t", stringsAsFactors=F) %>% mutate(set=s), simplify=F)  %>%  do.call("rbind",.) %>%
    mutate(completeness = Completeness..bp., contamination = (1-Purity..bp.),  score= completeness - 0.5*(contamination))

  all_bisco = sapply(0:9, function(s) read.table(paste0("sample_", s,"/sample_",s,".final.scores.out"), stringsAsFactors=F, head=T) %>% mutate(sample=paste0("sample_",s)) ,  simplify=F) %>%  do.call("rbind",.) %>%
    mutate(completeness = ifelse(score.gtdb_rel207_ar53==max & !is.na(score.gtdb_rel207_ar53), uniqueSCGs.gtdb_rel207_ar53/52, uniqueSCGs.gtdb_rel207_bac120/120),
          contamination = ifelse(score.gtdb_rel207_ar53==max & !is.na(score.gtdb_rel207_ar53), multipleSCGs.gtdb_rel207_ar53/52, multipleSCGs.gtdb_rel207_bac120/120),
          score=completeness - 0.5*(contamination), completeness - 0.5*(contamination))

  ## MISSING: all checkm
  all_checkm = sapply(0:9, function(s) read.table(paste0("sample_", s,"/sample_",s,".final.checkm.out"), head=T, comment.char="", stringsAsFactors=F, sep="\t"), simplify=F) %>%  do.call("rbind",.) %>%
    mutate(completeness= Completeness/100, contamination=Contamination/100, score=completeness - 0.5*(contamination))  %>% mutate(set=all_bins[match(Bin.Id, all_bins$V1),"V3"])

  all_out = list(timelogs=timelogs, all_amber=all_amber, all_bisco=all_bisco, all_checkm=all_checkm, all_bins = all_bins, all_counts = all_counts)

  saveRDS(all_out, "marine.all_out.Rds")
