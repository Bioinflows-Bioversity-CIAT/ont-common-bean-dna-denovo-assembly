rm(list = ls())

pks <- c('dplyr','stringr','argparse')

pks <- suppressPackageStartupMessages(sapply(pks, require, character.only = TRUE))

if (any(!pks)) stop('The package(s) ', names(pks)[!pks], ' is/are not available for load.')

pr <- ArgumentParser()

pr$add_argument("-blastn", type = "character", metavar = 'file',
                    help = "Path to blastp (fmt6) output file")
pr$add_argument("-cutoff", type = "integer", metavar = 'file',
                    help = "Allignment length cutoff")
pr$add_argument("-map", type = "character", metavar = 'file',
                    help = "Path to map with genetic positions file")
pr$add_argument("-out", type = "character", metavar = "path",
                    help = "Outputh path")                    
ar = pr$parse_args()

# Read genetic map [locus_id(str), lg(int), gen_pos(float)]
map <- read.table(ar$map, sep = '\t', header = T)
print(head(map))
# Read blastn  fmt6 output
# 'qseqid','sseqid','pident','length',
# 'mismatch','gapopen','qstart','qend',
# 'sstart','send','evalue','bitscore'
blastn <- read.table(ar$blastn, sep = '\t', header = F)
colnames(blastn) <- c('qseqid','sseqid','pident','length',
                      'mismatch','gapopen','qstart','qend',
                      'sstart','send','evalue','bitscore')
print(head(blastn))
# Select most significative and filter by alligment length
blast_filt <- blastn %>% 
  group_by(qseqid) %>% 
  slice_min(evalue) %>% 
  ungroup() %>% 
  distinct() %>% 
  filter(length > ar$cutoff) %>% 
  rename(locus_id = qseqid)

# Create the all_maps genetic map input file
# all_maps_input <- merge(map, blast_filt, by = 'locus_id') %>% 
#   rename(
#     ScaffoldID = sseqid,
#     Linkage = lg,
#     GeneticPosition = gen_pos) %>% 
#   mutate(
#     ScaffoldPosition = round((sstart + send)/2, 0)
#   ) %>%
#   select(ScaffoldID, ScaffoldPosition, Linkage, GeneticPosition)

print(head(blast_filt))


all_maps_input <- merge(blast_filt, map, by = "locus_id") %>% 
    mutate(
      lg_tag = glue::glue("lg{lg}:{gen_pos}"),
      avg_pos = round((sstart + send)/2, 0),
      i_pos = avg_pos,
      e_pos = avg_pos + 1,
      ) %>%
    group_by(lg, gen_pos) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    arrange(lg, gen_pos) %>%
    select(sseqid , i_pos, e_pos, lg_tag)

write.table(all_maps_input, file = ar$out,
    quote = F, row.names = F, col.names = F, sep = '\t')
