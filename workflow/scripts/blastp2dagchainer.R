rm(list = ls())

pks <- c('dplyr','stringr','argparse')

pks <- suppressPackageStartupMessages(sapply(pks, require, character.only = TRUE))

if (any(!pks)) stop('The package(s) ', names(pks)[!pks], ' is/are not available for load.')

pr <- ArgumentParser()

pr$add_argument("-blastp", type = "character", metavar = 'file',
                    help = "Path to blastp (fmt6) output file")

pr$add_argument("-ref_dic", type = "character", metavar = "file", 
                    help = "Path to reference gene dictionary")

pr$add_argument("-qry_dic", type = "character", metavar = "file", 
                    help = "Path to query gene dictionary")

pr$add_argument("-p", type = "character", metavar = "string",
                    help = "Reference gene id suffix to remove", default = ".1.p")

pr$add_argument("-out", type = "character", metavar = "path",
                    help = "Outputh path")                    

ar = pr$parse_args()


# Read the gene directory,
# group by Chr/contig, order by position,
# get the gene number along each group.
get_gene_order <- function(df, chrom, pos_i){
  df %>%
    group_by({{chrom}}) %>%
    arrange({{pos_i}}) %>%
    mutate(
      gene_ord = row_number()
    ) %>%
    ungroup()
}

# Merge blastp output with gene directory
# it generalizes the ref and query merge 
# to join the gene order and the e-value 
# into the same df
merge_blasp_gene <- function(blasp,
                             genes,
                             blastp_id,
                             genes_id,
                             contig_colname,
                             gene_ord_colname
                             ){
  
  blastp_merged <- merge(blasp,
                      genes[, c("contig", "gene_id", "gene_ord")],
                      by.x = blastp_id,
                      by.y = genes_id) %>%
    rename(
      {{contig_colname}} := contig,
      {{gene_ord_colname}} := gene_ord)
  
  return(blastp_merged)
}


blastp2dagchainer <- function(blastp,
                              r_genes,
                              q_genes,
                              ref_suffix = ".1.p"){
  
  # Replace string suffix of reference gene ids to match
  # the id in reference gene directory
  blastp <- blastp %>%
    mutate(
      sseqid = str_replace(sseqid, ref_suffix, "")
    )
  
  # merge blastp with reference gene directory
  blast_ref <- merge_blasp_gene(blasp = blastp,
                                genes = r_genes,
                                blastp_id = "sseqid",
                                genes_id = "gene_id",
                                contig_colname = "r_contig",
                                gene_ord_colname = "r_gene_ord")
  
  # merge blastp with query gene directory
  blast_ref_qry <- merge_blasp_gene(blasp = blast_ref,
                                    genes = q_genes,
                                    blastp_id = "qseqid",
                                    genes_id = "gene_id",
                                    contig_colname = "contig",
                                    gene_ord_colname = "gene_ord")
  
  # Duplicate gene order columns and order 
  # df as dagchainer expects
  dagchainer_fmt <- blast_ref_qry %>%
    mutate(
      r_gene_ord_2 = r_gene_ord,
      gene_ord_2 = gene_ord,
    ) %>%
    select(
      r_contig,
      sseqid,
      r_gene_ord,
      r_gene_ord_2,
      contig,
      qseqid,
      gene_ord,
      gene_ord_2,
      evalue
    ) %>%
    group_by(sseqid, qseqid) %>% 
    slice_min(evalue) %>%                
    ungroup() %>%
    distinct() %>%
    arrange(r_contig,
            r_gene_ord)
  return(dagchainer_fmt)
}

# Function call

# Query dictionary
q_df <- read.table(ar$qry_dic, sep = "\t", header = F)
colnames(q_df) <- c("contig", "pos_i", "pos_f", "gene_id")
q_df_dict <- get_gene_order(q_df, contig, pos_i)

# Reference dictionary
r_df <- read.table(ar$ref_dic, sep = "\t", header = F)
colnames(r_df) <- c("contig", "pos_i", "pos_f", "gene_id")
r_df_dict <- get_gene_order(r_df, contig, pos_i)

# Read blastp output
blastp <- read.table(ar$blastp, sep = "\t", header = F)
colnames(blastp) <- c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')

# Convert to DAGChainer format
out <- blastp2dagchainer(blastp, r_df_dict, q_df_dict, ref_suffix = ar$p)

write.table(out, ar$out, sep = '\t', row.names = F, col.names = F, quote = F)
