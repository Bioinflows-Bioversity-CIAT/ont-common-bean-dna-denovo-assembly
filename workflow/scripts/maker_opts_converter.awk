BEGIN {
	est = "est=" est_path;
	protein = "protein=" protein_path;
	genome = "genome=" genome_path;
	species = "augustus_species=" species_name;

}
{
	gsub("^est=", est);
	gsub("^protein=", protein);
	gsub("^model_org=all", "model_org= ");
	gsub("^genome=", genome);
	gsub("^augustus\_species\=", species);
	gsub("^keep\_preds\=0 ","keep\_preds\=1");
 	print
}
END {
	# pass
}
