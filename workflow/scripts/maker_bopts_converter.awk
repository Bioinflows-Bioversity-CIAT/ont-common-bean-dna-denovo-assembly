BEGIN {

}
{
	gsub("^depth_blastx\=", "depth_blastx\=10");
	gsub("^depth_blastn\=", "depth_blastn\=10");
  gsub("^depth_tblastx\=", "depth_tblastx\=10");
	
 	print
}
END {
	# pass
}
