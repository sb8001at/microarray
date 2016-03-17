def geneidentify(directory, filename)
	# Use text file made with expected-pattern test, that only have probes.
	input = directory + filename + ".txt"
	output = directory + filename + ".csv"
	file_in = File.open(input, "r") # input file
	file_out = File.open(output, "w") # output file, should be CSV
	

	# Taking AGI list from input file
	agi_ary = []
	while agi_list = file_in.gets
		agi_list.chomp!
		agi_list = agi_list.split(/\t/)
		agi_ary.push(agi_list[0])
	end

	file_in.close

	# Gene description import
	file1 = File.open("C:/Users/takao/Dropbox/gene_description_20140101.csv", "r")
	
	# Initializing several hash and arrays.
	gene_hash = Hash.new
	gene2_hash = Hash.new
	syb_hash = Hash.new
	id_agi = []

	# Taking information of genes and AGI codes, and make a hash.
	while id_agi = file1.gets
		id_agi.chomp!
		id_agi = id_agi.split(/,/)
		id_agi2 = id_agi[0]
		syb_hash.store(id_agi2[0..8], id_agi[2])
		gene_hash.store(id_agi2[0..8], id_agi[4])
		gene2_hash.store(id_agi2[0..8], id_agi[3])
	end

	file1.close

	# Taking informations of genes from each hash.
	agi_ary.each do |agi|
	if syb_hash[agi].nil? then
		syb = "no data"
	else
		syb = syb_hash[agi]
	end
	
	if gene_hash[agi].nil? then
		gene = "no data"
	else
		gene = gene_hash[agi]
	end
	
	if gene2_hash[agi].nil? then
		gene2 = "no data"
	else
		gene2 = gene2_hash[agi]
	end
	
	temp = [agi, syb, gene, gene2]
	
	file_out.print(agi,",",syb,",",gene, ",", gene2, "\n")
	
	end
	
	file_out.close

end
