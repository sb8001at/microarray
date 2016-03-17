# Use text file made with probes.
file2 = File.open("C:/Users/public/documents/microarray/PR excision/0_1hour/bZIP.txt", "r") # input file
file = File.open("C:/Users/public/documents/microarray/PR excision/0_1hour/bZIP_listseq.txt", "w") # output file, should be CSV

# Move AGI on left column and sequence on right colomn in Tair10_CDS.txt with Excel.
file1 = File.open("C:/Users/takao/Dropbox//Tair10_CDS.txt", "r")

# Initializing several hash and arrays.
agi_seq_hash = Hash.new
agi_ary = Array.new
id_agi2 = Array.new
id_gene2 = Array.new
id_syb2 = Array.new

# Taking information of probes and AGI codes, and make a hash.
while id_agi = file1.gets
	id_agi.chomp!
	id_agi2 = id_agi.split(/\t/)
	id_agi3 = id_agi2[0]
	id_agi4 = id_agi3[0..8]
	agi_seq_hash.store(id_agi4, id_agi2[1])
end

file1.close

# Taking AGI list from input file
while agi_list = file2.gets
agi_list.chomp!
agi_list = agi_list.split(/\t/)
agi_ary.push(agi_list[0])
end

file2.close

seq_ary = Array.new

# Taking informations of genes from each hash.
agi_ary.each do |agicode|
seq_ary.push(agi_seq_hash[agicode])
end

# Printing results on output file.
i = 0
while i < agi_ary.length
	if seq_ary[i] != nil
		file.print('>', agi_ary[i], "\n", seq_ary[i], "\n")
	end
i = i + 1
end

file.close

