=begin
With this program, we can take AGI and gene informations from set of gene symbols.
However, it's difficult to estimate AGI codes of all genes, because there are several
probes that address one gene. Small amount of information of AGI codes also affects
the difficulty.

Before running this program, visit http://www.arabidopsis.org/tools/bulk/microarray/index.jsp
and save information of probes, AGI, genes as text file.
=end

# Use text file made with expected-pattern test, that only have probes.
file4 = File.open("C:/users/takao/dropbox/salicylic acid related genes.txt", "r") # input file
file = File.open("C:/users/takao/dropbox/salicylic acid related genes.csv", "w") # output file, should be CSV

# Move AGI on left column and probe on right colomn in ATH1_AGI_rev.csv with Excel.
file1 = File.open("C:/Users/takao/Dropbox/Microarray_analysis/ATH1_AGI_rev.csv", "r")

# Save probe (left) and gene (or symbol information, right) in Excel as two csv files.
file2 = File.open("C:/Users/takao/Dropbox/Microarray_analysis/ATH1_gene.csv", "r")

# To take the list of symbols, conduct R program to aquire the information.
file3 = File.open("C:/Users/takao/Dropbox/Microarray_analysis/ATH1_syb.csv", "r")


# Initializing several hash and arrays.
agirev_hash = Hash.new
gene_hash = Hash.new
syb_hash = Hash.new
id_agi2 = []
id_gene2 = []
id_syb2 = []

# Taking information of probes and AGI codes, and make a hash.
while id_agi = file1.gets
id_agi.chomp!
id_agi2 = id_agi.split(/,/)
id_agi3 = id_agi2[0]
agirev_hash.store(id_agi3[0..8], id_agi2[1])
end

file1.close

# Taking information of probes and genes, and make a hash.
while id_gene = file2.gets
id_gene.chomp!
id_gene2 = id_gene.split(/,/)
gene_hash.store(id_gene2[0], id_gene2[1])
end

file2.close

# Taking information of probes and gene symbols, and make a hash.
while id_syb = file3.gets
id_syb.chomp!
id_syb2 = id_syb.split(/,/)
syb_hash.store(id_syb2[0], id_syb2[1])
end

file3.close

agi_ary = []

# Taking AGI list from input file
while agi_list = file4.gets
agi_list.chomp!
agi_list = agi_list.split(/\t/)
agi_ary.push(agi_list[0])
end

file4.close

syb_ary = Array.new
gene_ary = Array.new

# Taking informations of genes from each hash.
agi_ary.each do |agicode|
syb_ary.push(syb_hash[agirev_hash[agicode]])
gene_ary.push(gene_hash[agirev_hash[agicode]])
end

# Printing results on output file.
i = 0
while i < agi_ary.length
file.print(agi_ary[i],",",syb_ary[i],",",gene_ary[i],"\n")
i = i + 1
end

file.close

