import sys
import collections
from collections import Counter

def get_AG(site_file):
	c = collections.defaultdict(lambda :0)
	with open(site_file) as f:
		for line in f:
			line = line.strip().split('\t')
			chrID = line[0]
			position = line[1]
			ref_base = line[2]
			alt_base = line[3]
			coverage = line[4]
			strand = line[5]
			cleave_number = line[6]
			mut_number = line[7]
#				relative_pos = line[7]
#				read_length = line[8]
#				base_qual = line[9]
#				read_name = line[11]
#				strand = line[13]

#				if int(line[8]) - int(line[7]) == 1 and float(line[9]) > 26 and line[3] == 'A' and line[4] == 'G':
			if strand == 'pos' and ref_base == 'A' and alt_base == 'G':
				content = str(chrID)+'\t'+str(position)+'\t'+ str(ref_base)+'\t'+ str(alt_base)+ '\t'+ str(coverage)+'\t'+ "+\t" +cleave_number+"\t"+mut_number
			elif strand == 'neg' and ref_base == 'A' and alt_base == 'G':
				content = str(chrID)+'\t'+str(position)+'\t'+ 'T\tC\t'+ str(coverage)+'\t'+ "-\t" +cleave_number+"\t"+mut_number
			else:
				continue
			print(content)

print("chrom\tposition\tref_base\talt_base\tcoverage\tstrand\ttrunc_reads\tmut_reads")
#for i in sys.argv[1:]:
get_AG(sys.argv[1])
