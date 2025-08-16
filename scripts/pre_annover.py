import sys, os
import collections
def get_annvar_input(sites_file,out_file):
	c = collections.defaultdict(list)
	with open(sites_file) as f, \
	open(out_file,'w') as out:
		f.readline()
		for line in f:
			line = line.strip().split('\t')
			chr_id = line[0]
			chr_id.replace("chr", "")
			pos = line[1]
			strand = line[5]
			if strand == '+':
				ref_base = line[2]
				alt_base = line[3]
			elif strand == '-':
				ref_base = 'T'
				alt_base = 'C'
				
			row = chr_id + '\t' + pos+'\t'+pos+'\t'+ref_base+'\t'+alt_base+'\n'
			out.write(row)


sites_file = sys.argv[1]
out_file = sys.argv[2]
get_annvar_input(sites_file,out_file)
