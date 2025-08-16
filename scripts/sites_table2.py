import sys, os
import collections
import re

def get_final_sites(sites_file, variant_function):
	c = collections.defaultdict(list)
	with open(sites_file) as f:
		for line in f:
			line = line.strip().split('\t')
			chr_id = str(line[0])
			pos = str(line[1])
			coverage = str(line[4])
			strand = str(line[5])
			RT = str(line[6])
			mut = str(line[7])
			name = chr_id + '\t' + pos
			content = strand + '\t'+ coverage + '\t' + RT + "\t" + mut
			c[name].append(content)

	with open(variant_function) as f2:
		for line in f2:
			line = line.strip().split('\t')
			chr_id2 = str(line[2])
			pos2 = str(line[3])
			gene_name = str(line[1])
			region = str(line[0])
			name2 = chr_id2 +'\t' + pos2
			content2 = gene_name + '\t' + region
			c[name2].append(content2)

	with open(sites_file+'.tmp', 'w') as o:
		row = []
		for k,v in c.items():
			row.append(k+'\t'+'\t'.join(v))
		o.write('\n'.join(row))

sites_file = sys.argv[1]
variant_function = sys.argv[2]

#get_final_sites(sites_file, variant_function)

def step2_alu(sites_file,Alu_bed):
	input_file = sites_file+'.tmp'
	command1 = 'awk \'{print $1,$2,$2,$3,$4,$5,$6,$7,$8}\'  '+ input_file +' |tr \' \' \'\\t\' > '+sites_file+'.tmp2'
	command2 = 'bedtools intersect -a '+sites_file +'.tmp2 -b '+Alu_bed +' -loj|cut -f1-9,13-14|sort -k1,1 -k2,2n |uniq > '+sites_file+'.tmp3'
	os.system(command1)
	os.system(command2)

def uniq_alu(sites_file):
	row = []
	d = collections.defaultdict(list)
	input_file = sites_file +'.tmp3'
	with open(input_file) as f4:
		for line in f4:
			line = line.strip().split('\t')
			name3 = str(line[0])+'\t'+str(line[1])+'\t'+str(line[3])+'\t'+str(line[4])+'\t'+str(line[5])+'\t'+str(line[6])+'\t'+str(line[7]) + '\t'+str(line[8])
			content3 = str(line[9]) + '/' + str(line[10])
			d[name3].append(content3)
	with open(sites_file+'.tmp4', 'w') as out:
		for k,v in d.items():
			row.append(k + '\t'+'|'.join(v))
		out.write('\n'.join(row))
		out.write('\n')

Alu_bed = sys.argv[3]
#step2_alu(sites_file,Alu_bed)
#uniq_alu(sites_file)



import numpy as np
import pysam
from collections import Counter
values = 10

def complement(s):
	basecomplement = {
         "A":"T",
          "T":"A",
          "G":"C",
          "C":"G",
          "a":"t",
          "t":"a",
          "g":"c",
          "c":"g",
			"N":"N",}
	letters = list(s)
	letters = [basecomplement[base] for base in letters]
	letters = letters[::-1]

	return ''.join(letters)

def get_site_arround(sites_file, fasta_path, values):
	row = []
	fasta = pysam.FastaFile(fasta_path)
	sequences = []
	with open(sites_file+'.tmp4') as f4:
		for line in f4:
			line = line.strip().split('\t')
			chr_id = line[0]
			strand = line[2]
			start = int(line[1]) - values -1
			end = int(line[1]) + values
			seq = fasta.fetch(chr_id, start, end).upper()
			if '-' in strand:
				seq = complement(seq)
				
			else:
				seq = seq
			seq =str(seq)
			row.append('\t'.join(line)+'\t'+seq[0:int(values+1/2)]+'\t'+seq[int(values+1/2 +1):])
	with open(sites_file+'.tmp5', 'w') as f5:
		f5.write('\n'.join(row))
		f5.write('\n')

def get_encode_id(sites_file,id_file,output):
	input_file = sites_file+'.tmp5'
	id_name = id_file
	e = collections.defaultdict(list)	
	with open(id_name) as f6:
		for line in f6:
			line = line.strip().split('\t')
			e[line[1]].append(line[0])

	row = []
	g = []
	with open(input_file) as f7:
		for line in f7:
			line = line.strip().split('\t')
			if '(' in line[6]:
				genes = line[6].split("(")[0]
			else:
				genes = line[6].split(",")[0]
			if genes in e:
				encode = ','.join(e[genes])
			else:
				encode = '-'
	#		print(encode)
			row.append(encode+'\t'+'\t'.join(line))

	with open(output, 'w') as f8:
		f8.write('\n'.join(row))
		f8.write('\n')





sites_file = sys.argv[1]
variant_function = sys.argv[2]
Alu_bed = sys.argv[3]
genome = sys.argv[4]
id_file = sys.argv[5]
output = sys.argv[6]

get_final_sites(sites_file, variant_function)
step2_alu(sites_file,Alu_bed)
uniq_alu(sites_file)
get_site_arround(sites_file, genome, values)
get_encode_id(sites_file,id_file,output)
