import sys
import collections
from collections import Counter

def Merge(dict1, dict2): 
	res = dict()
	for k,v in dict1.items():
		if k in dict2:
			v = dict2[k] + v
			res[k] = v 
		else:
			v =  v 
			res[k] = v
	return res


def reverse_base(base):
	if base == 'A':
		rev_base = 'T'
	elif base == 'T':
		rev_base = 'A'
	elif base == 'C':
		rev_base = 'G'
	elif base == 'G':
		rev_base = 'C'

	elif base == 'N':
		rev_base = 'N'

	return rev_base

def get_AG(site_file):
	dict_I_low = dict()
	cal_3bp_sites = dict()
	filter_sites = dict()
	dict_rawmutation = dict()
	dict_I = dict()
	c = collections.defaultdict(lambda :0)
	d = collections.defaultdict(lambda :0)
	with open(site_file) as f:
		for line in f:
			line = line.strip().split('\t')
			if len(line) == 1:
				continue
			else:
				chrID = line[0]
				position = line[2]
				ref_base = line[3]
				alt_base = line[4]
				coverage = float(line[5])
				relative_pos = line[7]
				read_length = line[8]
				base_qual = line[9]
				read_name = line[11]
				strand = line[13]

				if int(line[8]) - int(line[7]) == 1 and coverage > 0 and float(line[9]) > 26:# and float(line[9]) > 10:# and line[3] == 'A' and line[4] == 'G':

					content = str(chrID)+'\t'+str(position)+'\t'+ str(ref_base)+'\t'+ str(alt_base)+ '\t'+ str(coverage)+'\t'+strand

					c[content] +=1
				if coverage > 0 and float(line[9]) > 26:

					content = str(chrID)+'\t'+str(position)+'\t'+ str(ref_base)+'\t'+ str(alt_base)+ '\t'+ str(coverage)+'\t'+strand
					d[content] +=1
						
#	dict3 = Merge(c, d)
	
#	print(dict3)
#	print(c)
#	print(d)
	for k,v in c.items():
		ref_base = k.split('\t')[2]
		alt_base = k.split('\t')[3]
		strand = k.split('\t')[5]
		position = k.split('\t')[1]
		chrID = k.split('\t')[0]
		coverage = k.split('\t')[4]
		new_content= str(chrID)+'\t'+str(position)
		new_value = str(ref_base)+'\t'+ str(alt_base)+ '\t' +str(coverage)+'\t'+strand+'\t'+ str(v)
		if v >= 6 and float(v)/float(k.split('\t')[4]) >= 0.1 and float(coverage) < 5000  and 'M' not in chrID:
			dict_I[new_content] = new_value
		if v <6 and v >=3 and float(v)/float(k.split('\t')[4]) >= 0.05 and float(coverage) < 5000 and 'M' not in chrID:
			dict_I_low[new_content] = new_value

	for k,v in d.items():
		ref_base = k.split('\t')[2]
		alt_base = k.split('\t')[3]
		strand = k.split('\t')[5]
		position = k.split('\t')[1]
		chrID = k.split('\t')[0]
		coverage = k.split('\t')[4]
		new_content= str(chrID)+'\t'+str(position)
		new_value = str(ref_base)+'\t'+ str(alt_base)+ '\t' +str(coverage)+'\t'+strand+'\t'+ str(v)

		if v >= 4 and float(v)/float(k.split('\t')[4]) >= 0.1:
			dict_rawmutation[new_content] = new_value
	
	for m,n in dict_I.items():
		chrID = m.split('\t')[0]
		position = m.split('\t')[1]
		p = 0
		for i in range(1,3):
			l1_position = int(position) - i
			r1_position = int(position) + i
			wrong_l1 = chrID + "\t" + str(l1_position)
			wrong_r1 = chrID + "\t" + str(r1_position)
			right_partten = n.split('\t')[0] + "\t" + n.split('\t')[1]
			if wrong_l1 in dict_rawmutation and float(dict_rawmutation[wrong_l1].split("\t")[4]) / float(dict_rawmutation[wrong_l1].split("\t")[2]) > 0.1:
				wrong_l1_partten =  dict_rawmutation[wrong_l1].split("\t")[0]+"\t"+dict_rawmutation[wrong_l1].split("\t")[1]
				if  wrong_l1_partten != right_partten:# or wrong_l1_partten == "A\tA":
					p = 1
			#	else:
			#		filter_sites[m]=n

			if wrong_r1 in dict_rawmutation and float(dict_rawmutation[wrong_r1].split("\t")[4]) / float(dict_rawmutation[wrong_r1].split("\t")[2]) > 0.1:
				wrong_r1_partten =  dict_rawmutation[wrong_r1].split("\t")[0]+"\t"+dict_rawmutation[wrong_r1].split("\t")[1]
				if	wrong_r1_partten != right_partten:# or wrong_r1_partten == "A\tA":
					p = 1
			
		if p == 0:
	#		filter_sites[m]=str(n) +"\t"+str(d[m])
			a = "\t".join((m+"\t"+n).split('\t')[:-1])
			filter_sites[m]=str(n) +"\t"+str(d[a])
	
	for k,v in filter_sites.items():
		chrID = k.split('\t')[0]
		position = k.split('\t')[1]
	#	print(chrID+"\t"+position)
		for i in range(1,3):
			l1_position = int(position) - i
			r1_position = int(position) + i
			wrong_l1 = chrID + "\t" + str(l1_position)
			wrong_r1 = chrID + "\t" + str(r1_position)
			if wrong_l1 in dict_I_low:
				b = wrong_l1 + "\t" + "\t".join(dict_I_low[wrong_l1].split('\t')[:-1])
				cal_3bp_sites[wrong_l1] = dict_I_low[wrong_l1] + "\t" +str(d[b])
			if wrong_r1 in dict_I_low:
				b = wrong_r1 + "\t" + "\t".join(dict_I_low[wrong_r1].split('\t')[:-1])
				cal_3bp_sites[wrong_r1] = dict_I_low[wrong_r1] + "\t" +str(d[b])
	with open(site_file+"_high_sites.txt",'w') as o1:
		for k,v in filter_sites.items():
			o1.write(k+"\t"+v+"\n")
	with open(site_file+"_low_sites.txt",'w') as o2:
		for k,v in cal_3bp_sites.items():
			o2.write(k+"\t"+v+"\n")

	with open(site_file+"_total.sites.txt",'w') as o3:
		for k,v in filter_sites.items():
			o3.write(k+"\t"+v+"\n")
		for k,v in cal_3bp_sites.items():
			o3.write(k+"\t"+v+"\n")

	o1.close()
	o2.close()
	o3.close()
#	print(dict_I_low)
#print("chrom\tposition\tref_base\talt_base\tcoverage\tstrand\ttrunc_reads")
for i in sys.argv[1:]:
	get_AG(i)
