import sys

def neg(file_neg,strand):
	#row = []
	with open(file_neg) as f:
		for line in f:
			line = line.strip().split('\t')
			if len(line) ==1:
				continue
			if line[3] == 'T':
				ref_base = 'A'
			elif line[3] == 'A':
				ref_base = 'T'
			elif line[3] == 'C':
				ref_base = 'G'
			elif line[3] == 'G':
				ref_base = 'C'
			elif line[3] == 'N':
				ref_base = 'N'
			if line[4] == 'T':
				alt_base = 'A'
			elif line[4] == 'A':
				alt_base = 'T'
			elif line[4] == 'C':
				alt_base = 'G'
			elif line[4] == 'G':
				alt_base = 'C'
			elif line[4] == 'N':
				alt_base = 'N'
			abs_pos = int(line[8])+1 - int(line[7])
			print(line[0] +'\t'+ line[1]+ '\t'+line[2]+'\t'+ref_base +'\t'+ alt_base +'\t'+line[5] +'\t'+line[6]+'\t'+str(abs_pos) +'\t'+line[8]+'\t'+line[9] +'\t'+line[10]+'\t'+line[11] +'\t'+line[12]+'\t'+strand)
#			print(row)

def pos(file_pos,strand):
	row = []
	with open(file_pos) as f:
		for line in f:
			line = line.strip().split('\t')
			if len(line) ==1:
				continue
			print('\t'.join(line)+'\t'+strand)


file_pos = sys.argv[1]
file_neg = sys.argv[2]

pos(file_pos,'pos')
neg(file_neg,'neg')

#print(row_pos)
#print(row_neg)
