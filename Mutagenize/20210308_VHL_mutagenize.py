#PROGRAM FOR GENERATING OLIGOS FROM 2021-03-08 VHL_sge_order text file
#text file at /camp/home/findlag/home/shared/projects/SGE/VHL/oligo_design/20210308_sge_order.txt
#to be run in directory: /camp/home/findlag/home/shared/projects/SGE/VHL/oligo_design
'''
python ~/home/users/findlag/bin/20210308_VHL_mutagenize.py /camp/home/findlag/home/shared/projects/SGE/VHL/oligo_design/20210308_sge_order.txt /camp/home/findlag/home/shared/projects/SGE/VHL/oligo_design/20210308_VHL_SGE_oligos.csv /camp/home/findlag/home/shared/projects/SGE/VHL/oligo_design/20210308_VHL_SGE_oligos_padded.csv
'''

#Second "padded" output file to be 200mers, padded with fixed sequence on the 5' end for you (note script is part hardcoded for adapters)
#provide the sge_order_worksheet file above as sys.argv[1]
#define an out.file as sys.argv[2]
#third output -- sys.argv[3] is the padded file

import sys
import random

#FUNCTIONS

def mutagenize(my_string):
	upper_string = my_string.upper()
	#the first element will be the wild-type sequence
	all_variants = [upper_string]
	for i in range (0,len(upper_string)):
		for j in ("A", "C", "G", "T"):
			if j != upper_string[i]:
				all_variants += [upper_string[:i]+j+upper_string[i+1:]]
	return all_variants

def codon_mutagenize(my_string):
	upper_string = my_string.upper()
	all_possible_codons = []
	for j in ("A", "C", "G", "T"):
		for k in ("A", "C", "G", "T"):
			for l in ("A", "C", "G", "T"):
				all_possible_codons+=[j+k+l]

	#initialize the list by adding the wt sequence
	all_codon_variants = [upper_string]

	for i in range (0,len(upper_string)):
		if i % 3 == 0:
			wt_codon = upper_string[i:i+3]
			all_possible_mutants = all_possible_codons[:]
			all_possible_mutants.remove(wt_codon)
			for mutant in all_possible_mutants:
				all_codon_variants += [upper_string[:i]+mutant+upper_string[i+3:]]
	return all_codon_variants

#DOES NOT ADD WT SEQUENCE (assumes it will come from elsewhere)
def five_bp_del_mutagenize(my_string):
	upper_string = my_string.upper()
	all_del_variants = []
	del_size = 5
	for i in range (0,len(upper_string)-del_size+1):
		all_del_variants += [upper_string[:i]+upper_string[i+del_size:]]
	#insert a 2nd wt oligo at the end of the list
	return all_del_variants

#DOES NOT ADD WT SEQUENCE (assumes it will come from elsewhere)
def one_bp_del_mutagenize(my_string):
	upper_string = my_string.upper()
	all_del_variants = []
	del_size = 1
	for i in range (0,len(upper_string)-del_size+1):
		all_del_variants += [upper_string[:i]+upper_string[i+del_size:]]
	#insert a 2nd wt oligo at the end of the list
	return all_del_variants

def three_bp_del_mutagenize(my_string):
	upper_string = my_string.upper()
	all_del_variants = []
	del_size = 3
	for i in range (0,len(upper_string)-del_size+1):
		all_del_variants += [upper_string[:i]+upper_string[i+del_size:]]
	#insert a 2nd wt oligo at the end of the list
	return all_del_variants

#This function will take a string and pad it to a specified length with random, non-repetitive sequence on the 5' end!!!:
def pad_to_length(my_string,final_len):
	upper_string = my_string.strip().upper()
	str_len = len(upper_string)
	final_len = int(final_len)
	pad_len = final_len - str_len
	random_pad_seq = 'AT'
	previous_two = ''
	for i in (0,pad_len):
		my_base = 'N'

		previous_two = random_pad_seq[-2:]
		print previous_two, random_pad_seq, i
		my_base = 'N'
		while (my_base == 'N') or (my_base == previous_two[0] and my_base == previous_two[1]):
			rand_float = random.randint(0,3)
			if rand_float == 0:
				my_base = 'A'
			elif rand_float == 1:
				my_base = 'C'
			elif rand_float	== 2:
				my_base = 'G'
			elif rand_float == 3:
				my_base = 'T'
		random_pad_seq += my_base

	padded_seq = random_pad_seq + upper_string

	return(padded_seq)


# A file parser to pick out adapters and mutation sequences
#RULES for worksheet formatting:
#1.  Both oligo 'adapters' must be listed before the sequence to be mutated.
#2.  Lines indicating what the next line is need to come immediately before the next line.
		#markings:  < == upstream adapter,  > == downstream adapter,  $ == sequence to mutate all SNVs,  ^ == all 5 bp dels, ~ == all 3 bp dels
#3.  Lines with sequence should not include anything else and will be stripped of the /n character


sge_worksheet = open(sys.argv[1], 'r')
sge_oligos = open(sys.argv[2], 'w')
oligo_number = 0
f_adapt_check = False
r_adapt_check = False
mut_check = False
codon_mut_check = False
del_mut_check = False
bp3_del_mut_check = False

for line in sge_worksheet:
	#check to see if the previous line indicated the next line is either an adapter or a sequence to mutagenize, or neither
	if f_adapt_check:
		f_adapt = line.strip()
		f_adapt_check = False
	elif r_adapt_check:
		r_adapt = line.strip()
		r_adapt_check = False
	elif mut_check:
		mut_seq = line.strip()
		#code here to generate mutagenized oligos (via function) with concatenated adapters
		mutation_list = mutagenize(mut_seq)
		for mutant in mutation_list:
			oligo_number += 1
			sge_oligos.write('GF'+str(oligo_number)+','+f_adapt+mutant+r_adapt+'\n')

		mut_check = False

	elif codon_mut_check:
		codon_mut_seq = line.strip()
		#code here to generate mutagenized codon oligos with concatenated adapters
		codon_mutation_list = codon_mutagenize(codon_mut_seq)
		for codon_mutant in codon_mutation_list:
			oligo_number += 1
			sge_oligos.write('GF'+str(oligo_number)+','+f_adapt+codon_mutant+r_adapt+'\n')

		codon_mut_check = False

	#this will make both 1 bp and 5 bp deletions!!!
	elif del_mut_check:
		del_mut_seq = line.strip()
		#code here to generate kmer deletions in oligos with concatenated adapters
		del_mutation_list = five_bp_del_mutagenize(del_mut_seq)
		del_mutation_list.extend(one_bp_del_mutagenize(del_mut_seq))
		for del_mutant in del_mutation_list:
			oligo_number += 1
			sge_oligos.write('GF'+str(oligo_number)+','+f_adapt+del_mutant+r_adapt+'\n')

		del_mut_check = False

	elif bp3_del_mut_check:
		bp3_del_mut_seq = line.strip()
		#code here to generate kmer deletions in oligos with concatenated adapters
		bp3_del_mutation_list = three_bp_del_mutagenize(bp3_del_mut_seq)
		for bp3_del_mutant in bp3_del_mutation_list:
			oligo_number += 1
			sge_oligos.write('GF'+str(oligo_number)+','+f_adapt+bp3_del_mutant+r_adapt+'\n')

		bp3_del_mut_check = False

	elif line.startswith('<'):
		f_adapt_check = True
	elif line.startswith('>'):
		r_adapt_check = True
	elif line.startswith('$'):
		mut_check = True
	elif line.startswith('%'):
		codon_mut_check = True
	elif line.startswith('^'):
		del_mut_check = True
	elif line.startswith('~'):
		bp3_del_mut_check = True

	else:
		f_adapt_check = False
		r_adapt_check = False		
		mut_check = False
		codon_mut_check = False
		del_mut_check = False
		bp3_del_mut_check = False


print str(oligo_number)+" oligos generated."

sge_oligos.close()

#hardcoded here -- this is from pUC19 vector, followed by PU1 sequence
adapter_seq = 'GATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCTAAATGGCTGTGAGAGAGCTCAG'
padded_out = open(sys.argv[3],'w')

oligo_count = 0
with open(sys.argv[2], 'r') as oligo_file:
	
	for line in oligo_file:
		oligo_count +=1
		oligo_info = line.strip().split(',')
		oligo_number = oligo_info[0]
		oligo_seq = oligo_info[1]
		oligo_len = len(oligo_seq)
		add_len = 200 - oligo_len
		add_seq = adapter_seq[-add_len:]
		final_seq = add_seq + oligo_seq
		padded_out.write(oligo_number+','+final_seq+'\n')

'''
#this writes them all out a second time
with open(sys.argv[2], 'r') as oligo_file:
	duplicate_oligo_count = 0
	for line in oligo_file:
		duplicate_oligo_count += 1
		oligo_info = line.strip().split(',')
		oligo_number = duplicate_oligo_count+oligo_count
		oligo_seq = oligo_info[1]
		oligo_len = len(oligo_seq)
		add_len = 200 - oligo_len
		add_seq = adapter_seq[-add_len:]
		final_seq = add_seq + oligo_seq
		padded_out.write('GF'+str(oligo_number)+','+final_seq+'\n')

#this write them all out a third time!
with open(sys.argv[2], 'r') as oligo_file:
	duplicate_oligo_count = 0
	for line in oligo_file:
		duplicate_oligo_count += 1
		oligo_info = line.strip().split(',')
		oligo_number = duplicate_oligo_count+oligo_count*2
		oligo_seq = oligo_info[1]
		oligo_len = len(oligo_seq)
		add_len = 200 - oligo_len
		add_seq = adapter_seq[-add_len:]
		final_seq = add_seq + oligo_seq
		padded_out.write('GF'+str(oligo_number)+','+final_seq+'\n')

'''








