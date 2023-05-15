#20220728_VHLdb_parsing.py
#usage in linux running python 2.7:  python /camp/lab/findlayg/home/users/findlag/bin/20220728_VHLdb_parsing.py 
import sys
import os
import subprocess
import time
from datetime import datetime
sys.path.insert(0,'/camp/lab/findlayg/home/users/findlag/bin')
import dict_tools_11192015 as dict_tools

VHLdb_file = '/camp/lab/findlayg/home/shared/projects/SGE/VHL/vhldbMuts_1581968329043_20220728.txt'
VHLdb_trans_simp_file = '/camp/lab/findlayg/home/shared/projects/SGE/VHL/VHLdb_transmission_simplifications.txt'


VHLdb_trans_simp_dict = {}
with open(VHLdb_trans_simp_file, 'r') as VHL_tc:
	for line in VHL_tc:
		sline = line.strip()
		pairing = sline.split('\t')
		VHLdb_trans_simp_dict[pairing[0].strip()] = pairing[1]
VHLdb_trans_simp_dict[''] = 'not_provided'

key_list = VHLdb_trans_simp_dict.keys()
for my_key in key_list:
	VHLdb_trans_simp_dict[my_key+' '] = VHLdb_trans_simp_dict[my_key]



header_line = ''
entry_count = -1


line_dol = {}
with open(VHLdb_file, 'r') as VHLdb:
	for line in VHLdb:
		if entry_count == -1:
			header_line = line
			entry_count+=1
		else:
			entry_count+=1
			line_info = line.strip().split('\t')
			variant = line_info[1]
			line8 = line_info[8].strip()
			line9 = line_info[9].strip('"')
			line_info.append(VHLdb_trans_simp_dict[line9])
			if len(line_info) != 13:
				print len(line_info)
				print line
				time.sleep(2)
			line_dol[entry_count] = line_info
		
'''
new_trans = {}
for my_key in line_dol:
	if line_dol[my_key][12] not in new_trans:
		new_trans[line_dol[my_key][12]] = 1
	else:
		new_trans[line_dol[my_key][12]] = new_trans[line_dol[my_key][12]]+1
print new_trans
'''
#Have read the file

unique_alleles = 0
unique_snvs = 0
variant_dict = {}
single_variants = 0
multi_variants = 0
line8_entries = {}
		
for i in range(1,entry_count+1):
	line_info = line_dol[i]
	variant = line_info[1]
	line8 = line_info[8].strip()
	line9 = line_info[9].strip('"')
	if variant not in variant_dict:
		variant_dict[variant] = line_info
		variant_dict[variant].append(1)
		#print variant_dict[variant]
		unique_alleles+=1
	else:
		#print 'already present'
		#print variant_dict[variant][-1]
		variant_dict[variant][-1] = variant_dict[variant][-1]+1
		variant_dict[variant][8] = variant_dict[variant][8]+': :'+line_info[8]
		variant_dict[variant][9] = variant_dict[variant][9]+': :'+line_info[9]
		variant_dict[variant][10] = variant_dict[variant][10]+': :'+line_info[10]
		variant_dict[variant][11] = variant_dict[variant][11]+': :'+line_info[11]
		variant_dict[variant][12] = variant_dict[variant][12]+': :'+line_info[12]
	if line8 in line8_entries:
			line8_entries[line8]+=1
	else:
			line8_entries[line8] = 1
			
#parse disease column for phenotypic information
total_VHL_count = 0
total_VHL_dis_vars = 0
variant_count = 0
for variant in variant_dict:
	variant_count +=1 
	if variant == '':
		continue
	elif len(variant_dict[variant]) != 14:
		print "error, variant has wrong number of entries:", len(variant_dict[variant])
		print variant_dict[variant]
		time.sleep(1)
	else:
		variant_data = variant_dict[variant]
		variant_disease_info = variant_data[8]
		variant_disease_list = variant_disease_info.split(': :')

		variant_trans_info = variant_data[12]
		variant_trans_list = variant_trans_info.split(': :')
		#count total disease labels
		vhl_disease_count = 0
		t1_count = 0
		t2_count = 0
		t2a_count = 0
		t2b_count = 0
		t2c_count = 0
		#count total tumor labels
		rcc_count = 0
		pheo_count = 0
		hb_count = 0
		r_cyst_count = 0
		p_cyst_count = 0
		pgl_count = 0
		ra_count = 0
		polycyt_count = 0

		germline_count = 0
		somatic_count = 0
		homozygous_count = 0

		for my_entry in variant_disease_list:
			if (("syndrome" in my_entry) or ("Syndrome"in my_entry) or ("VHL Disease" in my_entry) or ("VHL disease" in my_entry) or ("Lindau Disease" in my_entry) or ("Lindau disease" in my_entry) or ("VHl" in my_entry) or ("VHL" in my_entry)):
				vhl_disease_count +=1
				total_VHL_count +=1
			if (("Type 1" in my_entry) or ("type 1" in my_entry) or ("Type One" in my_entry) or ("type one" in my_entry)):
				t1_count +=1
			if (("Type 2" in my_entry) or ("type 2" in my_entry) or ("Type Two" in my_entry) or ("type II" in my_entry)):
				t2_count +=1
			if (("2a" in my_entry) or ("2A" in my_entry) or ("IIA" in my_entry)):
				t2a_count +=1
			if (("2b" in my_entry) or ("2B" in my_entry) or ("IIB" in my_entry)):
				t2b_count +=1
			if (("2c" in my_entry) or ("2C" in my_entry) or ("IIC" in my_entry)):
				t2c_count +=1
			if (("RCC" in my_entry) or ("rcc" in my_entry) or ("renal cell carcinoma" in my_entry)):
				rcc_count +=1
			if (("pheo" in my_entry)):
				pheo_count +=1
			if (("HB" in my_entry) or ("angioblastoma" in my_entry)):
				hb_count +=1
			if (("renal cysts" in my_entry) or ("Renal cysts" in my_entry)):
				r_cyst_count +=1
			if (("pancreatic cysts" in my_entry) or ("Pancreatic cysts" in my_entry) or ("PCT" in my_entry)):
				p_cyst_count +=1
			if (("PGL" in my_entry) or ("araganglioma" in my_entry)):
				pgl_count +=1
			if (("etinal angioma" in my_entry) or ("etinal Angioma" in my_entry)):
				ra_count +=1
			if (("olycythemia" in my_entry)):
				polycyt_count +=1

		for my_entry in variant_trans_list:
			if ("Germline" in my_entry):
				germline_count += 1
			if ("Somatic" in my_entry):
				somatic_count +=1
			if ("Homozygous" in my_entry):
				homozygous_count +=1

		var_quants = [vhl_disease_count,t1_count,t2_count,t2a_count,t2b_count,t2c_count,rcc_count,pheo_count,hb_count,r_cyst_count,p_cyst_count,pgl_count,ra_count,polycyt_count,germline_count,somatic_count,homozygous_count]

		type1_only = "False"
		type1_or_2b_only = "False"
		type2_only = "False"
		type2a_only = "False"
		type2b_only = "False"
		type2c_only = "False"
		rcc_only = "False"
		rcc_no_pheo = "False"
		pheo_only = "False"
		no_polycythemia = "False"
		rcc_any = "False"

		
		if (t1_count > 0) and (t2_count == 0) and (t2a_count == 0) and (t2b_count == 0) and (t2c_count == 0):
			type1_only = "True"
		if ((t1_count > 0) or (t2b_count >0)) and (t2a_count == 0) and (t2c_count == 0):
			type1_or_2b_only = "True"
		if (t2_count > 0) and (t1_count == 0):
			type2_only = "True"
		if (t1_count == 0) and (t2a_count > 0) and (t2b_count == 0) and (t2c_count == 0):
			type2a_only = "True"
		if (t1_count == 0) and (t2a_count == 0) and (t2b_count > 0) and (t2c_count == 0):
			type2b_only = "True"
		if (t1_count == 0) and (t2a_count == 0) and (t2b_count == 0) and (t2c_count > 0):
			type2c_only = "True"
		if (rcc_count > 0) and (t2_count == 0) and (pheo_count == 0) and (hb_count == 0) and (r_cyst_count == 0) and (p_cyst_count == 0) and (pgl_count == 0) and (ra_count == 0) and (polycyt_count == 0):
			rcc_only = "True"
		if (rcc_count > 0) and (t2_count == 0) and (pheo_count == 0):
			rcc_no_pheo = "True"
		if (rcc_count == 0) and (pheo_count > 0) and (hb_count == 0) and (r_cyst_count == 0) and (p_cyst_count == 0) and (pgl_count == 0) and (ra_count == 0) and (polycyt_count == 0):
			pheo_only = "True"
		if (polycyt_count == 0):
			no_polycythemia = "True"
		if (rcc_count > 0) or (t2b_count > 0) or (t1_count > 0):
			rcc_any = "True"


		var_bins = [type1_only,type1_or_2b_only,type2_only,type2a_only,type2b_only,type2c_only,rcc_only,rcc_no_pheo,pheo_only,no_polycythemia,rcc_any]
		variant_data.extend(var_quants)
		variant_data.extend(var_bins)

		if vhl_disease_count != 0:
			total_VHL_dis_vars+=1


print total_VHL_count, "total VHL Disease instances found."
print total_VHL_dis_vars, "total variants with at least 1 VHL Disease entry."

#0-7,12-end
new_header_line = "Codon\tVariant\tEffect\tType\tSurface\tDssp\tNeemo\tBluues"
header_ext1 = "VHLdb_transmission_simple,vhldb_entries,vhldb_syndrome_count,vhldb_t1_count,vhldb_t2_count,vhldb_t2a_count,vhldb_t2b_count,vhldb_t2c_count,vhldb_rcc_count,vhldb_pheo_count,vhldb_hb_count,vhldb_t_cyst_count,vhldb_p_cyst_count,vhldb_pgl_count,vhldb_ra_count,vhldb_polycyt_count,vhldb_germline_count,vhldb_somatic_count,vhldb_homozygous_count".split(',')
header_ext2 = "vhldb_t1_only,vhldb_t1_or_2b_only,vhldb_t2_only,vhldb_t2a_only,vhldb_t2b_only,vhldb_t2c_only,vhldb_rcc_only,vhldb_rcc_no_pheo,vhldb_pheo_only,vhldb_no_polycythemia,vhldb_rcc_any".split(',')
header_out1 = '\t'.join(header_ext1)+'\t'
header_out2 = '\t'.join(header_ext2)+'\n'
vhldb_header_out = new_header_line+'\t'+header_out1+header_out2
print vhldb_header_out


check_cHGVS_names = open('/camp/lab/findlayg/home/shared/projects/SGE/VHL/vhldb_cHGVS.txt', 'w')
with open('/camp/lab/findlayg/home/shared/projects/SGE/VHL/20220805_vhldb_processed.txt', 'w') as VHLdb_out_file:
	VHLdb_out_file.write(vhldb_header_out)
	for variant in variant_dict:
		if variant == '':
			continue
		else:
			cHGVS = variant_dict[variant][1]
			variant_out1= variant_dict[variant][0:8]
			variant_out2 = variant_dict[variant][12:]
			variant_out = variant_out1 + variant_out2
			if len(variant_out) != len(vhldb_header_out.split('\t')):
				print "Error, length doesn't match header string"
				print len(variant_out), len(vhldb_header_out.split('\t'))
				time.sleep(1)
			else:
				check_cHGVS_names.write(cHGVS+"\n")
				output_list_of_strings = [str(i) for i in variant_out]
				if len(output_list_of_strings) != 38:
					print "Error detected", len(output_list_of_strings)


				VHLdb_out_file.write('\t'.join(output_list_of_strings)+"\n")
				
check_cHGVS_names.close()


#
'''
Line 8 is "Disease"

Line 9 is "Transmission"

Line 10 is "Comments" - not important to do anything with this.

Line 11 is PubMed ID - simply aggregate this for reference

'''

###Editing here -- have merged the entries that match based on variant id (e.g. c.1A>T), grouping data from 8, 9, 10, 11 
'''
header_line += 'vhldb_entries\t'

sv_list = []
snv_list = []
snv_count = 0

with open('/camp/lab/findlayg/home/shared/projects/SGE/VHL/vhldb_alleles.txt', 'w') as VHLdb_alleles:
	for i in variant_dict:
		if len(i.split(' ')) >= 2:
			multi_variants+=1
			VHLdb_alleles.write(i+'\n')
			print i
		else:
			single_variants +=1
			sv_list.append(i)
			print i
	
	VHLdb_alleles.write("single variants -- not snv"+'\n')
	
	for sv in sv_list:
		if '>' not in sv:
			VHLdb_alleles.write(sv+'\n')
		else:
			snv_list.append(sv)
			snv_count +=1

	VHLdb_alleles.write("single variants -- snv"+'\n')
	for snv in snv_list:
		VHLdb_alleles.write(snv+'\n')
		print snv
		print variant_dict[snv]

#with open('/camp/lab/findlayg/home/shared/projects/SGE/VHL/vhldb_by_snv.txt', 'w') as VHLdb_by_snv):

	

print header_line
print entry_count, "entries parsed"
print unique_alleles, "unique alleles"
print single_variants, "single variant alleles"
print snv_count, "snv count"

unique_disease_entries = 0
with open('/camp/lab/findlayg/home/shared/projects/SGE/VHL/vhldb_diseases.txt', 'w') as VHLdb_diseases:
	for my_key in line8_entries:
		unique_disease_entries+=1
		VHLdb_diseases.write(my_key+":"+str(line8_entries[my_key])+"\n")
	print unique_disease_entries

#for my_key in line8_entries:
#	print my_key, " counts:", line8_entries[my_key]

'''



