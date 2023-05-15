#usage on computing cluster at the CRICK (CAMP):
#python /camp/lab/findlayg/home/users/findlag/bin/parse_cbioportal_snvs.py /camp/lab/findlayg/home/shared/projects/SGE/VHL/2022102_VHL_cBioPortal_SNVs.txt /camp/lab/findlayg/home/shared/projects/SGE/VHL/2022102_VHL_cBioPortal_SNVs_dedup.txt

#Script is simply to look at cBioPortal entries and aggregate all entries for the same variant into a single line to prevent repeat counting 

import sys
import os
import subprocess
import time
from datetime import datetime
startTime = datetime.now()

var_dict = {}

def header_index(header_list, category): #and category must be a string that defines a category
    var_dict_index = header_list.index(category)
    return var_dict_index

my_header_list = ''
line_count = 0
unique_pos_alts = []
unique_CTDs = []
unique_OCs = []
counted_samples = []
sample_mutation_dict = {}
master_ids_counted = []

with open(sys.argv[1],'r') as cb_file:
	for line in cb_file:   #not aggregating yet -- need to define all cancer tpyes CTDs you want to track, and then count those
		if line_count == 0:
			line_count+=1
			my_header_list = line.strip().split("\t")			
		else:
			line_data = line.strip().split("\t")
			pos_alt = line_data[header_index(my_header_list,'pos_alt')]
			#print str(pos)+alt

			if pos_alt not in unique_pos_alts:
				unique_pos_alts.append(pos_alt)
			else:
				pass

			CTD = line_data[header_index(my_header_list,'Cancer.Type.Detailed')]
			OC = line_data[header_index(my_header_list,'Oncotree.Code')]
			if CTD not in unique_CTDs:
				unique_CTDs.append(CTD)
			else:
				pass

			if OC not in unique_OCs:
				unique_OCs.append(OC)
			else: 
				pass

	print len(unique_CTDs), "unique CTDs detected."

line_count = 0
with open(sys.argv[1],'r') as cb_file:
	for line in cb_file:
		if line_count == 0:
			line_count+=1
			for cancer_type in unique_CTDs:
				my_header_list.append(cancer_type)
				print my_header_list
				print len(unique_CTDs), ". <<-- Added this many headers that are cancer types."
		else:
			line_count+=1
			print line_count
			line_data = line.strip().split("\t")
			#filter to ignore samples from the same patient:

			#This will be the same for multiple variants seen in the same tumour
			alt = line_data[header_index(my_header_list,'Var')]
			pos = line_data[header_index(my_header_list,'Start.Pos')]
			sample_id = line_data[header_index(my_header_list,'Sample.ID')]
			pos_alt = line_data[header_index(my_header_list,'pos_alt')]
			nospp = int(line_data[header_index(my_header_list,'Number.of.Samples.Per.Patient')])
			sample_type = line_data[header_index(my_header_list,'Sample.Type')]
			CTD = line_data[header_index(my_header_list,'Cancer.Type.Detailed')]
			OC = line_data[header_index(my_header_list,'Oncotree.Code')]
			patient_sex = line_data[header_index(my_header_list,'Sex')]
			patient_id = line_data[header_index(my_header_list,'Other.Patient.ID')]
			patient_age = line_data[header_index(my_header_list,'Diagnosis.Age')]
			#this will be the same for multiple samples drawn from the same tumour
			master_id = patient_id+'.'+patient_age+'.'+patient_sex+'.'+CTD+'.'+pos_alt
			print master_id
			#QC on format of data:
			if len(line_data) != 34:
				print "error re: line data: length", len(line_data)
			if (str(pos)+alt) != line_data[header_index(my_header_list,'pos_alt')]:
				print "error, mismatch in pos_alt"
				time.sleep(2)
			#print str(pos)+alt

			if sample_id in list(sample_mutation_dict.keys()):

				if pos_alt in sample_mutation_dict[sample_id]:
					pass
				else:
					sample_mutation_dict[sample_id].append(pos_alt)
			else:
				sample_mutation_dict[sample_id] = [pos_alt]

			if master_id in master_ids_counted:
				continue #this means no more data will be extracted from this sample because the unique case (patient-tumour-mutation) was already included
			else:
				master_ids_counted.append(master_id)
			
				if pos_alt not in var_dict:
					print "variant added"
					var_dict[pos_alt] = line_data
					#will add a counts entry for each type of cancer in the database (56 CTDs)
					for cancer_type in unique_CTDs:
						if CTD == cancer_type:
							var_dict[pos_alt].append(1)
						else:
							var_dict[pos_alt].append(0)
				elif pos_alt in var_dict:
					var_dict[pos_alt][header_index(my_header_list,CTD)]+=1
					var_dict[pos_alt][header_index(my_header_list,'Sample.ID')] = var_dict[pos_alt][header_index(my_header_list,'Sample.ID')]+';'+sample_id

	print "finished second loop."
	my_header_list.extend(['total_cBP_entries','total_cBP_rcc','total_cBP_other'])
	my_header_list.extend(['OncoKB','CIViC','MyCancerGenome','CancerHotspot', '3DHotspot','multiple_vars','all_vars_present','independent_samples'])
	for pos_alt in unique_pos_alts:
		total_var_cancer_count = 0
		total_var_rcc_count = 0
		total_var_cancer_other = 0
		for i in range(34,len(var_dict[pos_alt])):
			#define term or set of terms to mark all RCC cases (Renal seems to work for all RCC types):
			if ("Renal" in my_header_list[i]):
				total_var_cancer_count += var_dict[pos_alt][i]
				total_var_rcc_count += var_dict[pos_alt][i]
			elif ("Pheochromocytoma" in my_header_list[i]):
				pass
			else:
				total_var_cancer_count += var_dict[pos_alt][i]
				total_var_cancer_other += var_dict[pos_alt][i]
		var_dict[pos_alt].extend([total_var_cancer_count,total_var_rcc_count,total_var_cancer_other])
		var_annotations = var_dict[pos_alt][header_index(my_header_list,'Annotation')].strip('"').split(';')
		print var_annotations
		sample_id = var_dict[pos_alt][header_index(my_header_list,'Sample.ID')]
		multiple_vars = 'No'
		all_vars_presnet = ''
		independent_samples = 0
		sample_id_list = sample_id.split(";")
		#will only report on multiple variants per sample if present in the first sample ID listed.
		if len(sample_mutation_dict[sample_id_list[0]]) >1:
			multiple_vars = "Yes"
			all_vars_present = ";".join(sample_mutation_dict[sample_id_list[0]])
		else:
			multiple_vars = "No"
			all_vars_present = sample_mutation_dict[sample_id_list[0]]
		if ";" in var_dict[pos_alt][header_index(my_header_list,'Sample.ID')]:
			independent_samples = len(var_dict[pos_alt][header_index(my_header_list,'Sample.ID')].split(';'))
		else:
			independent_samples = 1
		var_dict[pos_alt].extend(var_annotations)
		var_dict[pos_alt].extend([multiple_vars,all_vars_present,independent_samples])

#write it all out
with open(sys.argv[2],'w') as cBP_out:
	header_final_length = len(my_header_list)
	cBP_out.write('\t'.join(my_header_list)+'\n')
	for pos_alt in var_dict:
		out_line_length = len(var_dict[pos_alt])

		if out_line_length != header_final_length:
			print "error:  outline number of entries does not equal the length of the header list.  Outline length:", out_line_length, "header_length:", header_final_length 

		out_line = '\t'.join(str(x) for x in var_dict[pos_alt])+'\n'
		cBP_out.write(out_line)

#"Sample.ID"
#,"End.Pos","Ref","Var"
#"External.Patient.ID"