
import json 
import os 
#import biopython
#import urllib2 
import sys
import ast
#import numpy as np  
import statistics
import math
#import xlsxwriter
#from sklearn.decomposition import PCA
#from sklearn.cluster import KMeans, MeanShift
#from sklearn.cluster import AgglomerativeClustering
#from Bio import pairwise2
import pandas as pd
import subprocess


import time
start_time = time.time()

PATH='./data/features/'



#Contig length
def length(contigs):
	
	contig_length = []  
	for item in contigs:
		contig_length.append(len(contigs[item]))

	return contig_length 


#Average contig length
def average_contig_length(contigs):
	
	total_contig_length = length(contigs)
	total_contig_sum = sum(total_contig_length)
	average_contig_len = total_contig_sum / len(total_contig_length)

	return average_contig_len 



# Number of reads & average read length 
def read_information(file_name):
	count = 0 
	total_read_length = 0
	fo = open(file_name, "r")
	for seq_record in fo.readlines():
		count = count + 1 
		total_read_length = total_read_length + len(seq_record.seq)

	average_read_length = total_read_length / count 

	return_list = []

	return_list.append(count)
	return_list.append(average_read_length)

	return return_list






def N50(contigs):
	contig_length_list = length(contigs)
	total_length = sum(contig_length_list)
	inspect_length = total_length * 0.5 
	#Sorted list of contigs
	sorted_contig_length_list = sorted(contig_length_list,reverse=True)

	temp_sum = 0 
	#Loop to find the least item among the contig set which makes upto 50% of the contig lengths
	for item in sorted_contig_length_list:
		temp_sum += item 
		if temp_sum >= inspect_length:
			answer = item 
			break 
	return answer 




def N90(contigs):
	contig_length_list = length(contigs)
	total_length = sum(contig_length_list)
	inspect_length = total_length * 0.9 
	sorted_contig_length_list = sorted(contig_length_list,reverse=True)

	temp_sum = 0 
	for item in sorted_contig_length_list:
		temp_sum += item 
		if temp_sum >= inspect_length:
			answer = item 
			break 


	return answer 




def L50(contigs):

	N50_measurement = N50(contigs)
	contig_length_list = length(contigs)
	index = contig_length_list.index(N50_measurement) + 1 

	return index





#Coverage ( G: Actual Genome Size )
def coverage(G,name):
	#Coverage = no.of reads*av. read length / original genome length
	core_coverage = (read_information(name)[0] * read_information(name)[1]) / G 

	return core_coverage 




def GC_content(contigs):
	GC_count =0
	for item in contigs:
		if 'G' in item: 
			GC_count+=1
		if 'C' in item: 
			GC_count+= 1
	GC_count = float((GC_count/total)*100)
	return GC_count
	
	



def N50_feature(contigs):
	N50_feature_set = []
	N50_value = N50(contigs)
	contig_length_list = length(contigs)
	for item in contig_length_list:
		N50_feature_set.append(float(item)/float(N50_value))


	return N50_feature_set 




def read_length_features(contigs):
	
	contig_length_list = length(contigs)
	max_length = max(contig_length_list)
	min_length = min(contig_length_list)
	median_length = statistics.median(contig_length_list)#Median Length
	length_features = {}
	length_features["max"] = max_length 
	length_features["min"] = min_length 
	length_features["median"] = median_length 
	length_features["average"] = float(sum(contig_length_list)) / len(contig_length_list)

	return length_features 




def individual_length(contig):
	if (len(contig) == 0):
		print (1)

	return len(contig)




if __name__ == '__main__':


	#Storing the data for the contigs in a dictionary 
	contigs = {}

	counter = 0
	contig_list = []
	#Parsing the scaffolds generated to store in the dictionary
	for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
		contigs[seq_record.id] = seq_record.seq
		contig_list.append((seq_record.id, seq_record.seq))
		if counter < 4:
			print (seq_record.id)
			print (len(seq_record.seq))
		counter +=1
	print ("-------------------------")



	#Core Features : 
	contig_core_features = read_length_features(contigs)



	feat_names = ['Core_length', 'Dev_from_max', 'Dev_from_min', 'Dev_from_avg']
	feat_names = feat_names+ ['Dev_from_med', 'Dev_from_N50', 'Dev_from_L50']
	# feat_names = feat_names + ['Nb_repeats', 'Avg_rp_exp', 'Avg_rp_period', 'Avg_rp_length', 'Avg_rp_err']

	feat_df = pd.DataFrame(columns=["ID"] + feat_names)

	tval=0
	for contig_item in contig_list:
		tval+=1

		contig = contig_item[1]
		#Initializing a temporary list for features
		temp = []

		temp.append(contig_item[0])

		#Feature 1  : Core Length
		temp.append(individual_length(contig))

		#Feature 2 : Deviation from max
		temp.append(contig_core_features['max'] - individual_length(contig))


		#Feature 3 : Deviation from min
		temp.append(individual_length(contig) - contig_core_features['min'])


		#Feature 4 : Deviation from average
		avg = average_contig_length(contigs)
		dev_from_avg = individual_length(contig) - avg 
		if dev_from_avg>=0:
			temp.append(dev_from_avg)
		else:
			temp.append(-1*dev_from_avg)

		#Feature 5: Deviation from median
		dev_from_med = individual_length(contig) - contig_core_features['median']
		if dev_from_med>=0:
			temp.append(dev_from_med)
		else:
			temp.append(-1*dev_from_med)

		#Feature 6 : Deviation from N50
		n50 = N50(contigs)
		temp.append(float(individual_length(contig))/float(n50))


		#Feature 7: Deviation from L50
		l50 = L50(contigs)
		temp.append(float(individual_length(contig))/float(l50))


		#### Repetition Features below ####

		# temp_contig_file = open('mreps/tempfile.fa', 'w')
		# temp_contig_file.write('>'+contig_item[0] + '\n')
		# temp_contig_file.write(str(contig) )
		# temp_contig_file.close()

		# mreps_out = subprocess.check_output(['./mreps/mreps','-res', '1', '-exp' ,'3' ,'-fasta', 'mreps/tempfile.fa'])
		# line_wise = mreps_out.split('\n')
		
		# Avg_rp_exp =0.0
		# Avg_rp_period =0.0
		# Avg_rp_length =0.0
		# Avg_rp_err =0.0

		# try:
		# 	required_out = line_wise[6:]
		# 	nb_repeats = int(required_out[-1])
		# 	meta_data = required_out[:-1]
		# 	nr = len(meta_data)

		# 	for j in meta_data:
		# 		item = j.split()
		# 		rp_len = int(item[4])
		# 		rp_per = float(item[5][1:-1])
		# 		rp_exp = float(item[6][1:-1])
		# 		rp_err = float(item[7])
				
		# 		Avg_rp_exp += rp_exp
		# 		Avg_rp_period += rp_per
		# 		Avg_rp_length += rp_len
		# 		Avg_rp_err += rp_err

		# 	Avg_rp_exp /= nr
		# 	Avg_rp_period /= nr
		# 	Avg_rp_length /= nr
		# 	Avg_rp_err /= nr
		# except:
		# 	nb_repeats=0


		# temp.append(nb_repeats)
		# temp.append(Avg_rp_exp)
		# temp.append(Avg_rp_period)
		# temp.append(Avg_rp_length)
		# temp.append(Avg_rp_err)
		if tval<5:
			print (temp)

		if tval%400 == 0:
			print (tval)
			print("--- %s seconds ---" % (time.time() - start_time))



		temp_df = pd.DataFrame([temp], columns=["ID"] + feat_names)
		feat_df = feat_df.append(temp_df)

		if tval == 5000:
			break


	# Writing csv file
	print (feat_df.head())

	feat_df.to_csv(str(PATH+ sys.argv[1].split('.')[0].split('/')[-1] )+'_all_features.csv', index=None)



	##########################################################
	print("--- %s seconds ---" % (time.time() - start_time))
	print ("###############")







