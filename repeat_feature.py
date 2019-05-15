

#Libraries
import json 
import os 
from Bio import SeqIO
import urllib2 
import sys
import ast
import numpy as np  
import statistics
import math
import xlsxwriter
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, MeanShift,Birch
from sklearn.cluster import AgglomerativeClustering
from Bio import pairwise2
import pandas as pd
import subprocess

import time
start_time = time.time()

PATH='./data/features/'

### Construction of Features ###

#Function to calculate the length of the contigs
def length(contigs):
	#Storing the lengths in a list 
	contig_length = []
	for item in contigs:
		contig_length.append(len(contigs[item]))

	return contig_length 


#Function to calculate average contig length
def average_contig_length(contigs):
	#Total Contig Length 
	total_contig_length = length(contigs)
	total_contig_sum = sum(total_contig_length)
	average_contig_len = total_contig_sum / len(total_contig_length)

	return average_contig_len 



#Function to calculate number of reads 
def read_information(file_name):
	count = 0 
	total_read_length = 0
	for seq_record in SeqIO.parse(file_name,"fasta"):
		count = count + 1 
		total_read_length = total_read_length + len(seq_record.seq)

	average_read_length = total_read_length / count 

	return_list = []

	return_list.append(count)
	return_list.append(average_read_length)

	return return_list


#Function to calculate the coverage  @Parameters -- G: Actual Genome Size

def coverage(G,name):
	#Defining the core coverage value 
	core_coverage = (read_information(name)[0] * read_information(name)[1]) / G 

	return core_coverage 


def N50(contigs):
	#List containing the length of contigs
	contig_length_list = length(contigs)
	#Total Length of the contigs
	total_length = sum(contig_length_list)
	#Length to be inspected 
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


def L50(contigs):
	#N50 measurement for the Contigs
	N50_measurement = N50(contigs)
	#Length List for the contig
	contig_length_list = length(contigs)
	#Finding the Position of the N50 measurement
	index = contig_length_list.index(N50_measurement) + 1 

	return index





def N50_feature(contigs):
	N50_feature_set = []
	N50_value = N50(contigs)
	contig_length_list = length(contigs)
	for item in contig_length_list:
		N50_feature_set.append(float(item)/float(N50_value))


	return N50_feature_set 

def read_length_features(contigs):
	#List comprising of the lengths
	contig_length_list = length(contigs)
	#Maximum Length
	max_length = max(contig_length_list)
	#Minimum Length
	min_length = min(contig_length_list)
	#Median Length
	median_length = statistics.median(contig_length_list)

	length_features = {}
	length_features["max"] = max_length 
	length_features["min"] = min_length 
	length_features["median"] = median_length 
	length_features["average"] = float(sum(contig_length_list)) / len(contig_length_list)

	return length_features 


def individual_length(contig):
	if (len(contig) == 0):
		print 1

	return len(contig)



#Size of window for dotplot
window_size = 6

#Function to calculate average repetition lengths using dot plot
def dotplot_count(contig):
	contig_len = len(contig)

	#1-D Array to hold the dotplot, n space
	dotplot = [0]*contig_len

	#Currently not used for final result. Retained for future additions
	sumofLen = 0

	#Returned value
	count = 0
	
	#Iteratively updates dot plot value
	for i in range(contig_len):

		#j is used to start from the end
		for j in range(contig_len-1,-1,-1):	

			#Checks for the diagonal of the dotplot
			if(i == j):
				dotplot[i] = 0
				continue

			#Line being formed is terminated here (last element)
			if(j == contig_len-1):
				if(dotplot[j]>=window_size):
					sumofLen = sumofLen + dotplot[j]
					count = count + 1

			#Match in dotplot
			if(contig[i] == contig[j]):
				if(j > 0):
					dotplot[j] = dotplot[j-1] + 1
				else:
					dotplot[j] = 1

			#There is no match. 
			else:
				if(j > 0):
					#Check if the line upto the diagonally previous cell is longer than window
					if(dotplot[j-1]>=window_size):
						sumofLen = sumofLen + dotplot[j-1]
						count = count + 1

				dotplot[j] = 0

	return count


import sys
if __name__ == '__main__':
	#Reading the file  --- Temporarily not required
	# file_genomic = open(sys.argv[1],"r")

	#contig_file_name = ""
	#read_file_name = "" 

	#Storing the data for the contigs in a dictionary 
	contigs = {}

	counter = 0
	contig_list = []
	#Parsing the scaffolds generated to store in the dictionary
	for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
		contigs[seq_record.id] = seq_record.seq
		contig_list.append((seq_record.id, seq_record.seq))
		if counter < 2:
			print (seq_record.id)
			print len(seq_record.seq)
		counter +=1
	print "-------------------------"

	#Storage of the contigs
	# main_contig = []
	# check =0

	# #List Append
	# for item in contigs:
	# 	check+=1
	# 	if check <15:
	# 		print len(contigs[item])
	# 		print item
	# 	main_contig.append(contigs[item])





	feat_names = ['Nb_repeats', 'Avg_rp_exp', 'Avg_rp_period', 'Avg_rp_length', 'Avg_rp_err']

	tval =0
	feat_df = pd.DataFrame(columns=["ID"] + feat_names)


	for contig_item in contig_list:
		tval+=1
		contig = contig_item[1]
		#Initializing a temporary list for features
		temp = []

		temp.append(contig_item[0])

		temp_contig_file = open('mreps/tempfile.fa', 'w')
		temp_contig_file.write('>'+contig_item[0] + '\n')
		temp_contig_file.write(str(contig) )
		temp_contig_file.close()

		mreps_out = subprocess.check_output(['./mreps/mreps','-res', '1', '-exp' ,'3' ,'-fasta', 'mreps/tempfile.fa'])
		line_wise = mreps_out.split('\n')
		
		Avg_rp_exp =0.0
		Avg_rp_period =0.0
		Avg_rp_length =0.0
		Avg_rp_err =0.0

		try:
			required_out = line_wise[6:]

			nb_repeats = int(required_out[-1])

			meta_data = required_out[:-1]



			nr = len(meta_data)

			for j in meta_data:
				item = j.split()
				rp_len = int(item[4])
				rp_per = float(item[5][1:-1])
				rp_exp = float(item[6][1:-1])
				rp_err = float(item[7])
				
				Avg_rp_exp += rp_exp
				Avg_rp_period += rp_per
				Avg_rp_length += rp_len
				Avg_rp_err += rp_err

			Avg_rp_exp /= nr
			Avg_rp_period /= nr
			Avg_rp_length /= nr
			Avg_rp_err /= nr
		except:
			nb_repeats=0


		temp.append(nb_repeats)
		temp.append(Avg_rp_exp)
		temp.append(Avg_rp_period)
		temp.append(Avg_rp_length)
		temp.append(Avg_rp_err)
		if tval<5:
			print temp

		#Feature 8: Repetitions :
		#temp.append(dotplot_count(contig))

		temp_df = pd.DataFrame([temp], columns=["ID"] + feat_names)
		feat_df = feat_df.append(temp_df)



	# Writing csv file
	print feat_df.head()
	feat_df.to_csv(str(PATH+ sys.argv[1].split('.')[0].split('/')[-1] )+'_repeatfeatures.csv', index=None)

	#Final Feature Array -- Numpy array on which computations will be done
	# feature_array = np.array(main_feature)

	# nf = 4
	# pca = PCA(n_components=nf)
	# pca.fit(feature_array)
	# X_new = pca.transform(feature_array)

	# #print X_new
	# #print(pca.explained_variance_ratio_) 



	# #Agglomerative Clustering#
	# ag = AgglomerativeClustering(n_clusters = 2, affinity = 'euclidean').fit(feature_array)

	# k =  ag.labels_
	# #print feature_array
	# print k
	# k = k.tolist()


	# print k.count(1)

	# print k.count(0)
	# ### END OF AGGLOMERATIVE CLUSTERING

	# kmeans = KMeans(n_clusters=2).fit(feature_array)
	# print kmeans.labels_
	# print kmeans.labels_.tolist().count(0)
	# print kmeans.labels_.tolist().count(1)

	# #########################################



	# brc = Birch(threshold=0.5, branching_factor=50, n_clusters=2, compute_labels=True, copy=True).fit(feature_array)

	# print brc.labels_


	# ############### FOR ALIGNMENT SCORES #####################
	# alignments = pairwise2.align.globalxx("ACCGT","ACG")

	# scores = []
	# for item in alignments:
	# 	scores.append(item[2])

	# scores.append(3.1)
	# print max(scores)	

	##########################################################
	print("--- %s seconds ---" % (time.time() - start_time))
	print "###############"


	# full_genome = []

	# for seq_record in SeqIO.parse("full_genome.fa","fasta"):
	# 	full_genome.append(seq_record.seq) 


	# reference_genome = full_genome[0]
	# print len(reference_genome)


	# total_scores = []

	# for item in contigs:
	# 	alignments = pairwise2.align.globalxx(reference_genome,contigs[item])
	# 	scores = []
	# 	for i in alignments:
	# 		scores.append(i[2])

	# 	total_scores.append(max(scores))


	# # Usage of Machine Learning to generate scores: 

	# #Labels for the training set 
	# Y = np.array(total_scores)

	# X = feature_array




