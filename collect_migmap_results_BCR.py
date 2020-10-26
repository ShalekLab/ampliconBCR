##Script to consolidate output from corrected MIGMAP alignments into consensus CDR3s for cells, B Cells
#Sam Kazer, 5/24/15

##---NOTES----
#1. Run in the directory of the corrected aligned reads, ie the tab delimited corected aligned files


#!/usr/bin/python # 2.7.14

import os
import collections
import glob
import pandas as pd # 0.22
import numpy as np # 1.15.2

#get the current directory, and change to it
curdir = os.getcwd()

##--Set up dataframe to fill with consensus from each cell, also counter to act as the row index
df = pd.DataFrame()
#thresholds for calling consensus, ie if top alignment is 2x greater (in freq) than the next, then its the consensus; min count thresh
fthresh = 2
cthresh = 25
cols_to_keep = ['Sample', 'chain', 'v', 'd', 'j', 'cdr3aa', 'cdr3nt', 'freq', 'count']

##--Open each alignment tab delimited file and look for the consensus TRA and TRB
for alignment in glob.glob(curdir + "/*corr_align.txt"):
	print('Working on ' + alignment.replace(curdir+'/','').lstrip('/'))
	
	#Set booleans for having recorded TRA and TRB, and variables containing the TRA and TRB (2 potential TRAs)
	IGH = False
	IGK = False
	IGL = False
	IGH_cons = pd.DataFrame()
	IGK_cons = pd.DataFrame()
	IGL_cons = pd.DataFrame()
	
	#no need for with structure when using read_table, so just grab all of the data and put into dataframe
	table = pd.read_table(alignment)
	
	#FIRST CHECK TO SEE IF WE MEASURED IGH, IGK, IGL IN THIS CELL, if so set them TRUE and set no consensus
	if (table['v'].str.contains('IGHV').sum(numeric_only=False) == 0):
		IGH = True
		IGH_cons = pd.Series(['No IGH measured'], index=['v'])
	
	if (table['v'].str.contains('IGKV').sum(numeric_only=False) == 0):
		IGK = True
		IGK_cons = pd.Series(['No IGK measured'], index=['v'])
		
	if (table['v'].str.contains('IGLV').sum(numeric_only=False) == 0):
		IGL = True
		IGL_cons = pd.Series(['No IGL measured'], index=['v'])	
	
	#while loop to keep searching for consensus if both TRA and TRB haven't been found
	while (IGH == False or IGK == False or IGL == False):	
		#Now look for the IGH alignment (similar to code in collect_migmap_results_TCR.py)
		if (IGH == False):
			index_IGH = table['v'].str.contains('IGHV')
			table_IGH = table[index_IGH].reset_index(drop=True)
			table_IGH = table_IGH.replace('\?', '0', regex=True)
			
			groupby_cdr3aa = table_IGH[['freq','count']].groupby(table_IGH['cdr3aa'])
			groupby_cdr3aa_table = groupby_cdr3aa.sum().sort_values('count', ascending=False)
			unique_cdr3 = len(groupby_cdr3aa_table.index)
			
			##---Determine the consensus IGH gene
			#we have 2 or more cdr3s measured
			if (unique_cdr3 > 1):
				if (groupby_cdr3aa_table.ix[0,'count'] >= fthresh*groupby_cdr3aa_table.ix[1,'count'] and groupby_cdr3aa_table.ix[0,'count'] >= cthresh):
					consensus_cdr3 = groupby_cdr3aa_table.index.values[0]
					IGH_cons = pd.Series(table_IGH.ix[table_IGH[table_IGH['cdr3aa'].str.contains(consensus_cdr3) == True].index.values[0]])
					IGH_cons.ix['freq'] = groupby_cdr3aa_table.ix[0,'freq']
					IGH_cons.ix['count'] = groupby_cdr3aa_table.ix[0,'count']
					IGH = True
			
				else: 
					IGH_cons = pd.Series(['No Consensus'], index=['v'])
					IGH = True
					
			#we have 1 cdr3 measured		
			elif (unique_cdr3 == 1):
				if (groupby_cdr3aa_table.ix[0,'count'] >= cthresh):
					consensus_cdr3 = groupby_cdr3aa_table.index.values[0]
					IGH_cons = pd.Series(table_IGH.ix[table_IGH[table_IGH['cdr3aa'].str.contains(consensus_cdr3) == True].index.values[0]])
					IGH_cons.ix['freq'] = groupby_cdr3aa_table.ix[0,'freq']
					IGH_cons.ix['count'] = groupby_cdr3aa_table.ix[0,'count']
					IGH = True
			
				else: 
					IGH_cons = pd.Series(['No Consensus'], index=['v'])
					IGH = True
		
		#IGL Alignment
		if (IGL == False):
			index_IGL = table['v'].str.contains('IGLV')
			table_IGL = table[index_IGL].reset_index(drop=True)
			table_IGL = table_IGL.replace('\?', '0', regex=True)
			
			groupby_cdr3aa = table_IGL[['freq','count']].groupby(table_IGL['cdr3aa'])
			groupby_cdr3aa_table = groupby_cdr3aa.sum().sort_values('count', ascending=False)
			unique_cdr3 = len(groupby_cdr3aa_table.index)
			
			##---Determine the consensus IGL gene
			#we have 2 or more cdr3s measured
			if (unique_cdr3 > 1):
				if (groupby_cdr3aa_table.ix[0,'count'] >= fthresh*groupby_cdr3aa_table.ix[1,'count'] and groupby_cdr3aa_table.ix[0,'count'] >= cthresh):
					consensus_cdr3 = groupby_cdr3aa_table.index.values[0]
					IGL_cons = pd.Series(table_IGL.ix[table_IGL[table_IGL['cdr3aa'].str.contains(consensus_cdr3) == True].index.values[0]])
					IGL_cons.ix['freq'] = groupby_cdr3aa_table.ix[0,'freq']
					IGL_cons.ix['count'] = groupby_cdr3aa_table.ix[0,'count']
					IGL = True
			
				else: 
					IGL_cons = pd.Series(['No Consensus'], index=['v'])
					IGL = True
					
			#we have 1 cdr3 measured		
			elif (unique_cdr3 == 1):
				if (groupby_cdr3aa_table.ix[0,'count'] >= cthresh):
					consensus_cdr3 = groupby_cdr3aa_table.index.values[0]
					IGL_cons = pd.Series(table_IGL.ix[table_IGL[table_IGL['cdr3aa'].str.contains(consensus_cdr3) == True].index.values[0]])
					IGL_cons.ix['freq'] = groupby_cdr3aa_table.ix[0,'freq']
					IGL_cons.ix['count'] = groupby_cdr3aa_table.ix[0,'count']
					IGL = True
			
				else: 
					IGL_cons = pd.Series(['No Consensus'], index=['v'])
					IGL = True
		
		#IGK Alignment
		if (IGK == False):
			index_IGK = table['v'].str.contains('IGKV')
			table_IGK = table[index_IGK].reset_index(drop=True)
			table_IGK = table_IGK.replace('\?', '0', regex=True)
			
			groupby_cdr3aa = table_IGK[['freq','count']].groupby(table_IGK['cdr3aa'])
			groupby_cdr3aa_table = groupby_cdr3aa.sum().sort_values('count', ascending=False)
			unique_cdr3 = len(groupby_cdr3aa_table.index)
			
			##---Determine the consensus IGK gene
			#we have 2 or more cdr3s measured
			if (unique_cdr3 > 1):
				if (groupby_cdr3aa_table.ix[0,'count'] >= fthresh*groupby_cdr3aa_table.ix[1,'count'] and groupby_cdr3aa_table.ix[0,'count'] >= cthresh):
					consensus_cdr3 = groupby_cdr3aa_table.index.values[0]
					IGK_cons = pd.Series(table_IGK.ix[table_IGK[table_IGK['cdr3aa'].str.contains(consensus_cdr3) == True].index.values[0]])
					IGK_cons.ix['freq'] = groupby_cdr3aa_table.ix[0,'freq']
					IGK_cons.ix['count'] = groupby_cdr3aa_table.ix[0,'count']
					IGK = True
			
				else: 
					IGK_cons = pd.Series(['No Consensus'], index=['v'])
					IGK = True
					
			#we have 1 cdr3 measured		
			elif (unique_cdr3 == 1):
				if (groupby_cdr3aa_table.ix[0,'count'] >= cthresh):
					consensus_cdr3 = groupby_cdr3aa_table.index.values[0]
					IGK_cons = pd.Series(table_IGK.ix[table_IGK[table_IGK['cdr3aa'].str.contains(consensus_cdr3) == True].index.values[0]])
					IGK_cons.ix['freq'] = groupby_cdr3aa_table.ix[0,'freq']
					IGK_cons.ix['count'] = groupby_cdr3aa_table.ix[0,'count']
					IGK = True
			
				else: 
					IGK_cons = pd.Series(['No Consensus'], index=['v'])
					IGK = True

	##----Now that we have the consensus (or lack there of) IGH/IGL/IGK, add some final touches (ie sample and IGH/IGL/IGK) and then add to our df
	#Add the type of TCR to the top of the series for each
	IGH_cons = pd.Series(['IGH'], index=['chain']).append(IGH_cons)
	IGL_cons = pd.Series(['IGL'], index=['chain']).append(IGL_cons)
	IGK_cons = pd.Series(['IGK'], index=['chain']).append(IGK_cons)
	
	#Add the sample name to the beginning of the series for each
	sample_name = alignment.replace('_L001_corr_align.txt', '_comb').replace(curdir+'/','').lstrip('/')
	IGH_cons = pd.Series([sample_name], index=['Sample']).append(IGH_cons)
	IGL_cons = pd.Series([sample_name], index=['Sample']).append(IGL_cons)
	IGK_cons = pd.Series([sample_name], index=['Sample']).append(IGK_cons)
	
	#Append our final df with the consensus from this alignment
	df = df.append(IGH_cons, ignore_index=True)
	df = df.append(IGL_cons, ignore_index=True)
	df = df.append(IGK_cons, ignore_index=True)
	
##---Clean up the df, grab the columns we want, and write out to a tab delimiited txt
#subset only the columns we want to keep for output
df_clean = df[cols_to_keep]
	
#write out to a text delimited file
df_clean.to_csv(curdir + '/consensus_alignments.txt', sep='\t', index=False)
	
