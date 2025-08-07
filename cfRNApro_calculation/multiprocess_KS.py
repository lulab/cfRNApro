import numpy as np
import pandas as pd
import argparse
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from threading import Thread
import warnings
import time
warnings.filterwarnings('ignore')


def squareSum(n):
	n = n / 2
	return  (n*(n+1)*(2*n+1)/3)


def text_save(filename, data, geneName):
    file = open(filename,'w')
    for i in range(len(data)):
    	file.write(geneName[i]+'\t')
    	file.write('\t'.join([str(j) for j in data[i]]))
    	file.write('\n')
    file.close()
    print("save successfully!") 

def process(Args,_ref,_bed):
	abundance = []
	res = []
	pdf = []
	for i in range(_ref.shape[0]):
#		if i % 1000 == 0:
#			print(str(i) + '...')
		_CDF = np.zeros(_ref.iloc[i,2] - _ref.iloc[i,1],dtype = int) # initialize cdf array
		if _ref.iloc[i,3] == '+': # pos and neg strands are processed differently
			_start = _ref.iloc[i,1] # ref gene start
			_end = _ref.iloc[i,2] # ref gene end
			_valid = _bed.iloc[np.array(_bed['start'] <= _end) & np.array(_bed['end'] >= _start),:] # fragments fall in ref gene
			_valid['start'] = (_valid.iloc[:,1] - _start).clip(lower=0).fillna(0) # relative start position for CDF
			_valid['end'] = (_valid.iloc[:,2] - _start).clip(upper=_end - _start).fillna(_end - _start) # relative end position for CDF
			for i in range(_valid.shape[0]): # calculate pdf
				_CDF[max(0,_valid.iloc[i,1]-1):max(0,_valid.iloc[i,2]-1)] += _valid.iloc[i,3]
			pdf.append(_CDF)
			for i in range(len(_CDF)-1,0,-1): # calculate cdf
				_CDF[i-1] += _CDF[i]
			_CDF = _CDF[::-1] # reverse pos strand
			
		else:
			_start = _ref.iloc[i,1] # ref gene start
			_end = _ref.iloc[i,2] # ref gene end
			_valid = _bed.iloc[np.array(_bed['start'] <= _end) & np.array(_bed['end'] >= _start),:] # fragments fall in ref gene
			_valid['start'] = (_valid.iloc[:,1] - _start).clip(lower=0).fillna(0) # relative start position for CDF
			_valid['end'] = (_valid.iloc[:,2] - _start).clip(upper=_end - _start).fillna(_end - _start) # relative end position for CDF
			for i in range(_valid.shape[0]): # calculate pdf
				_CDF[_valid.iloc[i,1]:_valid.iloc[i,2]] += _valid.iloc[i,3]
			pdf.append(_CDF)
			for i in range(1,len(_CDF)): # calculate cdf
				_CDF[i] += _CDF[i-1]
		if _CDF[-1] != 0 and len(_CDF) != 0:
			abundance.append(_CDF[-1]/(len(_CDF)))
			#areaScore = (   sum(_CDF) - _CDF[-1] * (len(_CDF)-1)/2 ) / np.sqrt(   (np.sum(np.square(_CDF - np.mean(_CDF)))) + squareSum(len(_CDF)-1)   )
			# areaScore = ( sum(_CDF) - _CDF[-1] * (len(_CDF)-1)/2 ) / (_CDF[-1] * len(_CDF))
			# res.append(areaScore)
			pos = max(_CDF / _CDF[-1] - np.array(range(len(_CDF)))/len(_CDF))
			neg = min(_CDF / _CDF[-1] - np.array(range(len(_CDF)))/len(_CDF))
			if pos + neg >= 0:
				ksScore = pos
			else:
				ksScore = neg
			res.append(ksScore)

		else:
			abundance.append(0)
			res.append(0)

	return abundance,res#,pdf

def main():
	parser = argparse.ArgumentParser(description='generate ks.txt from cdf')
	parser.add_argument('--ref',  	type=str, default='/BioII/lulab_b/huangkeyun/hechang/fragmentome/Degradation/ref/long_RNA.gencode.bed', help='ref file')
	parser.add_argument('--bed',  	type=str, default='/BioII/lulab_b/baiyilan/project/DegradationScore/Pico/bedgraph/CRC-2124325.bedgraph', help='bed file')
	parser.add_argument('--out',  	type=str, default='/BioII/lulab_b/baiyilan/project/DegradationScore/Pico/ks/CRC-2124325.ks.txt', help='output file')
	Args = parser.parse_args()
	chr_list = ['chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrM','chrX','chrY']
	reference = pd.read_csv(Args.ref,sep='\t',header = None,names = ['chr','start','end','name','value','strand'],index_col=3)
	reference.drop('value',axis = 1,inplace = True)
	bedgraph = pd.read_csv(Args.bed,sep='\t',header = None,names = ['chr','start','end','value'])
	geneName = list(reference.index)

	bed_dict = {}
	for i in chr_list:
		bed_dict[i] = bedgraph[bedgraph['chr'] == i].reset_index(drop = True)
	
	ref_dict = {}
	for i in chr_list:
		ref_dict[i] = reference[reference['chr'] == i]


	pool = ProcessPoolExecutor(max_workers=25)
	features = {}
	for i in chr_list:
		features[i] = pool.submit(process, Args,ref_dict[i],bed_dict[i])
	pool.shutdown()
	
	abundance,res = [], []
	#pdf = []
	for i in chr_list:
		#print(i + ': ' + str(features[i].done()))
		abundance += features[i].result()[0]
		res += features[i].result()[1]
		#pdf += features[i].result()[2]

	df = pd.DataFrame()
	df['#name'] = reference.index
	df['length'] = (reference['end'] - reference['start']).values
	df['abundance'] = abundance
	df['max_dev'] = res
	df['ks_statistics'] = [0] * df.shape[0]
	df.set_index('#name',inplace = True)
	df.to_csv(Args.out,sep='\t')

	# file = open(Args.out,'w')
	# for i in range(len(pdf)):
	# 	file.write('\t'.join([str(j) for j in pdf[i]]))
	# 	file.write('\n')
	# file.close()

#	text_save('smallref.cdf',ans,geneName) # save results



if __name__ == '__main__':
    main()


