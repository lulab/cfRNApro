import numpy as np
import pandas as pd
import argparse
import math
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from threading import Thread
import warnings
import time
warnings.filterwarnings('ignore')

def bias(x):
    L = len(x)
    if sum(x[(L//2):]) * 2 < sum(x):
        return -1
    else:
        return 1

def mine(x):
    x = x / sum(x)
    return (np.power(np.linalg.norm(x,3),3) - np.power(np.linalg.norm(x,2),4)) /  np.power(np.linalg.norm(x,2),4)

def maxKL(x):
	KL= 0
	L = len(x)
	x = x / sum(x)
	x = np.array([float(i) for i in x if i > 0])
	M = max(x)
	for i in range(len(x)):
		KL += x[i] * np.log(x[i] / M)
	return KL
def KL(x):
	KL= 0
	L = len(x)
	x = x / sum(x)
	x = np.array([float(i) for i in x if i > 0])
	for i in range(len(x)):
		KL += x[i] * np.log(x[i] / (1/L))
	return KL

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

def shannon_entropy(arg):
	'''
	calculate shannon's H = -sum(P*log(P)). arg is a list of float numbers. Note we used
	natural log here.
	'''
	lst_sum = sum(arg)
	entropy = 0.0
	for i in arg:
		entropy += (i/lst_sum) * math.log(i/lst_sum)
	if entropy == 0:
		return 0
	else:
		return -entropy

def process(Args,_ref,_bed):
	abundance = []
	kl = []
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
			if sum(_CDF) == 0:
				kl.append(0)
			else:
				kl.append(KL(_CDF))
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
			if sum(_CDF) == 0:
				kl.append(0)
			else:
				kl.append(KL(_CDF))
			for i in range(1,len(_CDF)): # calculate cdf
				_CDF[i] += _CDF[i-1]
		if _CDF[-1] != 0 and len(_CDF) != 0:
			abundance.append(_CDF[-1]/(len(_CDF)))
		else:
			abundance.append(0)

	return abundance,kl

def main():
	parser = argparse.ArgumentParser(description='generate kl.txt from cdf')
	parser.add_argument('--ref',  	type=str, default='/BioII/lulab_b/huangkeyun/hechang/fragmentome/Degradation/ref/long_RNA.gencode.bed', help='ref file')
	parser.add_argument('--bed',  	type=str, default='/BioII/lulab_b/taoyuhuan/GSE68086/bedgraph/SRR1982584.bedgraph', help='bed file')
	parser.add_argument('--out',  	type=str, default='/BioII/lulab_b/baiyilan/project/DegradationScore/TEP2015/kl/SRR1982584.kl.txt', help='output file')
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
	
	abundance,kl = [], []
	for i in chr_list:
		#print(i + ': ' + str(features[i].done()))
		abundance += features[i].result()[0]
		kl += features[i].result()[1]

	df = pd.DataFrame()
	df['#name'] = reference.index
	df['length'] = (reference['end'] - reference['start']).values
	df['abundance'] = abundance
	df['kl'] = kl
	df.set_index('#name',inplace = True)
	df.to_csv(Args.out,sep='\t')


if __name__ == '__main__':
    main()


 