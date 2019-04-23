import sys
import logging
import itertools
import operator
import vcf
from cyvcf2 import VCF, Writer
from collections import defaultdict
import collections
from .customcontainer import DefaultOrderedDict
import gc
import math
from sys import getsizeof, stderr
import numpy as np
from difflib import SequenceMatcher
try:
	from reprlib import repr
except ImportError:
	pass
cimport cython
cimport cpp

def create_blocks(target_file, out_blockends_file, variant_set, all_variant_set):
	vcf_reader = vcf.VcfReader(target_file)._vcf_reader
	#vcf_reader = vcf.Reader(open(target_file, 'r'))
	i = 0
	blocks, ps = [], []
	intE, E, pairs = [], [], []
	count = 0
	phasedcount = 0
	haplo = {}
	
	for record in vcf_reader:
		for sample in vcf_reader.samples:
			#if record.genotype(sample).phased: #and record.POS in variant_list ? (von compute_ref vorher übergeben)
			if (record.POS in variant_set):		
				#phasing information for each variant position is stored in dict
				haplo[record.POS] = (record.genotype(sample)['GT'][0],record.genotype(sample)['GT'][2])
				phasedcount += 1
				#variants without a PS tag get an "unknown" value (this should be equal to the unphased ones)
			#	if not record.genotype(sample)['PS']:
				if not (record.genotype(sample).phased ):	
					count += 1	
					haplo[record.POS] = (record.genotype(sample)['GT'][0],record.genotype(sample)['GT'][2])
					#homozygous sites are set to the phase set of the last variant position that was phased				
					if (record.genotype(sample).gt_type == 2 and len(ps) != 0 and phasedcount >= 1):
						index = -1
						while not isinstance(ps[index], int):
							index -= 1
						if (index >= -len(ps)):
							ps.append(ps[index])	
				else:
					phase_set = record.genotype(sample)['PS']
					ps.append(phase_set)
####  homozygous positions are currently left out #######
#			if not (record.genotype(sample).phased ):
#				haplo[record.POS] = (record.genotype(sample)['GT'][0],record.genotype(sample)['GT'][2])
#				#homozygous sites are set to the phase set of the last variant position that was phased				
#				if (record.genotype(sample).gt_type == 2 and len(ps) != 0 and phasedcount >= 1):
#					index = -1
#					while not isinstance(ps[index], int):
#						index -= 1
#					if (index >= -len(ps)):
#						ps.append(ps[index])	
#				#unphased sites are set to "unknown"				
#				else:
#					count += 1
#					ps.append('U'+str(count))	

	vcf_reader = vcf.VcfReader(target_file)._vcf_reader
	count = 0
	positionslist = []
	for record in vcf_reader:
		for sample in vcf_reader.samples:
			if (record.POS in all_variant_set):
				count+=1
				positionslist.append(record.POS)
				haplo[record.POS] = (record.genotype(sample)['GT'][0],record.genotype(sample)['GT'][2])
	print("positionslist: ", positionslist[:50])
	print("count:", count)
	print("len all_variant set: ", len(all_variant_set))
	print("len haplo: ", len(haplo.keys()))
	print("not phased : ", count)
	ends = DefaultOrderedDict(list)
	starts = DefaultOrderedDict(list)
	#The end positions of each phase set are computed and stored as block ending positions
	for i in range(len(ps)):
		phase_set = ps[i]	
		if (i == 0):
			starts[phase_set].append(i)
		if (i > 0 and phase_set != blocks[i-1]):
			ends[ps[i-1]].append(i-1)
			starts[phase_set].append(i)
		if (i == len(ps)-1):
			ends[phase_set].append(i)
		blocks.append(phase_set)
	#pairs is a list that contains one list of tuples per phase set ("key")
	for key in ends.keys():
		if (len(ends[key]) == 1):
			pairs.append([(starts[key][0], ends[key][0])])
		else:
			pairlist = []
			for i in range(0,len(ends[key])):
				pairlist.append((starts[key][i], ends[key][i]))
			pairs.append(pairlist)
	#create list E containing block boundaries
	#add first block boundary (if part of an intermediate block, remove later)  
	E.append(pairs[0][-1][1])
	#internal block boundaries are found if present
	for i in range(1,len(pairs)):
		internal = False
		for j in range(0, i):
			if (pairs[i][-1][1] < pairs[j][-1][1]):
				internal = True
			#detect whether any items need to be removed since they are part of nested blocks
			if (len(pairs[i]) > 1):
				if (pairs[i][0][0] > pairs[j][0][0] and pairs[i][0][0] < pairs[j][-1][0] and pairs[i][-1][0] > pairs[j][-1][1]):
					if (pairs[j][-1][1] in E):
						E.remove(pairs[j][-1][1])
						intE.append((pairs[j]))	
		if internal:
			intE.append(pairs[i])
		else:
			E.append(pairs[i][-1][1])
	intE.sort()
	#add ending position of intermediate blocks (for nested blocks: of the last intermediate block) to the boundary set
	E2 = []
	for pairlist in intE:
		E2.append(pairlist[-1][1])
	E_whole = E+E2
	E_whole.sort()
	
#	E_whole = []
#	for e in range(len(variant_set)):
#		E_whole.append(e)
#	intE = []
#	E = []
#	for e in range(len(variant_set)):
#		E.append(e)
	
	E_file = open(out_blockends_file,'w')
	for e in E_whole:
		if (e == E_whole[len(E_whole)-1]):
			E_file.write(str(e))
		else:
			E_file.write(str(e))
			E_file.write(',')
			
	return(E_whole,haplo, intE)
	
#reads a VCF file and extracts the haplotype information, both haplotypes are written into a file as strings
def compute_haplotypes(newfile, resultfilename, variant_set):
	targetset = set()
	vcf_target = VCF(newfile)
	haplo1 = ""
	haplo2 = ""
	usedPos = []
	for v in vcf_target:
		if (v.POS in variant_set):
			targetset.add(v.POS)
			comp1 = 0
			comp2 = 0
			if (len(v.gt_bases[0].split('|')) > 1):
				leftbase = v.gt_bases[0].split('|')[0]
				rightbase = v.gt_bases[0].split('|')[1]
			elif(len(v.gt_bases[0].split('/')) > 1):
				leftbase = v.gt_bases[0].split('/')[0]
				rightbase = v.gt_bases[0].split('/')[1]
			if (leftbase == v.REF):
				comp1 = 0
			elif (leftbase == v.ALT[0]):
				comp1 = 1
			if (rightbase == v.REF):
				comp2 = 0
			elif (rightbase == v.ALT[0]):
				comp2 = 1
			haplo1 += str(comp1)	
			haplo2 += str(comp2)
			usedPos.append(v.POS)
	originalhaplofile = open(resultfilename,'w')
	originalhaplofile.write(haplo1)
	originalhaplofile.write('\n')	
	originalhaplofile.write(haplo2)
	return(haplo1,haplo2,usedPos)

#tests the nested blocks
def update_haplotypes(haplofile, E, pathfile, intE,unknown_variants, testlist1, testlist2, usedPos):
	#read the haplotypes from file and store them in strings
	H_A, H_B = "",""
	with open(haplofile) as filestream:
		for i,line in enumerate(filestream):
			if (i==0):
				H_A = line.strip('\n')
			if (i==1):
				H_B = line.strip('\n')

	#read the paths from file and store them in lists		
	helppath1 = []
	helppath2 = []
	with open(pathfile) as filestream:
		for i,line in enumerate(filestream):
			if (i==0):
				for num in line.strip().split(' '):				
					helppath1.append(num)
			elif (i==1):
				for num in line.strip().split(' '):
					helppath2.append(num)
	helppath1.append(helppath1[len(helppath1)-1])
	helppath2.append(helppath2[len(helppath2)-1])
	path1, path2 = [], []
	for i in range(0, len(helppath1)-1):
		path1.append(int(float(helppath1[i])))
	for i in range(0, len(helppath2)-1):
		path2.append(int(float(helppath2[i])))
	unique_refs = set(path1)
	print("unique refs: ", unique_refs)
	print("number of unique refs: ", len(unique_refs))
	ha = list(H_A)
	hb = list(H_B)

	test1, test2 = [], []
	map_old_to_new = [i for i in range(len(ha))]
	for e in E:
		test1.append(testlist1[e])
	for i in sorted(unknown_variants.keys()):
		ha.insert(i,unknown_variants[i][0])
		hb.insert(i,unknown_variants[i][1])
		usedPos.insert(i,unknown_variants[i][2])
		path1.insert(i, "nan")
		path2.insert(i, "nan")
		for e in range(len(E)):
			if E[e] >= i:
				E[e] += 1
		for el in range(len(map_old_to_new)):				
			if map_old_to_new[el] >= i:
				map_old_to_new[el] += 1
#	print(E)
	print("map old to new 6: ", map_old_to_new[6])
	for e in E:
		test2.append(testlist2[e])
	for i in range(len(test1)):
		if (test1[i] != test2[i]):
			print("not equal: ", i)
			
	unknown_is_blockend = 0
	for e in E:
		if e in unknown_variants.keys():
			unknown_is_blockend += 1
	print("unknown variant is at block end: ", unknown_is_blockend)
	
	H_A = "".join(str(el) for el in ha)
	H_B = "".join(str(el) for el in hb)
	#resolve any switches present in the paths and update the haplotypes accordingly
	(newhaplo1, newhaplo2, newpath1, newpath2, cutpositions) = improve_paths(list(path1),list(path2),E,intE, H_A,H_B, map_old_to_new)

#	(newhaplo1, newhaplo2, newpath1, newpath2, cutpositions) = improve_paths(list(path1),list(path2),E,intE, H_A,H_B)
	return((newhaplo1, newhaplo2), (newpath1, newpath2), usedPos, cutpositions, H_A)

#tests nested blocks
def improve_paths(path1, path2, E, intE, H_A, H_B, map_old_to_new):
#def improve_paths(path1, path2, E, intE, H_A, H_B):
	print("in improve paths, number of internal block ends: ", len(intE))
	newpath1 = path1[:]
	newpath2 = path2[:]
	#store original paths for later use to not lose this information when updating the paths
	oldpath1 = newpath1[:]
	oldpath2 = newpath2[:]
	
	newhaplo1 = list(H_A)[:]
	newhaplo2 = list(H_B)[:]
	oldhaplo1 = newhaplo1[:]
	oldhaplo2 = newhaplo2[:]
	switchcounter = 0
	before = "".join(newhaplo1)
#	print("newhaplo before: ", before)
	#at every block ending positions, it is checked whether the paths cross around that position
	cutpositions = []
	for i in E:
		if (i != len(path1)-1 and i < len(path1)):
			old_e = map_old_to_new.index(i)
			prec = map_old_to_new[old_e+1]
		#	print("path i, i+1: ", path1[i],path1[prec])
		#	prec =i+1
			#when a crossing between paths is found at position i, the paths are swapped from position 0 to i, the haplotypes accordingly
		#	if (path1[i] != path1[i+1] and path2[i] != path2[i+1] and (path1[i] == path2[i+1] or path2[i] == path1[i+1])):
			if (path1[i] != path1[prec] and path2[i] != path2[prec] and (path1[i] == path2[prec] or path2[i] == path1[prec])):
				cutpositions.append((i,prec))
				switchcounter += 1
				newpath1[0:i+1] = oldpath2[0:i+1]
				newpath2[0:i+1] = oldpath1[0:i+1]
				oldpath1 = newpath1[:]
				oldpath2 = newpath2[:]
				newhaplo1[0:i+1] = oldhaplo2[0:i+1]
				newhaplo2[0:i+1] = oldhaplo1[0:i+1]
#				newpath1[0:prec] = oldpath2[0:prec]
#				newpath2[0:prec] = oldpath1[0:prec]
#				oldpath1 = newpath1[:]
#				oldpath2 = newpath2[:]
#				newhaplo1[0:prec] = oldhaplo2[0:prec]
#				newhaplo2[0:prec] = oldhaplo1[0:prec]
				oldhaplo1 = newhaplo1[:]
				oldhaplo2 = newhaplo2[:]
	after = "".join(newhaplo1)
#	print("newhaplo after: ", after)
	print("switch happened in %s times: ", switchcounter)
	for pairlist in intE:
		switch = check_switching(pairlist,newpath1, newpath2)
		if switch:
			for pair in pairlist:				
				i = pair[0]
				j = pair[1]
				newpath1[i:j+1] = oldpath2[i:j+1]
				newpath2[i:j+1] = oldpath1[i:j+1]
				oldpath1 = newpath1[:]
				oldpath2 = newpath2[:]
				newhaplo1[i:j+1] = oldhaplo2[i:j+1]
				newhaplo2[i:j+1] = oldhaplo1[i:j+1]
				oldhaplo1 = newhaplo1[:]
				oldhaplo2 = newhaplo2[:]
	return("".join(newhaplo1),"".join(newhaplo2), newpath1, newpath2, cutpositions)

#tests nested blocks
def check_switching(pairlist, path1,path2):
	switch = 0
	to_switch = False
	for pair in pairlist:
		i = pair[0] 
		j = pair[1]
		if (path1[j] != path1[j+1] and path2[j] != path2[j+1] and (path1[j] == path2[j+1] or path2[j]==path1[j+1])):
			switch += 1
		if (path1[j] == path1[j+1] or path2[j]==path2[j+1]):
			switch -= 1
		if (i != 0):
			if (path1[i] != path1[i-1] and path2[i] != path2[i-1] and (path1[i] == path2[i-1] or path2[i]==path1[i-1])):
				switch += 1
			if (path1[i] == path1[i-1] or path2[i] == path2[i-1]):
				switch -= 1
	if (switch > 0):
		to_switch = True
	return(to_switch)
	
def compute_referencepanel(ref_file, target_file, variant_set, ref_output, samples_to_use, variants_to_use):
	#compute size of ref_matrix and initialize it with zeros
	#width=number of variants between first and last relevant variant
	#height=twice the number of samples, as each sample offers two haplotypes

	height = 2*len(samples_to_use)

	strlist = []
	string = ''
	for i in range(0,height):
		strlist.append(string)
	vcf_ref = VCF(ref_file, lazy=True)
	print(vcf_ref.samples[:10])
	doubleset = set()
	index = -1
	used_sample_names = []
	for variant in variants_to_use:
		if (variant.POS not in doubleset):
			index += 1
			doubleset.add(variant.POS)
			counter = 0
			for baseindex in samples_to_use:
				
				genotype = variant.genotypes[baseindex]
				strlist[counter] += str(genotype[0])		
				strlist[counter+1] += str(genotype[1])
				counter += 2	
	sampleindex = 0
	sample_map = {}
	for baseindex in samples_to_use:	
		used_sample_names.append(vcf_ref.samples[baseindex])	
		sample_map[sampleindex] = vcf_ref.samples[baseindex]
		sampleindex += 1
	print("number of used samples: ", len(used_sample_names))
	print("used sample names: ", used_sample_names)
	output = open(ref_output, 'w')	
	for string in strlist:
		output.write(string+'\n')
	output.close()
	return (output)
	
def compute_referencepanel_heuristic(ref_file, target_file, variants_to_use, ref_output, start, end):
	height=(end-start)*2

	strlist = []
	string = ''
	for i in range(0,height): 
		strlist.append(string)
	vcf_ref = VCF(ref_file, lazy=True)
	doubleset = set()
	
	for variant in variants_to_use:
		if (variant.POS not in doubleset):
			doubleset.add(variant.POS)
			counter = 0
			for genotype in variant.genotypes[start:end]:			
			#	if (base[0] == variant.REF):
			#		strlist[counter]+="0"
			#	else:		
			#		strlist[counter]+="1"
			#	if (base[(2)] == variant.REF):
			#		strlist[counter+1]+="0"
			#	else:
			#		strlist[counter+1]+="1"
				strlist[counter] += str(genotype[0])			
				strlist[counter+1] += str(genotype[1])			
				counter += 2
	output = open(ref_output, 'w')
	#check whether strlist has length 1000 (or 2504?)
	for string in strlist:
		output.write(string+'\n')
	#compute hamming distances and use the indices of best values
	output.close()
	return (output)

def hamming(s0, s1):
	assert len(s0) == len(s1)
	return sum(c0 != c1 for c0, c1 in zip(s0, s1))

def hamming_homozyg(s0,s1):
	counter = 0
	for c0,c1 in zip(s0,s1):
		if (int(c0) ==1):
			counter +=1
	print("count: ", counter)
	errors = sum(c0 != c1 for c0,c1 in zip(s0,s1) if int(c0)==1)
	return(errors/counter)

def lcs(X, Y): 
	# find the length of the strings 
	m = len(X) 
	n = len(Y) 

	# declaring the array for storing the dp values 
	L = [[None]*(n + 1) for i in xrange(m + 1)] 

	"""Following steps build L[m + 1][n + 1] in bottom up fashion 
	Note: L[i][j] contains length of LCS of X[0..i-1] 
	and Y[0..j-1]"""
	for i in range(m + 1): 
		for j in range(n + 1): 
			if i == 0 or j == 0 : 
				L[i][j] = 0
			elif X[i-1] == Y[j-1]: 
				L[i][j] = L[i-1][j-1]+1
			else: 
			#	L[i][j] = max(L[i-1][j], L[i][j-1]) 
				L[i][j] = L[i-1][j-1]
				
	# L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1] 
	return L[m][n] 


def compute_referencepanel_heuristic2(ref_file, target_file, variants_to_use, ref_output, hap1,hap2):
	vcf_ref = VCF(ref_file, lazy=True)
	height = 2*len(vcf_ref.samples)

	strlist = []
	genlist = []
	string = ''
	for i in range(0,height):
		strlist.append(string)
	for i in range(0, len(vcf_ref.samples)):	
		genlist.append(string)
		
	doubleset = set()
	index = -1
	for variant in variants_to_use:
		if (variant.POS not in doubleset):
			index += 1
			doubleset.add(variant.POS)
			counter = 0
			for genotype in variant.genotypes:						
			#	genlist[counter] += str(genotype[0]+genotype[1])
			#	counter += 1
				strlist[counter] += str(genotype[0])	
				strlist[counter+1] += str(genotype[1])		
				counter += 2

	print("length stringlist: ", len(strlist))
	print("length genlist: ", len(genlist))
	print("length genotype: ", len(genlist[0]))
	#test1: hamming distance der switch sequences
	#test2: hamming distance der haplos direkt (minimum aus HD zu hap1 und hap2)
	switch_hap1 = ''.join(('0' if hap1[i-1] == hap1[i] else '1') for i in range(1, len(hap1)))
	target_gen = ''.join(str(int(hap1[i])+int(hap2[i])) for i in range(0, len(hap1)))
#	print("length target genotype: ", len(target_gen))
#	print("target gen 0: ", target_gen[0])
	hamming_values = {}
	for i in range(len(strlist)):
		string = strlist[i]
		assert(len(string)==len(hap1))
		switch_string = ''.join(('0' if string[i-1] == string[i] else '1') for i in range(1, len(string)))
		value = hamming(switch_string, switch_hap1)
	#	hd = value/len(switch_string)
	#	hamming_values[i] = hd
		hamming_values[i] = value
	hd_values_sorted = sorted(hamming_values.items(), key=operator.itemgetter(1))
	print("hamming distances: ", hd_values_sorted)
	samples_to_use = [i[0] for i in hd_values_sorted[:600]]
	print("samples to use: ", samples_to_use)
#	samples_to_use = sorted(list(set([i[0]//2 for i in hd_values_sorted[:500]])))
	samples_to_use = [i[0]//2 for i in hd_values_sorted[:600]]	
	print("samples to use: ", samples_to_use)
	print("number of samples: ", len(samples_to_use))

#	for i in range(len(genlist)):
#		string = genlist[i]
#		assert(len(string) == len(target_gen))
#		value = hamming_homozyg(string, target_gen)
#		hamming_values[i] = value
#	hd_values_sorted = sorted(hamming_values.items(), key=operator.itemgetter(1))
#	print("hamming distances: ", hd_values_sorted)
#	samples_to_use = [i[0] for i in hd_values_sorted[:500]]
#	print("samples to use: ", samples_to_use)
#	samples_to_use = [i[0]//2 for i in hd_values_sorted[:500]]
#	print("samples to use: ", samples_to_use)
#	print("number of samples: ", len(samples_to_use))
	
	return (samples_to_use)
