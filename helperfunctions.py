# -*- coding: utf-8 -*-
from __future__ import print_function
import sys
import logging
import itertools
from itertools import chain
import vcf
from cyvcf2 import VCF, Writer
from collections import defaultdict
import collections
from collections import deque
from operator import itemgetter
import customcontainer
import gc
import math
from sys import getsizeof, stderr
import numpy as np
try:
    from reprlib import repr
except ImportError:
    pass


#Reads the phase set information from the VCF file and creates a haplotype block structure by outputting block ending positions. Haplotype positions and their phasing information is stored in a dictionary	
def create_blocks(target_file, out_blockends_file):
	vcf_reader=vcf.Reader(open(target_file, 'r'))
	i = 0
	blocks, ps = [], []
	intE, E, pairs = [], [], []
	count = 0
	phasedcount = 0
	haplo = {}
	for record in vcf_reader:
		for sample in vcf_reader.samples:
			if record.genotype(sample).phased:
				#phasing information for each variant position is stored in dict
				haplo[record.POS] = (record.genotype(sample)['GT'][0],record.genotype(sample)['GT'][2])
				phasedcount += 1
				#variants without a PS tag get an "unknown" value 
				if not record.genotype(sample)['PS']:
					count += 1
					ps.append('U'+str(count))
				else:
					phase_set = record.genotype(sample)['PS']
					ps.append(phase_set)	
			if not (record.genotype(sample).phased ):
				haplo[record.POS] = (record.genotype(sample)['GT'][0],record.genotype(sample)['GT'][2])
				#homozygous sites are set to the phase set of the last variant position that was phased				
				if (record.genotype(sample).gt_type == 2 and len(ps) != 0 and phasedcount >= 1):
					index = -1
					while not isinstance(ps[index], int):
						index -= 1
					if (index >= -len(ps)):
						ps.append(ps[index])	
				#unphased sites are set to "unknown"				
				else:
					count += 1
					ps.append('U'+str(count))		
	ends = customcontainer.DefaultOrderedDict(list)
	starts = customcontainer.DefaultOrderedDict(list)
	#The end positions of each phase set are computed and stored as block ending positions
	for i in range(len(ps)):
		phase_set = ps[i]	
		if (i == len(ps)-1):
			ends[phase_set].append(i)
		if (i == 0):
			starts[ps[i]].append(i)
		if (i > 0 and phase_set != blocks[i-1]):
			ends[ps[i-1]].append(i-1)
			starts[ps[i]].append(i)
		blocks.append(phase_set)
	for key in ends.keys():
		if (len(ends[key]) == 1):
			pairs.append((starts[key][0], ends[key][0]))
		else:
			pairs.append((starts[key][0], ends[key][len(ends[key])-1]))
	E.append(pairs[0][1])
	#internal block boundaries are found if present
	for i in range(1,len(pairs)):
		internal = False
		for j in range(0, i):
			if (pairs[i][1] < pairs[j][1]):
				internal = True			
		if internal:
			intE.append(pairs[i])
		else:
			E.append(pairs[i][1])

#use if the internal block boundaries are needed as well	
	E2 = []
	for pair in intE:
		E2.append(pair[0])
	E_whole = E+E2
	E_whole.sort()
	E_file = open(out_blockends_file,'w')
	for e in E_whole:
		if (e == E_whole[len(E_whole)-1]):
			E_file.write(str(e))
		else:
			E_file.write(str(e))
			E_file.write(',')
	return(E_whole, haplo, intE)	

#for later evaluation: Extract strings from the given haplotype block sequence and from the original haplotype file and store them in haplotypes.txt	
def compute_original(originalfile, newfile, resultfilename):
	targetlist = []
	vcf_target = VCF(newfile)
	haplo1 = ""
	haplo2 = ""
	#write haplotype information from the target VCF file into strings
	for v in vcf_target:
		targetlist.append(v.POS)
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
	
	#write haplotype information from the original VCF file into strings. Only positions also present in the target are used.
	vcf = VCF(originalfile)
	counter = -1
	errors = 0
	enderrors = 0
	errors1 = []
	originalhaplo1 = ""
	originalhaplo2 = ""
	for v in vcf:
		comp1 = 0
		comp2 = 0
		if v.POS in targetlist:
			counter += 1
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
			originalhaplo1 += str(comp1)	
			originalhaplo2 += str(comp2)
	#store all created strings in output file
	originalhaplofile = open(resultfilename,'w')
	originalhaplofile.write(haplo1)
	originalhaplofile.write('\n')	
	originalhaplofile.write(haplo2)
	originalhaplofile.write('\n')
	originalhaplofile.write(originalhaplo1)	
	originalhaplofile.write('\n')	
	originalhaplofile.write(originalhaplo2)
	return((haplo1,haplo2),(originalhaplo1,originalhaplo2))

#Used for blockparser_inference: Creates haplotype strings where each position that is missing (in comparison to an original haplotype file) gets an unknown haplotype (1/1 or 0/0 by default)
def compute_incomplete_haplos(originalfile, newfile, E_file, haplofilename):
	vcf = VCF(newfile)
	#store the haplotype information from the target VCF file inside a dictionary mapping the variants' positions to the haplotypes
	variantlist = {}
	for v in vcf:
		comp1, comp2 = 0,0
		if (len(v.gt_bases[0].split('|')) > 1):
			leftbase = v.gt_bases[0].split('|')[0]
			rightbase = v.gt_bases[0].split('|')[1]
		elif(len(v.gt_bases[0].split('/')) > 1):
			leftbase = v.gt_bases[0].split('/')[0]
			rightbase = v.gt_bases[0].split('/')[1]
		if (leftbase != '.' and rightbase != '.'):
			if (leftbase == v.REF):
				comp1 = 0
			elif (leftbase == v.ALT[0]):
				comp1 = 1
			if (rightbase == v.REF):
				comp2 = 0
			elif (rightbase == v.ALT[0]):
				comp2 = 1
		variantlist[v.POS] = (comp1, comp2)
	#create haplotype strings, where for each position that is in the original, but missing in the target, an unknown genotype is inserted
	orvcf = VCF(originalfile)
	originallist, intE = [],[]
	haplo1, haplo2 = "",""
	index = -1
	for v in orvcf:
		index += 1
		originallist.append(v.POS)
		if (v.POS in variantlist.keys()):
			haplo1 += str(variantlist[v.POS][0])	
			haplo2 += str(variantlist[v.POS][1])
		else:
			if (len(v.gt_bases[0].split('/')) > 1):
				haplo1 += str(1)
				haplo2 += str(1)
			else:
				haplo1 += str(0)
				haplo2 += str(1)
			#missing positions are treated as internal blocks
			intE.append(index)
	haplofile = open(haplofilename,'w')
	haplofile.write(haplo1)
	haplofile.write('\n')	
	haplofile.write(haplo2)
	e_file = open(E_file,'w')
	for e in intE:
		if (e == intE[len(intE)-1]):
			e_file.write(str(e))
		else:
			e_file.write(str(e)+',')
	return(haplo1, haplo2, intE)

#reads a VCF file and extracts the haplotype information, both haplotypes are written into a file as strings
def compute_haplos(newfile, resultfilename):
	targetlist = []
	vcf_target = VCF(newfile)
	haplo1 = ""
	haplo2 = ""
	for v in vcf_target:
		targetlist.append(v.POS)
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
	originalhaplofile = open(resultfilename,'w')
	originalhaplofile.write(haplo1)
	originalhaplofile.write('\n')	
	originalhaplofile.write(haplo2)
	return((haplo1,haplo2))

#crossovers within the path pair are found and resolved, while the haplotypes are updated accordingly whenever path parts are swapped
def improve_paths(path1, path2, E, H_A, H_B):
	newpath1 = path1[:]
	newpath2 = path2[:]
	#store original paths for later use to not lose this information when updating the paths
	oldpath1 = newpath1[:]
	oldpath2 = newpath2[:]
	
	#store haplotypes in strings and copy into lists to store for later use
	ha_s = str(H_A)
	hb_s = str(H_B)
	oldhaplo1 = ha_s[:]
	oldhaplo2 = hb_s[:]
	twosidedswitches = 0
	onesidedswitches = 0
	noswitch = 0
	samerow = 0
	bothrowssame = 0
	onesidedswitch_list = []
	bothrowssamelist = []
	twosidedswitchlist = []
	#at every block ending positions, it is checked whether the paths cross around that position
	for i in E:
		if (i != len(path1)-1 and i < len(path1)):
			#when a crossing between paths is found at position i, the paths are swapped from position 0 to i, the haplotypes accordingly
			if (path1[i] != path1[i+1] and path2[i] != path2[i+1] and (path1[i] == path2[i+1] or path2[i] == path1[i+1])):
				if (path1[i] == path2[i+1] and path2[i] == path1[i+1]):
					twosidedswitches += 1
					twosidedswitchlist.append(i)
				elif ((path1[i] == path2[i+1] and path2[i] != path1[i+1]) or (path1[i] != path2[i+1] and path2[i] == path1[i+1])):
					onesidedswitches += 1
					onesidedswitch_list.append(i)
				newpath1[0:i+1] = oldpath2[0:i+1]
				newpath2[0:i+1] = oldpath1[0:i+1]
				ha_s = oldhaplo2[0:i+1] + oldhaplo1[i+1:]
				hb_s = oldhaplo1[0:i+1] + oldhaplo2[i+1:]
				oldpath1 = newpath1[:]
				oldpath2 = newpath2[:]
				oldhaplo1 = ha_s[:]
				oldhaplo2 = hb_s[:]
			#crossing is incomplete; for testing purposes only
			elif ((path1[i] == path1[i+1] and path2[i]!=path2[i+1]) or (path1[i] != path1[i+1] and path2[i]==path2[i+1])):
				samerow += 1
				onesidedswitch_list.append(i)
			#no crossing occurs; for testing purposes only
			elif (path1[i] == path1[i+1] and path2[i]==path2[i+1]):
				bothrowssame += 1	
				bothrowssamelist.append(i)
			else:
				noswitch+= 1
	return((newpath1,newpath2), (ha_s, hb_s), twosidedswitches, onesidedswitches, noswitch, samerow, bothrowssame, onesidedswitch_list, bothrowssamelist, twosidedswitchlist)

#Used in the case that nested blocks occur. Crossovers within the path pair are found and resolved, while the haplotypes are updated accordingly whenever path parts are swapped
def improve_paths2(path1, path2, intE, H_A, H_B):
	newpath1 = path1[:]
	newpath2 = path2[:]
	#store original paths for later use to not lose this information when updating the paths
	oldpath1 = newpath1[:]
	oldpath2 = newpath2[:]
	
	#store haplotypes in strings and copy into lists to store for later use
	ha_s = str(H_A)
	hb_s = str(H_B)
	oldhaplo1 = ha_s[:]
	oldhaplo2 = hb_s[:]

	twosidedswitches = 0
	onesidedswitches = 0
	noswitch = 0
	samerow = 0
	bothrowssame = 0
	onesidedswitch_list = []
	bothrowssamelist = []
	twosidedswitchlist = []
	
	#for every internal block (when dealing with nested blocks), switches are found around the starting position of the internal block
	for pair in intE:
		if (pair[0] != len(path1)-1 and pair[1] < len(path1)):
			i = pair[0]
			j = pair[1]
			#a switch is found when both paths cross around position i. Then, the corresponding path parts are swapped between i and j, thus, within the internal block only. 
			#When paths are swapped, the corresponding haplotype fragment is also swapped to ensure the correct updating of the haplotypes
			if (path1[i-1] != path1[i] and path2[i-1] != path2[i] and (path1[i-1] == path2[i] or path2[i-1] == path1[i])):
				newpath1[i:j+1] = oldpath2[i:j+1]
				newpath2[i:j+1] = oldpath1[i:j+1]
				ha_s = oldhaplo1[0:i] + oldhaplo2[i:j+1] + oldhaplo1[j+1:]
				hb_s = oldhaplo2[0:i] + oldhaplo1[i:j+1] + oldhaplo2[j+1:]
				oldpath1 = newpath1[:]
				oldpath2 = newpath2[:]
				oldhaplo1 = ha_s[:]
				oldhaplo2 = hb_s[:]
	return((newpath1,newpath2), (ha_s, hb_s), twosidedswitches, onesidedswitches, noswitch, samerow, bothrowssame, onesidedswitch_list, bothrowssamelist, twosidedswitchlist)

#Used in blockparser_inference. Updates the haplotypes accordingly to the paths and resolves any existing switches. Returns two haplotypes that are optimal for the pair of paths.
def StrandSeq_improvepaths(pathfile, haplofile, intEfile):
	#reads the paths from file and stores them in lists
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
	path1, path2 = [], []
	for i in range(0, len(helppath1)-1):
		path1.append(int(float(helppath1[i])))
	for i in range(0, len(helppath2)-1):
		path2.append(int(float(helppath2[i])))	

	#reads the haplotypes from file and stores them in strings
	H_A, H_B = "",""
	with open(haplofile) as filestream:
		for i,line in enumerate(filestream):
			if (i==0):
				H_A = line.strip('\n')
			if (i==1):
				H_B = line.strip('\n')

	#reads the number of internal block positions from file into a list
	intE = []
	with open(intEfile) as filestream:
		for i,line in enumerate(filestream):
			if (i==0):
				for num in line.strip().split(','):					
					intE.append((int(num), int(num)))
	#resolve any existing switches within the paths and update the haplotypes accordingly
	(newpath1, newpath2) = improve_paths2(list(path1),list(path2), intE, H_A, H_B)[0]
	(newhaplo1, newhaplo2) = improve_paths2(list(path1),list(path2), intE, H_A, H_B)[1]
	return((newpath1, newpath2),(newhaplo1, newhaplo2))

#comparison to the original haplotypes: Switch error rate within blocks and at block ends is computed.
def compute_switcherrors(haplo1, originalhaplo1, haplo2, originalhaplo2, E_file):
	#reads the set of block ending positions into lists. 
	E, oldE = [], []
	with open(E_file) as filestream:
		for i,line in enumerate(filestream):
			if (i==0):
				for num in line.strip().split(','):					
					E.append(int(num))
					oldE.append(int(num))
	haplo = {}
	index = {}

	switchsequence_target,switchsequence_original = [], []
	switchsequence_target2,switchsequence_original2 = [], []
	switchnumber = 0
	switchnumber_blockends = 0
	testswitch = 0
	flipnumber = 0
	flipnumber2 = 0
	switcherrorlist = []	
	fliperrorlist = []
	withinlist = []
	newhaplo1, newhaplo2, neworhaplo1, neworhaplo2 = "","","",""
	oldindex, index = [], []

	#for the switch sequences, only heterozygous sites are relevant. These are created in newhaplo1,newhaplo2. For later use, the original ending positions are stored in the dictionary haplo, mapping the original end positions to the new ones.
	for i in range(0, len(haplo1)):
		if (haplo1[i] == haplo2[i] and originalhaplo1[i] == originalhaplo2[i]):
			for j in range(len(E)):
				if (oldE[j]>=i):				 
					E[j] = E[j]-1
					haplo[oldE[j]] = E[j]	
		else:
			newhaplo1 += str(haplo1[i])
			newhaplo2 += str(haplo2[i])
			neworhaplo1 += str(originalhaplo1[i])
			neworhaplo2 += str(originalhaplo2[i])

	#switch sequences are created for both haplotypes and both original haplotypes. 
	for i in range(0, len(newhaplo1)-1):
		switchsequence_target.append((int(newhaplo1[i])^int(newhaplo1[i+1])))
		switchsequence_target2.append((int(newhaplo2[i])^int(newhaplo2[i+1])))
		switchsequence_original.append((int(neworhaplo1[i])^int(neworhaplo1[i+1])))
		switchsequence_original2.append((int(neworhaplo2[i])^int(neworhaplo2[i+1])))
	#by comparing the switch sequences of original and target (one of each is sufficient since only heterozygous sites are used), the number of switch errors is computed
	for i in range(0,len(switchsequence_target)):
		if (i not in E):
			if (switchsequence_target[i] != switchsequence_original[i]):
				switchnumber += 1
				withinlist.append(i)
		if (i in E):
			if (switchsequence_target[i] != switchsequence_original[i]):
				switchnumber_blockends += 1		
				switcherrorlist.append(i)
	liste = switcherrorlist
			
	#resulting error numbers are computed
	SER = (switchnumber/float(len(haplo1)))*100
	SER_B = (switchnumber_blockends/float(len(oldE)))*100
	return(SER, SER_B,switchnumber, switchnumber_blockends, testswitch, switcherrorlist, flipnumber, withinlist)	

def computeindex(haplo, liste):
	newlist= []
	doublelist = []
	for i in liste:
		for key in list(haplo.keys()):
			if (i not in doublelist): 
				if (haplo[key]==i):
					newlist.append(key)
					doublelist.append(i)
	return(newlist)

#Used in blockparser. Updates the haplotypes accordingly to the paths and resolves any existing switches. Returns two haplotypes that are optimal for the pair of paths. Difference to
#new_haplotypes: In the improve_paths function to resolve any switches, internal blocks are included.
def new_haplos(haplofile, E, pathfile, intE):
	#read the haplotypes from file and store them in strings
	H_A, H_B = "",""
	with open(haplofile) as filestream:
		for i,line in enumerate(filestream):
			if (i==0):
				H_A = line.strip('\n')
			if (i==1):
				H_B = line.strip('\n')

	#read the optimal DP paths from file and store in lists		
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
	
	#resolve any switches that exist within the paths and update the haplotypes accordingly
	(newhaplo1, newhaplo2) = improve_paths2(list(path1),list(path2), intE, H_A, H_B)[1]
	return((newhaplo1, newhaplo2))
	
#Used in blockparser_evaluate. Updates the haplotypes accordingly to the paths and resolves any existing switches. Returns two haplotypes that are optimal for the pair of paths.
def new_haplotypes(haplofile, E, pathfile, intE):
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
	
	#resolve any switches present in the paths and update the haplotypes accordingly
	(newhaplo1, newhaplo2) = improve_paths(list(path1),list(path2),E, H_A, H_B)[1]
	return((newhaplo1, newhaplo2))
