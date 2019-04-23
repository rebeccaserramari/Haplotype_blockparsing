"""
Use a reference panel to connect the haplotype blocks

"""
from __future__ import absolute_import
import logging
from .core import compute_referencepanel, create_blocks, compute_haplotypes, update_haplotypes, compute_referencepanel_heuristic, compute_referencepanel_heuristic2
from .customcontainer import DefaultOrderedDict
import sys
from cyvcf2 import VCF
import platform
import resource
from .core import scoring_computation
from collections import defaultdict
from operator import itemgetter

import collections
import math
from . import __version__
from .timer import StageTimer
from .utils import detect_file_format
logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('variant_file', metavar='BLOCKS', help='VCF file with phased haplotype blocks for one chromosome')
	arg('reference_file', metavar='REFERENCE', help='VCF file that serves as reference panel. This is used to connect haplotype blocks. ')
	arg('output_file', metavar="OUTPUT",help='Output file in VCF format')
	arg('--include-unphased',dest = 'include_unphased', default=False, action='store_true', help = 'Include unphased variant positions for inferring the haplotypes at these sites (Default: Use phased variants only)')
	arg('--heuristic',dest = 'heuristic', default=False, action='store_true', help = 'Use a heuristic to choose the most likely reference haplotypes (Default: Use all given samples as references)')
	arg('--heuristic2',dest = 'heuristic2', default=False, action='store_true', help = 'Use a heuristic to choose the most likely reference haplotypes (Default: Use all given samples as references)')

	
def validate(args, parser):
	pass

	
	
def included_variants(ref_file, variant_file, include_unphased):
	ref_set = set()
	chromosome_reference_set = set()
	refvcf = VCF(ref_file)
	for variant in refvcf:
		ref_set.add(variant.POS)
		chromosome_reference_set.add(variant.CHROM)
	#for multiple chromosomes, raise an error
	if (len(chromosome_reference_set)>1):
		logger.error("The given reference file is supposed to contain variants from one chromosome only!")
		sys.exit(1)
	#create the list of variants to be considered	
	variant_set = set()
	chromosome_variant_set = set()
	vcf = VCF(variant_file)
	haplo1,haplo2 = "",""
	variant_index = -1
	unknown_variants = {}	
	unknown_variantlist = []	
	all_variant_set = set()
	testlist1, testlist2 = [], []
	for v in vcf:
		#use variant if it is present in the panel, heterozygous and phased
	#	if (include_unphased==False) and (v.POS in ref_set) and (v.gt_types[0] == 1) and (v.gt_phases[0] == True) and (v.is_snp):	
		if (include_unphased==False) and (v.gt_types[0] == 1) and (v.gt_phases[0] == True) and (v.is_snp):	
			variant_index += 1
			if (v.CHROM not in chromosome_reference_set):
				logger.error("Reference and target file must contain variants from the same chromosome!")
				sys.exit(1)
			chromosome_variant_set.add(v.CHROM)
			if (v.POS in ref_set):
				variant_set.add(v.POS)
				testlist1.append(v.POS)
			else:
				unknown_variants[variant_index] = (v.genotypes[0][0],v.genotypes[0][1], v.POS)
				unknown_variantlist.append(v.POS)
			all_variant_set.add(v.POS)
			testlist2.append(v.POS)
		elif include_unphased and (v.POS in ref_set) and (v.gt_types[0] == 1) and (v.is_snp):
	#	elif include_unphased and (v.POS in ref_set) and (v.is_snp):
			if (v.CHROM not in chromosome_reference_set):
				logger.error("Reference and target file must contain variants from the same chromosome!")
				sys.exit(1)
			variant_set.add(v.POS)
			chromosome_variant_set.add(v.CHROM)
	print("variant index: ", variant_index)
	print("len unknown variants: ", len(unknown_variants))
	print(unknown_variants.keys())
	if (len(chromosome_variant_set)>1):
		logger.error("The given variant file is supposed to contain variants from one chromosome only!")
		sys.exit(1)
	if not variant_set:
		logger.error("No positions found in both the reference and the target file.")
		sys.exit(1)
	return(variant_set, unknown_variants, all_variant_set, testlist1, testlist2, unknown_variantlist)

def run_connect_heuristic(variant_file,reference_file,variant_set,output_file, include_unphased, heuristic, vcf_ref_map, variants_to_use, all_variant_set):
	print("going to run connect heuristic works")
	vcf_ref = VCF(reference_file)
	iterations = len(vcf_ref.samples)//25
	most_frequent_refs = []
	most_frequent_refs2 = []
	samples_to_use = set()
	map_result_samples = []
	print("iterations: ", iterations)
	for it in range(iterations+1):
		logger.info("Iteration %s running", it)
		fileformat = detect_file_format(reference_file)
		if fileformat == 'VCF':
			logger.info("Computing reference panel.")
			start = it*25
			end = (it+1)*25
			if (it==iterations):
				end = len(vcf_ref.samples)
			print("start and end: ",start, end)
			compute_referencepanel_heuristic(reference_file,variant_file,variants_to_use, "ref.txt", start,end)
		else:
			logger.error("Wrong file format. Please enter a reference file in VCF format.")

		#reads the phase set information from the VCF and returns a dictionary containing variant positions and its phase set, as well as block ending positions and internal blocks 
		logger.info("Reading phased blocks")
		(E_whole, haplo, intE) = create_blocks(variant_file, "blockends.txt", variant_set, all_variant_set)
		#writes the original haplotypes from the target into strings
	#	compute_haplotypes(variant_file, "haplotypes.txt", variant_set)
		(haplo1,haplo2,usedPos) = compute_haplotypes(variant_file, "haplotypes.txt", variant_set)

		#performs the computation of the scoring matrix (dynamic programming) and the backtracing to find the optimal pair of paths through the panel. Paths and corresponding costs are stored in files.
		logger.info("Dynamic programming: Performing haplotype block parsing")
		scoring_computation("haplotypes.txt", "blockends.txt", "ref.txt", 8, "costs.txt", "paths.txt")
		
		#compute the likeliest references: Those that appear for the longest sequences in the paths	
		helppath1, helppath2 = [], []
		with open("paths.txt") as filestream:
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
		print("length of path: ", len(path1))
		

		print("length path1: ", len(path1))
		print("length path2: ", len(path2))
		ref_path1, ref_path2 = defaultdict(), defaultdict()				
		for element in path1:
			if (element+start) not in ref_path1:
				ref_path1[element+start] = 0
			else:
				ref_path1[element+start] += 1
		sorted_refs_path1 = sorted(ref_path1.items(), key=itemgetter(1), reverse=True)
		for element in path2:
			if (element+start) not in ref_path2:
				ref_path2[element+start] = 0
			else:
				ref_path2[element+start] += 1
		sorted_refs_path2 = sorted(ref_path2.items(), key=itemgetter(1), reverse=True)
		print("number of references in path1: ", len(sorted_refs_path1))
		print("number of references in path2: ", len(sorted_refs_path2))
	#	most_frequent_refs.extend(sorted(sorted_refs_path1+sorted_refs_path2, key=itemgetter(1), reverse=True))
		most_frequent_refs.extend(sorted(sorted_refs_path1, key=itemgetter(1), reverse=True))
		print("number of most frequent references: ", len(most_frequent_refs))
		
		sums_path1, sums_path2 = defaultdict(), defaultdict()
		lengths = []
		counter = 0
		for index, element in enumerate(path1):
			if (index == 0):
				counter=1
				continue
			else:
				if (element == path1[index-1]):
					counter += 1
				else:
					lengths.append((path1[index-1], counter))
					counter = 1
			if (index == len(path1)-1):
				lengths.append((path1[index], counter))
		for element in lengths:
			if (element[0]+start) not in sums_path1:
				sums_path1[element[0]+start] = 0
			else:
				sums_path1[element[0]+start] += element[1]
		sorted_sums_path1 = sorted(sums_path1.items(), key=itemgetter(1), reverse=True)
		
		lengths_path2 = []
		counter_path2 = 0
		for index, element in enumerate(path2):
			if (index == 0):
				counter_path2=1
				continue
			else:
				if (element == path2[index-1]):
					counter_path2 += 1
				else:
					lengths_path2.append((path2[index-1], counter_path2))
					counter_path2 = 1
			if (index == len(path2)-1):
				lengths_path2.append((path2[index], counter_path2))
		for element in lengths_path2:
			if (element[0]+start) not in sums_path2:
				sums_path2[element[0]+start] = 0
			else:
				sums_path2[element[0]+start] += element[1]
		sorted_sums_path2 = sorted(sums_path2.items(), key=itemgetter(1), reverse=True)	
			
		most_frequent_refs2.extend(sorted(sorted_sums_path1, key=itemgetter(1),reverse=True))
		print("number of most frequent references 2: ", len(most_frequent_refs2))
	print("number of refs: ", len(most_frequent_refs))
	most_frequent_refs.sort(key=itemgetter(1), reverse=True)
	print("first references: ", most_frequent_refs[:20])
	most_frequent_refs2.sort(key=itemgetter(1), reverse=True)
	print("first references: ", most_frequent_refs2[:20])
#	samples_to_use = list(set([i[0]//2 for i in most_frequent_refs[:500]]))
	samples_to_use = list(set([i[0]//2 for i in most_frequent_refs2][:800]))
	counter = 0
	tmp1 = [i[0] for i in most_frequent_refs]
	tmp2 = [i[0] for i in most_frequent_refs2]
	for i in tmp1:
		if i not in tmp2:
			counter += 1
	print("differences: ", counter)
	print("number of samples: ", len(samples_to_use))
	return(samples_to_use)


			
def run_connect(variant_file,reference_file,output_file, include_unphased, heuristic, heuristic2):
	var = VCF(variant_file)
	logger.info("This is WhatsHap (reference-based phasing) %s running under Python %s", __version__, platform.python_version())

	timers = StageTimer()
	timers.start('overall')
	with timers("variant_list"):
		logger.info("Choosing suitable variant positions.")
		(variant_set, unknown_variants, all_variant_set, testlist1, testlist2, unknown_variantlist) = included_variants(reference_file, variant_file,include_unphased)

	vcf = VCF("Test/NA12878_chr22_original.vcf")
	vcf_ref = VCF(reference_file)
	vcf_ref_map = defaultdict()
	for variant in vcf_ref:
		vcf_ref_map[variant.POS] = variant
	print("number of all variants: ", len(vcf_ref_map.keys()))
	variants_to_use = [value for key,value in vcf_ref_map.items() if key in variant_set]
	print("number of variants to use: ", len(variants_to_use))	
	print("number of all variants regardless of ref: ", len(all_variant_set))
	if heuristic:
		samples_to_use = run_connect_heuristic(variant_file,reference_file,variant_set, output_file, include_unphased, heuristic, vcf_ref_map, variants_to_use, all_variant_set)
	elif heuristic2:
#		refvcf = VCF(reference_file)
#		ref_set = set()
#		for variant in refvcf:
#			ref_set.add(variant.POS)
#			
#		variant_heuristic = set()
#		vcf = VCF(variant_file)
#		for v in vcf:
#			if (v.POS in ref_set) and v.is_snp:
#				variant_heuristic.add(v.POS)
#				
#		refvcf = VCF(reference_file)
#		vcf_ref_map = defaultdict()
#		for variant in refvcf:
#			vcf_ref_map[variant.POS] = variant
#		variants_heur = [value for key,value in vcf_ref_map.items() if key in variant_heuristic]
#		
#		vcf_target = VCF(variant_file)
#		hap1 = ""
#		hap2 = ""
#		for v in vcf_target:
#			if (v.POS in variant_heuristic):
#				comp1 = 0
#				comp2 = 0
#				if (len(v.gt_bases[0].split('|')) > 1):
#					leftbase = v.gt_bases[0].split('|')[0]
#					rightbase = v.gt_bases[0].split('|')[1]
#				elif(len(v.gt_bases[0].split('/')) > 1):
#					leftbase = v.gt_bases[0].split('/')[0]
#					rightbase = v.gt_bases[0].split('/')[1]
#				if (leftbase == v.REF):
#					comp1 = 0
#				elif (leftbase == v.ALT[0]):
#					comp1 = 1
#				if (rightbase == v.REF):
#					comp2 = 0
#				elif (rightbase == v.ALT[0]):
#					comp2 = 1
#				hap1 += str(comp1)	
#				hap2 += str(comp2)
#		print("number of variants_heur: ", len(variants_heur))
#		print("length of hap1: ", len(hap1))
#		samples_to_use = compute_referencepanel_heuristic2(reference_file, variant_file, variants_heur, output_file, hap1,hap2)
	
		hap1, hap2, positions = compute_haplotypes(variant_file, "haplotypes.txt", variant_set)
		samples_to_use = compute_referencepanel_heuristic2(reference_file, variant_file, variants_to_use, output_file, hap1,hap2)
	else:
		vcf_ref = VCF(reference_file)
		samples_to_use = list(range(len(vcf_ref.samples)))


	with timers("ref_panel"):
		fileformat = detect_file_format(reference_file)
		if fileformat == 'VCF':
			logger.info("Computing reference panel.")
			compute_referencepanel(reference_file,variant_file,variant_set, "ref.txt", samples_to_use, variants_to_use)						
		else:
			logger.error("Wrong file format. Please enter a reference file in VCF format.")
	#reads the phase set information from the VCF and returns a dictionary containing variant positions and its phase set, as well as block ending positions and internal blocks 
	with timers("creating_blocks"):
		logger.info("Reading phased blocks")
		
		(E_whole, haplo, intE) = create_blocks(variant_file, "blockends.txt", variant_set, all_variant_set)
		print("number of blocks: ", len(E_whole))
		print("number of variants: ", len(variant_set))
		#writes the original haplotypes from the target into strings
		(haplo1,haplo2,usedPos) = compute_haplotypes(variant_file, "haplotypes.txt", variant_set)

	#performs the computation of the scoring matrix (dynamic programming) and the backtracing to find the optimal pair of paths through the panel. Paths and corresponding costs are stored in files.
	with timers("dynamic_programming"):
		logger.info("Dynamic programming: Performing haplotype block parsing")
		scoring_computation("haplotypes.txt", "blockends.txt", "ref.txt", 8, "costs.txt", "paths.txt")

	#update the haplotypes accordingly to the paths and resolve any existing switches. Returns two haplotypes that are optimal for the pair of paths.
	with timers("creating_haplotypes"):
		logger.info("Assembling the resulting haplotypes")
		((newhaplo1, newhaplo2), (newpath1,newpath2), usedPos, cutpositions, H_A) = update_haplotypes("haplotypes.txt", E_whole, "paths.txt", intE,unknown_variants, testlist1, testlist2, usedPos)
	print("cut positions: ", cutpositions)

	#org = before switching
	org = ""
	for i in range(len(H_A)):
		if i not in unknown_variants.keys():
			org += H_A[i]
#	print("before switching, shortened: ", org)

	control = ""
	for i in range(len(newhaplo1)):
		if i not in unknown_variants.keys():
			control += newhaplo1[i]
#	print("haplo after switching, shortened: ", control)	
#for testing purposes:
	#write out the original haplotypes and the new ones
	vcf = VCF("Test/NA12878_chr22_original.vcf")
	counter = -1
	originalhaplo1, originalhaplo2 = "", ""
	originalPos = []
	for v in vcf:
		comp1, comp2 = 0,0
	#	if v.POS in variant_set:
		if v.POS in all_variant_set:
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
			originalPos.append(v.POS)
#	print("original: ", originalhaplo1)
	cut_starts = [i[0] for i in cutpositions]
	cut_ends = [i[1] for i in cutpositions]
	target = VCF(variant_file)
	varcounter = -1
	cut_startpos, cut_endpos = [], []
	test = ""
	for v in target:
		if v.POS in all_variant_set:
			varcounter += 1
			if varcounter in cut_starts:
				cut_startpos.append(v.POS)
			if varcounter in cut_ends:
				cut_endpos.append(v.POS)
		if v.POS in variant_set:
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
			test += str(comp1)
#	print("test: ", test)
	print("length new haplos: ", len(newhaplo1))
	print("length original haplos: ", len(originalhaplo1))

#	print("originalhaplo1: ", originalhaplo1)
#	print("newhaplo1: ", newhaplo1)
	print("newhaplo1:      ", newhaplo2[:60])
	print("usedPos: ", usedPos[:60])
	print("originalpos: ", originalPos[:60])
	print("length usedPos and originalPos: ", len(usedPos),len(originalPos))
	
	
#	for i in range(len(usedPos)):
#		if (usedPos[i] != originalPos[i]):
#			print("DIFFERENCE at : ", usedPos[i], i)
	#create switch sequences for both
	switchsequence_target,switchsequence_original = [], []
	switchnumber = 0
	for i in range(0, len(newhaplo1)-1):
		switchsequence_target.append((int(newhaplo1[i])^int(newhaplo1[i+1])))
		switchsequence_original.append((int(originalhaplo1[i])^int(originalhaplo1[i+1])))
	
	switchsequence_target = ''.join(('0' if newhaplo1[i-1] == newhaplo1[i] else '1') for i in range(1, len(newhaplo1)))
	switchsequence_original = ''.join(('0' if originalhaplo1[i-1] == originalhaplo1[i] else '1') for i in range(1, len(originalhaplo1)))
	#compute switch errors and output them	
	cut_positions = []
	for i in range(0,len(switchsequence_target)):
		if (switchsequence_target[i] != switchsequence_original[i]):
			cut_positions.append(i)
			switchnumber += 1	
	SER = (switchnumber/len(switchsequence_target))*100
	print("SER: ", SER)
	print("switchnumber: ", switchnumber)
#	print("cut positions: ", cut_positions)
	with timers("output"):
		logger.info("Writing output file. This action may take a few minutes.")
		new_haplo_dict = DefaultOrderedDict()
		i = 0
		keylist = []
		for key in sorted(haplo.keys()):
			keylist.append(key)
			new_haplo_dict[key] = (newhaplo1[i],newhaplo2[i])
			i += 1
		key_dict = DefaultOrderedDict()
		print("keylist: ", keylist[:50])
		if (16057248 in keylist):
			print("is in keylist, position: ", sorted(keylist).index(16057248))
		if (16059753 in keylist):
			print("is in keylist, position: ", sorted(keylist).index(16059753))
		print("variant set: ", sorted(list(all_variant_set))[:50])
		print("number of new haplo dict: ", len(new_haplo_dict.keys()))
		print("number of haplo dict: ", len(sorted(haplo.keys())))
	#	print("variant set: ", variant_set)
		with open(variant_file) as f:
			with open(output_file, 'w') as of:
				linecounter = 0
				counter = 0
				current = 0
				for line in f:
					linecounter += 1			
					key_found = False
					if line.startswith('#'):
						of.write(line)
					else:
						for key in list(new_haplo_dict.keys())[current:]:
					#	for key in list(new_haplo_dict.keys()):							
							if (str(key) in line.split('.')[0]):	
								current = list(new_haplo_dict.keys()).index(key)
								key_found = True
								counter += 1
					
								if ("0|1" in line):
								#	of.write(line.replace("0|1",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))
									if key in unknown_variantlist:
										of.write(line.rsplit(':',1)[0].replace("0|1",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1])+':'+str(key)+'\n')
									else:					
										of.write(line.replace("0|1",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))
									#	of.write(line.rsplit(':',1)[0].replace("0|1",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1])+':16057248'+'\n')			
								elif ("1|0" in line):
								#	of.write(line.replace("1|0",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))
									if key in unknown_variantlist:
										of.write(line.rsplit(':',1)[0].replace("1|0",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1])+':'+str(key)+'\n')
									else:
										of.write(line.replace("1|0",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))
									#	of.write(line.rsplit(':',1)[0].replace("0|1",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1])+':16057248'+'\n')
								break
						if (key_found == False):
						#	of.write(line.rsplit(':',1)[0]+':'+str(linecounter)+'\n')
							of.write(line)
	print("key found in cases: ", counter)
	logger.info('\n== SUMMARY ==')
	timers.stop("overall")
	if sys.platform == 'linux':
		memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		logger.info('Maximum memory usage: %.3f GB', memory_kb / 1E6)
	logger.info('Time spent choosing suitable variant positions:                      %6.1f s', timers.elapsed('variant_list'))
	logger.info('Time spent creating reference panel:                      %6.1f s', timers.elapsed('ref_panel'))
	logger.info('Time spent creating haplotype blocks:                      %6.1f s', timers.elapsed('creating_blocks'))
	logger.info('Time spent finding optimal paths through the reference panel:                      %6.1f s', timers.elapsed('dynamic_programming'))
	logger.info('Time spent assembling the connected haplotypes:                      %6.1f s', timers.elapsed('creating_haplotypes'))
	logger.info('Time spent writing output file:                      %6.1f s', timers.elapsed('output'))
	logger.info('Total elapsed time:                      %6.1f s', timers.elapsed('overall'))

def main(args):
	run_connect(**vars(args))
	#run_connect_heuristic(**vars(args))
		