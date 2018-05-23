Files in this folder:
The scripts blockparser_evaluate.py, blockparser.py and blockparser_inference.py call the necessary functions to execute the block parsing process
either with test data and corresponding results or without the testing, or to be used with StrandSeq data. Functions for the computation of blocks,
updating of paths and haplotypes and evaluation are stored in helperfunctions.py. 
blockparser.pyx calls the dynamic programming step for the computation of the scoring matrix. This is performed within DP_matrix.cpp. 
setup.py contains the necessary data for compilation of the Cython code. 

1. Compiling:
	python setup.py build_ext --inplace

2. Block parsing with subsequent evaluation and computation of switch errors:
- call: python blockparser_evaluate.py [inputlist]
- [inputlist] must contain:
	1. [target.vcf] The target file containing the phased haplotype data in a block structure (PS tag needed)
	2. [original.vcf] A file containing the ground truth haplotypes used for evaluation
	3. [refpanel.txt] The reference panel
	4. [mutation] The mutation parameter (float)
	5. [output] Filename for the updated VCF file
	
	Test files are present in the folder "Test":
	1. NA12878_chr22.vcf
	2. NA12878_chr22_original.vcf
	3. reffile_chr22_CEU.txt
	
	Example call:
	python blockparser_evaluate.py Test/NA12878_chr22.vcf Test/NA12878_chr22_original.vcf Test/reffile_chr22_CEU.txt 1 NA12878_chr22_updated.vcf
	
	The number of switch errors is written to standard output.


3. Block parsing without testing: 
- call: python blockparser.py [inputlist]
- [inputlist] must contain:
	1. [target.vcf] The target file containing the phased haplotype data in a block structure (PS tag needed)
	2. [refpanel.txt] The reference panel
	3. [mutation] The mutation parameter (float)
	4. [output] Filename for the updated VCF file
	
	Test files are present in the folder "Test":
	1. NA12878_chr22.vcf
	2. reffile_chr22_CEU.txt

	Example call:
	python blockparser.py Test/NA12878_chr22.vcf Test/reffile_chr22_CEU.txt 1 NA12878_chr22_updated.vcf
	

4. Block parsing for using sparse data from Strand-Seq 
- call: python blockparser_inference.py [inputlist]
- [inputlist] must contain: 
		1. [original.vcf] For testing purposes: The original file containing sites that are missing in the target
		2. [target.vcf] The target file with the sparse haplotype data
		3. [refpanel.txt] The reference panel
		4. [mutation] The mutation parameter (float)
		5. [output] Filename for the updated VCF file
	
	Test files are present in folder "Test_StrandSeq": 
	1. NA19240_chr22_original.vcf: VCF containing the ground truth phasing 
	2. NA19240_chr22_target.vcf: VCF that needs haplotype inference
	3. refpanel_chr22_YRI.txt: Reference panel for the YRI population
	
	Example call: 
	python blockparser_inference.py Test_StrandSeq/NA19240_chr22_original.vcf Test_StrandSeq/NA19240_chr22_target.vcf Test_StrandSeq/refpanel_chr22_YRI.txt 1 Test_StrandSeq/NA19240_chr22_newlyphased.vcf
	

5. For the computation of the reference panel, the following script was used:
- call: python compute_ref.py [inputlist]
Where [inputlist]  contains:
1: The VCF file to be used as a reference panel
2: The target VCF with the phased data that is to be updated
3: filename for output TXT file 

