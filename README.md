# Consensus Extractor

ConsensusExtractor - writtten by Bioinformatics Group CDAC (Centre For Development Of Advanced Computing), Pune 
in collboration with Pirbright Institute, UK

 _**Authors - Neeraj Bharti, Yatish B. Patil, Sunitha Manjari, Dr. Jan Kim**_

 Consesnsus Extractor - gives user a consensus sequence of a selected region from bam files 
 in comaprison with reference genome file, where it integrates data from mpileup, vcf file to fasta format
 which is limited to substitution, insertion sequences and the output is in consensus fasta format.

#  Dependencies required
*  bash shell(/bin/bash)
*  perl(linux)
*  samtools-0.1.19
*  vcftools_0.1.11
*  tabix-0.2.6 also include bgzip

# Inputs
 * Please mention user Input path for .bam  - must contain sorted bam files with index files.
 * Please mention user Input path for reference genome file in .fa format < Input path  for reference genome file name>.fa (.fa format only)

# Output 
Output will be created in current working directory of script - output name will be as given by user.

The conflicting variants message will be displayed if there is an insertion in one sequence and the same region is deleted in other sequence. 

User can use this consensus output in fasta format for alignment analysis.

# Note
after downloading scripts, please make it executable and follow the command below

`Usage:`
` ./consensusExtractor.sh <Full Input_Path for bam files> <Full Input path for genome reference file> <chr:Start_coordinate-End_Coordinate> <Output Filename> `

# Commands Used as example

User must give &lt;chromosome number:start_coordinate-End_coordinate&gt; in ENSEMBL or IGV format(with or without comma seperated for start and end coordinates)

# Examples

` ./consensusExtractor.sh /input/bam_dir/ /input/ref/ref_genome.fa 3:14,534-3,45,987 output_consensus.fasta `

` ./consensusExtractor.sh /input/bam_files/ /input/test/ref/ref_genome.fa 3:143567-234569 xyz_file `

