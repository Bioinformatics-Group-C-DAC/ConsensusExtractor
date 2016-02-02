#!/bin/bash

##########################################################################
#                                                                        #
#                      consensusExtractor.sh                             #
#                           Version 0.91                                 #
##########################################################################


input_path_bam="${1:?[Usage: consensusExtractor.sh <Input_Path for bam files> <ref_genome.fa> <chr:start_Position-End_Position> <Output_filename>]}"

ref_filename="${2:?[Usage: consensusExtractor.sh <Input_Path for bam files> <ref_genome.fa> <chr:start_Position-End_Position> <Output_filename>]}"

chr_start_end="${3:?[Usage: consensusExtractor.sh <Input_Path for bam files> <ref_genome.fa> <chr:start_Position-End_Position> <Output_filename>]}"

output_filename="${4:?[Usage: consensusExtractor.sh <Input_Path for bam files> <ref_genome.fa> <chr:start_Position-End_Position> <Output_filename>]}"

ref="${ref_filename##*/}"
name_ref="${ref%.*}"
#echo "$name_ref"


if [[ "$chr_start_end" =~ ":" && "$chr_start_end" =~ "-" ]]; then

chr=$(echo "$chr_start_end" | perl -pe 's/(.*?):(.*?)-(.*)/$1/g')

start_co=$(echo "$chr_start_end" | perl -pe 's/(.*?):(.*?)-(.*)/$2/' | sed -e 's/,//g')

end_co=$(echo "$chr_start_end" | perl -pe 's/(.*?):(.*?)-(.*)/$3/' | sed -e 's/,//g')

else
        echo "Check the chr:start-end coordinate format as ENSEMBL or IGV"
        exit 1
fi
#echo -e "$chr $start_co $end_co"

if [[ $start_co -gt $end_co  ]]; then
        echo "Check the Start Position - Start Position should be less than End Position"
        exit 1
else
        pos_chr_start_end=$(echo "$chr":"$start_co"-"$end_co")
        name_chr_start_end=$(echo "$chr"_"$start_co"_"$end_co")
#        echo "$name_chr_start_end"
#        echo "$pos_chr_start_end"

fi

echo "Your command given as "
echo " sh consensusExtractor.sh "$input_path_bam" "$ref_filename" "$pos_chr_start_end" "$4"_"$name_chr_start_end" "

input_bam_dir=$(echo "$input_path_bam")
#echo "$input_bam_dir"
mkdir -p "$input_bam_dir/tmp/"
mkdir -p "$input_bam_dir/consensus_outputs/"
cd $input_bam_dir

for i in $( ls *.bam)
do

filename="${i%.*}"

samtools view -b "$filename".bam $pos_chr_start_end > "$filename"_"$name_chr_start_end".bam

samtools sort "$filename"_"$name_chr_start_end".bam "$filename"_"$name_chr_start_end"_sorted
samtools index "$filename"_"$name_chr_start_end"_sorted.bam

samtools mpileup -u -d 100000 -f "$name_ref".fa "$filename"_"$name_chr_start_end"_sorted.bam > "$filename"_"$name_chr_start_end".upileup

bcftools view -cg "$filename"_"$name_chr_start_end".upileup > "$filename"_"$name_chr_start_end".vcf

bgzip "$filename"_"$name_chr_start_end".vcf

tabix -p vcf "$filename"_"$name_chr_start_end".vcf.gz

samtools faidx "$name_ref".fa $pos_chr_start_end > "$name_chr_start_end"_"$name_ref"

cat "$name_chr_start_end"_"$name_ref" | vcf-consensus "$filename"_"$name_chr_start_end".vcf.gz | sed -e "s/^>.*/>${filename}_${name_chr_start_end}/g" >> $input_bam_dir/consensus_outputs/$4_"$name_chr_start_end"


mv "$filename"_"$name_chr_start_end".bam "$filename"_"$name_chr_start_end"_sorted.bam "$filename"_"$name_chr_start_end"_sorted.bam.bai $input_bam_dir/tmp/

mv "$name_chr_start_end"_"$name_ref" $input_bam_dir/consensus_outputs/

done

mv *.upileup *.vcf.gz *.tbi $input_bam_dir/tmp/

