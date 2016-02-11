#!/bin/bash

##########################################################################
#                                                                        #
#                      consensusExtractor.sh                             #
#                           Version 0.9.1                                #
##########################################################################


input_path_bam="${1:?[Usage: consensusExtractor.sh <Input_Path for bam files> <Input_Path for ref_genome.fa> <chr:start_Position-End_Position> <Output_filename>]}"

ref_filename="${2:?[Usage: consensusExtractor.sh <Input_Path for bam files> <Input_Path for ref_genome.fa> <chr:start_Position-End_Position> <Output_filename>]}"

chr_start_end="${3:?[Usage: consensusExtractor.sh <Input_Path for bam files> <Input_Path for ref_genome.fa> <chr:start_Position-End_Position> <Output_filename>]}"

output_filename="${4:?[Usage: consensusExtractor.sh <Input_Path for bam files> <Input_Path for ref_genome.fa> <chr:start_Position-End_Position> <Output_filename>]}"


bam_files=$(ls $input_path_bam)
#echo "$bam_files"

if [[ "$bam_files" =~ ".bam" ]]; then

bam_count=$(ls $input_path_bam/*.bam | wc -l )
#echo "$bam_count"

else
	echo "Please check input path for bam files"
	exit 1

fi
bam_count=$(ls $input_path_bam/*.bam | wc -l )
if [[ "$bam_files" =~ ".bam" ]]; then
bai_count=$(ls $input_path_bam/*.bam.bai | wc -l )
#echo "$bai_count"
	if [[ $bam_count == $bai_count  ]]; then
        echo ".bam and index file are there" > /dev/null
else 
       echo "Please check respective index .bai file for .bam file "
       exit 1
		fi
	fi

ref="${ref_filename##*/}"
name_ref="${ref%.*}"
#echo "$name_ref"

if [[ -f "$ref_filename" ]]; then
        echo "reference genome .fa is there" > /dev/null
else
        echo "Please check input path and name for reference genome in .fa format in reference genome directory"
        exit 1
fi


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
echo " sh consensusExtractor.sh "$input_path_bam" "$ref_filename" "$pos_chr_start_end" "$4" "

input_bam_dir=$(echo "$input_path_bam")
#echo "$input_bam_dir"
temp=$(mktemp -d $name_chr_start_end.XXXXXX)

for i in $( ls $input_bam_dir/*.bam)
do

bam="${i##*/}"
#echo "$bam"

filename="${bam%.*}"
#echo "$filename"

samtools view -b "$i" $pos_chr_start_end > "$filename"_"$name_chr_start_end".bam

samtools sort "$filename"_"$name_chr_start_end".bam "$filename"_"$name_chr_start_end"_sorted

samtools index "$filename"_"$name_chr_start_end"_sorted.bam

samtools mpileup -u -d 100000 -f "$ref_filename" "$filename"_"$name_chr_start_end"_sorted.bam > "$filename"_"$name_chr_start_end".upileup

bcftools view -cg "$filename"_"$name_chr_start_end".upileup > "$filename"_"$name_chr_start_end".vcf

bgzip "$filename"_"$name_chr_start_end".vcf

tabix -p vcf "$filename"_"$name_chr_start_end".vcf.gz

samtools faidx "$ref_filename" $pos_chr_start_end > "$name_chr_start_end"_"$name_ref"

cat "$name_chr_start_end"_"$name_ref" | vcf-consensus "$filename"_"$name_chr_start_end".vcf.gz | sed -e "s/^>.*/>${filename}_${name_chr_start_end}/g" >> "$4"


mv "$filename"_"$name_chr_start_end".bam "$filename"_"$name_chr_start_end"_sorted.bam "$filename"_"$name_chr_start_end"_sorted.bam.bai $temp

done

cat "$name_chr_start_end"_"$name_ref" "$4" >> "$4" 2> /dev/null
sed -i "s/^>${pos_chr_start_end}/>${name_ref}_${pos_chr_start_end}/g" "$4"

mv *.upileup *.vcf.gz *.tbi "$name_chr_start_end"_"$name_ref" $temp
rm -rf $temp
