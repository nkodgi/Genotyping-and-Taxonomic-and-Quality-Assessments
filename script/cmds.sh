#!/bin/bash
# Created a conda environment for genomics tools for read cleaning, assembly, and fastANI
conda create -n genomics3 -y
conda activate genomics3

#installed packages
conda install -c bioconda -c conda-forge entrez-direct sra-tools fastqc fastp spades fastANI pigz tree seqtk -y
#created a directory to pull our fastq reads
mkdir -pv ~/hw3/raw_reads
cd ~/hw3/raw_reads

#prefetched all my reads
prefetch SRR27160580 SRR27160579 SRR27160578

#used fastq-dump and a for loop to run all three
for sra in SRR27160580 SRR27160579 SRR27160578
do
    fastq-dump $sra --outdir ~/hw3/raw_reads --split-files --skip-technical
done

#checked file sizes and gunzipped the files
du -sh ~/hw3/raw_reads/*.fastq 
pigz -9fv ~/hw3/raw_reads/*.fastq

#read cleaned using fastp
fastp -i ~/hw3/raw_reads/SRR27160580_1.fastq.gz -I ~/hw3/raw_reads/SRR27160580_2.fastq.gz -o ~/hw3/raw_reads/SRR27160580_clean_1.fastq.gz -O ~/hw3/raw_reads/SRR27160580_clean_2.fastq.gz
fastp -i ~/hw3/raw_reads/SRR27160579_1.fastq.gz -I ~/hw3/raw_reads/SRR27160579_2.fastq.gz -o ~/hw3/raw_reads/SRR27160579_clean_1.fastq.gz -O ~/hw3/raw_reads/SRR27160579_clean_2.fastq.gz
fastp -i ~/hw3/raw_reads/SRR27160578_1.fastq.gz -I ~/hw3/raw_reads/SRR27160578_2.fastq.gz -o ~/hw3/raw_reads/SRR27160578_clean_1.fastq.gz -O ~/hw3/raw_reads/SRR27160578_clean_2.fastq.gz
ls ~/hw3/raw_reads #made sure fastp worked, I also looked at all the graphs in the html file to see what the genomes looked like
mkdir ~/hw3/spades
cd ~/hw3/spades

#used spades for assembly - all three files took me around 20 minutes bc I have a 2 core computer
spades.py -1 ~/hw3/raw_reads/SRR27160580_clean_1.fastq.gz -2 ~/hw3/raw_reads/SRR27160580_clean_2.fastq.gz -o ~/hw3/spades/SRR27160580 --only-assembler 1> ~/hw3/spades/SRR27160580_spades.stdout.txt 2> ~/hw3/spades/SRR27160580_spades.stderr.txt
spades.py -1 ~/hw3/raw_reads/SRR27160579_clean_1.fastq.gz -2 ~/hw3/raw_reads/SRR27160579_clean_2.fastq.gz -o ~/hw3/spades/SRR27160579 --only-assembler 1> ~/hw3/spades/SRR27160579_spades.stdout.txt 2> ~/hw3/spades/SRR27160579_spades.stderr.txt
spades.py -1 ~/hw3/raw_reads/SRR27160578_clean_1.fastq.gz -2 ~/hw3/raw_reads/SRR27160578_clean_2.fastq.gz -o ~/hw3/spades/SRR27160578 --only-assembler 1> ~/hw3/spades/SRR27160578_spades.stdout.txt 2> ~/hw3/spades/SRR27160578_spades.stderr.txt

#renamed all my contigs files to their SRA numbers
mv ~/hw3/spades/SRR27160578/contigs.fasta ~/hw3/spades/SRR27160578/SRR27160578_contigs.fasta
mv ~/hw3/spades/SRR27160579/contigs.fasta ~/hw3/spades/SRR27160579/SRR27160579_contigs.fasta
mv ~/hw3/spades/SRR27160580/contigs.fasta ~/hw3/spades/SRR27160580/SRR27160580_contigs.fasta
mv ~/hw3/spades/SRR27160578/SRR27160578_contigs.fasta ~/hw3/spades/SRR27160579/SRR27160579_contigs.fasta ~/hw3/spades/SRR27160580/SRR27160580_contigs.fasta ~/hw3/spades
grep -c '>' ~/hw3/spades/*.fasta #checked all the fasta file sizes to see how many contigs were inside!

#grabbing the type strain genome for comparison for FastANI
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/735/925/GCF_020735925.1_ASM2073592v1/GCF_020735925.1_ASM2073592v1_genomic.fna.gz -o ~/hw3/spades/GCF_020735925.1_ASM2073592v1_genomic.fna.gz
gunzip -kv ~/hw3/spades/GCF_020735925.1_ASM2073592v1_genomic.fna.gz

#did some post-assembly trimming by filtering out all contigs less than 500bps before fastANI
seqtk seq -L 500 ~/hw3/spades/SRR27160578_contigs.fasta > ~/hw3/spades/SRR27160578_filtered_contigs.fna
seqtk seq -L 500 ~/hw3/spades/SRR27160579_contigs.fasta > ~/hw3/spades/SRR27160579_filtered_contigs.fna
seqtk seq -L 500 ~/hw3/spades/SRR27160580_contigs.fasta > ~/hw3/spades/SRR27160580_filtered_contigs.fna

#fastANI for comparison to reference genomes
fastANI --query ~/hw3/spades/SRR27160578_filtered_contigs.fna --ref ~/hw3/spades/GCF_020735925.1_ASM2073592v1_genomic.fna --output ~/hw3/spades/SRR27160578_FastANI_Output.tsv
fastANI --query ~/hw3/spades/SRR27160579_filtered_contigs.fna --ref ~/hw3/spades/GCF_020735925.1_ASM2073592v1_genomic.fna --output ~/hw3/spades/SRR27160579_FastANI_Output.tsv
fastANI --query ~/hw3/spades/SRR27160580_filtered_contigs.fna --ref ~/hw3/spades/GCF_020735925.1_ASM2073592v1_genomic.fna --output ~/hw3/spades/SRR27160580_FastANI_Output.tsv

#added the alignment percent and alignment lengths for each tsv file and combined them - then I added a header for this fastANI table
awk '{alignment_percent = $4/$5*100} {alignment_length = $4*3000} {print $0 "\t" alignment_percent "\t" alignment_length}' ~/hw3/spades/SRR27160580_FastANI_Output.tsv > ~/hw3/spades/SRR27160580_FastANI_Output_With_Alignment.tsv
awk '{alignment_percent = $4/$5*100} {alignment_length = $4*3000} {print $0 "\t" alignment_percent "\t" alignment_length}' ~/hw3/spades/SRR27160579_FastANI_Output.tsv > ~/hw3/spades/SRR27160579_FastANI_Output_With_Alignment.tsv
awk '{alignment_percent = $4/$5*100} {alignment_length = $4*3000} {print $0 "\t" alignment_percent "\t" alignment_length}' ~/hw3/spades/SRR27160578_FastANI_Output.tsv > ~/hw3/spades/SRR27160578_FastANI_Output_With_Alignment.tsv
cat ~/hw3/spades/SRR27160578_FastANI_Output_With_Alignment.tsv ~/hw3/spades/SRR27160579_FastANI_Output_With_Alignment.tsv ~/hw3/spades/SRR27160580_FastANI_Output_With_Alignment.tsv > ~/hw3/spades/Combined_FastANI_Output_With_Alignment.tsv

#added a header to my tsv file and renamed it to the regex match requested
sed "1i Query\tReference\t%ANI\tNum_Fragments_Mapped\tTotal_Query_Fragments\t%Query_Aligned\tBasepairs_Query_Aligned" ~/hw3/spades/Combined_FastANI_Output_With_Alignment.tsv > ~/hw3/spades/fastani.tsv
mkdir ~/hw3/mlst
cd ~/hw3/mlst

#time for mlst - I included the genomic reference sample in this table as a comparison
conda install -c bioconda mlst -y #installing mlst in this environment, I know it downgraded some previously installed dependencies but that's fine for this specific exercise

#created links to preserve some disk space 
ln -sv ~/hw3/spades/GCF_020735925.1_ASM2073592v1_genomic.fna ~/hw3/mlst
ln -sv ~/hw3/spades/SRR27160578_filtered_contigs.fna ~/hw3/mlst
ln -sv ~/hw3/spades/SRR27160579_filtered_contigs.fna ~/hw3/mlst
ln -sv ~/hw3/spades/SRR27160580_filtered_contigs.fna ~/hw3/mlst
mlst ~/hw3/mlst/*.fna > ~/hw3/mlst/mlst.tsv #ran mlst and created a tsv table
cat ~/hw3/mlst/mlst.tsv
sed -i '1i Sample Name\tSpecies\tST\tAllele1\tAllele2\tAllele3\tAllele4\tAllele5\tAllele6\tAllele7' ~/hw3/mlst/mlst.tsv #header
cat ~/hw3/mlst/mlst.tsv
conda deactivate

#created a new busco environment to assess genome completeness and contamination
conda create -n busco_env -c conda-forge -c bioconda busco -y 
conda activate busco_env
mkdir ~/hw3/busco
cd ~/hw3/busco
busco -i ~/hw3/spades/SRR27160578_filtered_contigs.fna -o SRR27160578_BUSCO -l bacteria_odb12 -m genome --cpu 2 
#after research the best match was the new Bacteria Ortho Database 12 on Busco-I did not find a specific Bordetella parapertussis lineage

cd ~/hw3/busco/SRR27160578_BUSCO #checked out the outputs and saw a summary file that reported data

#summary output reported a 98.3% completeness and 1.8% percent of fragmented and missing buscos which I reported in the contamination column since there was no explicit contamination percentage included in busco reporting
cat ~/hw3/busco/SRR27160578_BUSCO/short_summary.specific.bacteria_odb12.SRR27160578_BUSCO.txt
echo -e "Completion\tContamination" > ~/hw3/busco/quality.tsv
grep "C:" ~/hw3/busco/SRR27160578_BUSCO/short_summary.specific.bacteria_odb12.SRR27160578_BUSCO.txt | \
sed 's/.*C:\([0-9.]*\)%\[.*F:\([0-9.]*\)%\,M:\([0-9.]*\)%.*$/\1\t\2\t\3/' | \
awk '{print $1 "\t" ($2 + $3)}' >> ~/hw3/busco/quality.tsv #parsed and extracted my values and added my F and M values for my contamination values

cat ~/hw3/busco/quality.tsv #checking my quality tsv output
mv ~/hw3/busco/quality.tsv ~/hw3
mv ~/hw3/mlst/mlst.tsv ~/hw3
mv ~/hw3/spades/fastani.tsv ~/hw3 #moving everything to prepare for gunzipping and submission
conda deactivate
