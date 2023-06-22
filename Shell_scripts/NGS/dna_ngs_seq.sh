

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.analysis.sequence.index

#_# --> Alternative code
# --> Simple comments

pat_id=$1 #This uses a file with the names of individuals of interest
organism=$2 #This uses the species of the individual (e.x Homo_sapiens)

mv *.index full_seq.index

#_#awk -F"\t" '{print $1"\t"$10"\t"$20"\t"$24"\t"$25}' full_seq.index | awk -F"\t" '$3 != "" {print}' > fixed_seq.index #Alternative way keeping more columns



#We keep columns 20 (link for download) and 25 (BP counts) and we remove non-paired rows with the pipe
awk -F"\t" '{print $3"\t"$20"\t"$25}' full_seq.index | awk -F"\t" '$2 != "" {print}' > fixed_seq.index


for i in $( cat ${1} ); do
    awk -v var="$i" '$0 ~ var' fixed_seq.index | sort -k 3 -n -r | head -2 | awk '{print $2}' > ${i}_highest_bp_link.txt
    for j in $( cat ${i}_highest_bp_link.txt ); do
	echo "Downloading ${i} files"
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/${j}
	echo "Finished Downloading ${i} files"
    done
done

#Inserting the genome download incase these want to be ran seperately

echo "Downloading genome index file"
rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt ./

tail -289715 assembly_summary_refseq.txt | sed 's/# //1' | awk 'BEGIN{FS=OFS="\t"} {gsub(/ /, "_", $8)} 1' > assembly_fixed.txt #This will remove the first row that is instructions, fix the first row so it will be tab delimited and replace spaces of the organism names with _ so they can be used as argument from bash


awk -F"\t" '{print $7"\t"$8}' assembly_fixed.txt > organism_name_id.txt #This file can be used to search for an organism id when encountering issues with the argument

grep -iE "${2}" assembly_fixed.txt | cut -f 20 > ftp.txt

head -1 ftp.txt > ftp_first.txt

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_first.txt | sed 's/https/rsync/g' > fna_link.txt

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_first.txt | sed 's/https/rsync/g' > gff_link.txt #this gff file might not be necessary for the bwa

echo "Downloading ${2} Fasta file"
source fna_link.txt
echo "Finished Downloading ${2} Fasta file"

echo "Downloading ${2} Annotation Genome"
source gff_link.txt
echo"Finished Downloading ${2} Annotation Genome"

mkdir all_dl_links
mv *_link.txt all_dl_links/

echo "Unpacking files"
gunzip -v *.gz #Unpack everything together. This is optional for bwa but not for genome
echo "Finished Unpacking"

mv *.fna ${2}.fna #Rename to simpler terms
mv *.gff ${2}.gff


#Generate the unique experiment ID (SRR / ERR)
ls *.fastq | sed 's/_[12].*//g' | uniq > pat_uniq_rr_id.txt


#Genome indexing with bwa
echo "Starting genome indexing"
bwa index -a bwtsw ${2}.fna #Bwtsw algorithm works for human genome. -a i is another algorithm that works for smaller

#Μapping with bwa. -t 8 = 8 core, then we use the ref genome, the paired end files and the -o output
for pat in $( cat pat_uniq_rr_id.txt ); do
    echo "Generating SAM files for ${pat}"
    bwa mem -t 8 ./${2}.fna ${pat}_1*.fastq ${pat}_2*.fastq -o ${pat}.sam
    #_#touch ${pat}.bwa.end
done
echo "Done generating SAM files"


for fx in $( cat pat_uniq_rr_id.txt ); do
    samtools fixmate -O bam ${fx}.sam ${fx}.bam #Convert sam to bam
    samtools sort -o ${fx}.sort.bam ${fx}.bam #Sorts the files
    echo "Starting base filtering with mpileup for ${fx}, this may take a while"
    bcftools mpileup -g 10 -Oz -o ${fx}.gvcf.gz -f ${2}.fna ${fx}.sort.bam #Mpileup is base filtering. Maybe more is necessary?
    echo "Base filtering done. Starting indexing (output is .csi file)"
    bcftools index ${fx}.gvcf.gz #Indexes vcf: Output .csi file
done


echo "Almost there. Now merging files, this may take some time."

bcftools merge -Oz --gvcf ${2}.fna --merge all -o merged.vcf.gz *.gvcf.gz #Merges the files


bcftools call -mv merged.vcf.gz -o merged.vcf.gz.mv.vcf 


#Tidying up the directory
mkdir index_files
mv *.index index_files/

mkdir original_fastq
mv *.fastq original_fastq/

mkdir sam_files
mv *.sam sam_files/

mkdir bam_files
mv *.bam bam_files/

echo " ===== All done. Final file is: merged.vcf.gz.mv.vcf ====== "

#Since giannis asked this, what other filtering should be done besides mpileup? Also should it be done pre or post vcf?
#Found this page but im not sure if this appropriate : https://www.htslib.org/workflow/filter.html
