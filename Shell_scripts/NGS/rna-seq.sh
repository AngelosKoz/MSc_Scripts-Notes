id_file=$1 #Uses a file with the individuals of interest so as to download the files 
organism=$2 #Uses the species of the individuals (e.x Homo_sapiens)


#_# --> Means this is an alternative code
# --> Simple Comments

#Get File with all the information. I found this as a general file and wanted to work on it to get the necessary information.

wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-197/E-MTAB-197.sdrf.txt

#The file seems to not print columns properly, so we replace spaces with "_" on the first row. Tried doing it with 2 different methods, one was replacing everything in the file with sed as seen below

#_#sed 's/ /_/g' E-MTAB-197.sdrf.txt > e-mtab_fixed.txt

#Or by using awk. Went with this one further on.

awk '{gsub (/ /, "_", $0); print $0}' E-MTAB-197.sdrf.txt > awk_sub_e-mtab.txt

#Doing a whole file replace does not seem to affect the url, or the patient ID which i am interested in, located in column 31 and 1 respectively.


#The second which gives me some weird results was trying to only replace the first row spaces with "_" but it gave me some inputs where i wasnt expecting ones when printing columns with awk. I did not opt for this method as seen below

#_#head -1 E-MTAB-197.sdrf.txt | sed 's/ /_/g' > header.txt
#_#tail -160 E-MTAB-197.sdrf.txt > e-mtab-tail.txt
#_#cat header.txt e-mtab-tail.txt > e-mtab_test.txt

#Since there are many parameters to try and automate this, nor parameters distinctiv to the 5 patients i had to create a list with the patient IDs of interest.
#I wanted to use awk for a change in the for loop and found this method. I checked that column 31 contained the urls and used it to extract them for the download.
#The pat_id.txt was created using nano which i suppose would be a more automated way since the user can create the file for different patients.
#One way to approach it by passing a variable would be the below, where the user would indicate the file name containing the patient id's of interest.


for i in $( cat $1 ); do
    awk -v var="$i" '$0 ~ var' awk_sub_e-mtab.txt | awk '{print $30}' >> patient_err_id.txt
    awk -v var="$i" '$0 ~ var' awk_sub_e-mtab.txt | awk '{print $31}' > patient_url.txt
    for url in $( cat patient_url.txt ); do
	echo "Downloading ${i} fastq"
	wget ${url}
    done
done

rm patient_url.txt #Deleting the url file
uniq patient_err_id.txt > patient_err_uniq_id.txt
rm patient_err_id.txt

#Another way to approach it would be by using grep on the result of the awk like so:
#_# awk -v var="$i" '$0 ~ var' awk_sub_e-mtab.txt | grep -o "ftp://.*fastq.gz" > patient_url.txt
#This way there is no need to check the column, we would rather have to check the start and end of the provided link.

echo "All files have been downloaded. Now unzipping"
gunzip *.fastq.gz #Unzip the files

#This step is for fast QC
#The patient_err_uniq_id.txt contained characteristic id's for the patients.
#Making a few directories so the files are all tidy
mkdir fasta_qc
mkdir fastq_untrimmed_pat_files
mkdir fasta_qc/sickle_single
mkdir fasta_qc/untrimmed_fastqc
mkdir fasta_qc/trimmed_fastqc


#Trimming here is done with the suggested parameters 12 and 15 (since this is paired end). The most proper way would be to check the qc files and go on based on those.
for j in $( cat patient_err_uniq_id.txt ); do
    echo "Running Quality Control on ${j}"
    fastqc ${j}*.fastq
    mv ${j}*fastqc* fasta_qc/untrimmed_fastqc/
    echo "Trimming ${j}"
    sickle pe -f ${j}_1.fastq -r ${j}_2.fastq -t sanger -o trimmed_${j}_1.fastq -p trimmed_${j}_2.fastq -s ${j}_sickle_single.fastq -q 12 -l 15
    echo "Running verification Quality Control on trimmed files"
    fastqc trimmed*
    mv trimmed*fastqc* fasta_qc/trimmed_fastqc/
    mv *_sickle_single* fasta_qc/sickle_single/
    mv ${j}_1.fastq fastq_untrimmed_pat_files/
    mv ${j}_2.fastq fastq_untrimmed_pat_files/
done

cd fasta_qc/untrimmed_fastqc
echo "Running multiQC on original files"
multiqc . #Multi qc on untrimmed for comparison

cd ../trimmed_fastqc/ #Moving into trimmed_fastqc dir
echo "Running multiQC on trimmed files"
multiqc .

cd ../../ #Moving to original dir


#One way of downloading ref genome: This does not automate things much and i did not opt for this.
#_#wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
#_#wget wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz


#Instructions for this on https://www.metagenomics.wiki/tools/fastq/ncbi-ftp-genome-download
#The way this is structured, both taxonomy id and organism name can be used as arguments

echo "Downloading genome index file"
rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt ./

tail -289715 assembly_summary_refseq.txt | sed 's/# //1' | awk 'BEGIN{FS=OFS="\t"} {gsub(/ /, "_", $8)} 1' > assembly_fixed.txt #This will remove the first row that is instructions, fix the first row so it will be tab delimited and replace spaces of the organism names with _ so they can be used as argument from bash


awk -F"\t" '{print $7"\t"$8}' assembly_fixed.txt > organism_name_id.txt #This file can be used to search for an organism id when encountering issues with the argument

grep -iE "${2}" assembly_fixed.txt | cut -f 20 > ftp.txt

head -1 ftp.txt > ftp_first.txt

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_first.txt | sed 's/https/rsync/g' > fna_link.txt

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_first.txt | sed 's/https/rsync/g' > gff_link.txt

echo "Downloading ${2} Fasta file"
source fna_link.txt
echo "Finished Downloading ${2} Fasta file"

echo "Downloading ${2} Annotation Genome"
source gff_link.txt
echo"Finished Downloading ${2} Annotation Genome"

echo "Unpacking files"
gunzip -v *.gz #Unpack everything together. This is optional for bwa but not for genome. -v shows progress
echo "Finished Unpacking"



mv *.fna ${2}.fna #Rename to simpler terms
mv *.gff ${2}.gff


#Moving on to star indexing
mkdir star_ind
cd star_ind/ #I am now at Dir: <starting dir>/star_align
#_#ln -s ../${2}.gtf . #This is needed if we want to run with gtf format
ln -s ../${2}.gff .
ln -s ../${2}.fna .

#Indexing with gtf converted file
#Using this method to automate the process, we would need the gffread tool to convert it to gtf so that STAR and HTSeq can use it.
#_#gffread *.gff -T -o ${2}.gtf #Self note: try to find an executable to pass along with the pipeline
#_#STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ${2}.fna --sjdbGTFfile ${2}.gtf --sjdbOverhang 100 --runThreadN 12

#Indexing with gff
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ${2}.fna --sjdbGTFfile ${2}.gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --runThreadN 12

fastqdir="../"
for map in $( cat ../patient_err_uniq_id.txt ); do
    echo "Mapping ${map}"
    STAR --genomeDir ./ --readFilesIn ${fastqdir}trimmed_${map}_1.fastq ${fastqdir}trimmed_${map}_2.fastq --runThreadN 12 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${map}
done

echo "Finished mapping all files"


###### Here i am trying to get counts from the STAR itself without using HTSeq, by adding the option --quantMode GeneCounts. I get this error:
#error#./rna-seq.sh: line 161:  1113 Segmentation fault      STAR --genomeDir ./ --readFilesIn ${fastqdir}trimmed_${map_2}_1.fastq ${fastqdir}trimmed_${map_2}_2.fastq --runThreadN 8 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix star_count_${map_2}

#mkdir star_with_counts
#cd star_with_counts/
#ln -s ../${2}.gff .
#ln -s ../${2}.fna .
#STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ${2}.fna --sjdbGTFfile ${2}.gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --runThreadN 12

#fastqdir="../"
#for map_2 in $( cat ../patient_err_uniq_id.txt ); do
#    echo "Mapping ${map_2}"
#    STAR --genomeDir ./ --readFilesIn ${fastqdir}trimmed_${map_2}_1.fastq ${fastqdir}trimmed_${map_2}_2.fastq --runThreadN 8 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix star_count_${map_2}
#done  



echo "Use samtools to view the BAM file like so: samtools view ERR009159Aligned.sortedByCoord.out.bam | less -S"

cd ../ #I am now at Dir : <starting dir>/
mkdir htseq
cd htseq/
ln -s ../${2}.gff .
ln -s ../star_ind/*Aligned*.bam .


#Sam tools indexing necessary before HTSeq
for ind in $( ls *Aligned*.bam ); do
    echo "Indexing ${ind}"
    samtools index ${ind}
done

echo "Finished Indexing"
    
#Alternative way for htseq with xargs. Does all 5 parallel.
#_#ls *Aligned*.bam | xargs -n1 -P5 -I {} sh -c 'htseq-count -s no -r pos -t exon -i pacid -f bam "$1" ./*.gff > "$1.counts"' -- {}


#Htseq for counts. -s is for strand specific (no) -r is for sorted by (position) 
for bm in $( ls *Aligned*.bam ); do
    echo "Running HTSeq on ${bm}"
    htseq-count -s no -r pos -t exon -i gene_id -f bam ${bm} ./${2}.gtf > ${bm}.counts #change this to gff / gtf accordingly. Gff does not work so far, try to do filtering by removing some rows
done

echo "Finished generating counts. Files are generated with ending of .counts"

#Trying to get counts from gff by Dbxref, after filtering with: grep -iE "gene=" Homo_sapiens.gff > fixed_${2}.gff. This was necessary since some of the rows did not include input for HTSeq. the gene= pattern was not seen in those lines. This can be checked with grep -iEv "gene=" Homo_sapiens.gff
#_#for bm in $( ls *Aligned*.bam ); do
#_#    echo 'Running HTSeq on ${bm}' #careful of ' instead of double
#_#    htseq-count -s no -r pos -t exon -i Dbxref -f bam ${bm} ./fixed_Homo_sapiens.gff > new_test_${bm}.counts
#_#done


