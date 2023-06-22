# NGS - RNASeq Analysis ------ Gia to idio eidos gonidiwma

#Stoxos na brw th pososthta poy paragetai apo ena gonidio, kai oxi to gonidio kathauto
#Mporw na prosdiorisw thn posothta ekfrashw

#Ginetai cDNA kai sth synexeia allhloyxish --> Mapping (


#BWA -- Gia allhloyxish DNA

#Epeidh sugrinw me idioy mhkous allhloyxies, den xreiazetai na kanw kanonikopoihsh


# 1) SRR Files - Download. Auta uparxoyn apla gia na ta metatrepsw se fastq. 
#sra-->fastq. Einai binary

#Apo NCBI einai kalo na xrhsimopoihsw to prefetch poy dinoyn gia katebasma twn arxeiwn

#--# cat srr.list | xargs -n1 -P6 -I {} curl "https://sra-pub-run-odp.s3.amazonaws.com/sra/{}/{}" -o {} .sra
#xargs xrhsimopoiei to output cat sth sunexeia.
# -n1 == kathe fora na xrhsimopoihsei 1 th fora gia na mhn ta mperdepsei
# -P6 == Tha ta xrhsimopoihsei se mia 6ada
# -I == gia na dwsei ta arguements
# -o einai ths curl gia na mporesei na moy dwsei to arxeio

# 2) Metatroph se fastq

# 3) Trimming -- Sickle -- Diabazei poiothta diabasmatos. Kobei apo ta akra an to akro einai mikrotero apo kapoio threshold. An einai megalo auto pou prepei na kopsei. Exei epishs kapoio threshold gia to megethos poy tha kopsei
#cat srr.list | xargs -n1 -P6 -I{} sickle se -f $dir{}.fastq -t sanger -o ${dir}trimmed_{}.fastq -q 35 -l 45 (to q kai to l einai ta 2 threshold)


# 4) QC (quality control) fastqc -- Pairnei apotelesma apo fastqc (qc per file) kai ta kanei na fainontai omorfa.
#Sth synexeia an thelw ta apotelesmata gia kathe ena apo auta na ta exv mazi tote kanw --multiqc--
#Sta akra synithws xanw poiothta 

#   -- MultiQC -- : pip install multi

# 5)Mapping (STAR) 
# BAM (binary of SAM - pou exei antistoixhthei ena gonidiwma) ή SAM files --> Antistoixh kathe read sth perioxh toy. Kanei to allignment - mapping

# 6)BAM --> counts (annotation) Dhladh posa reads exw se kathe perioxh  (HTSeq-count)

# 7) R Differential expression analysis

#import subprocess, argparse, sys

parser = argparse.ArgumentParser(description = '')
parder.add_arguement()
parder.add_arguement()
parder.add_arguement()

args = parser.parse_args()


#prefetch SR142019

#fasterqdump --> Metatroph arxeiwn SRA se fastq


#Single End Read: Kanw seq mono th mia akrh

#Paired End Read: Kanw seq kai stis dyo akres toy tmhmatos poy tha allhloyxisw. Briskei Alternative splicing




#ln -s == Softlink
#rename 's/.fastq.small.fastq/.small.fastq/' *.fastq

#For tmux if there is an issue
echo #DISPLAY --> DEIXNEI LOCAL HOST
DISPLAY=localhost:11.0


#cat srr.list | xargs -n1 -P6 -I{} fastqc {}.small.fastq
#Parallhla kai gia ta 6 kanei to fastq

##### PHYTOZOME #######

#rna-star --- gia na to katebasw
#STAR -- trexw etsi

# prepei na ginei to indexing prwtoy treksw to STAR
