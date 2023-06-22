#x=NM_001276351
x=$1
tr=$2

mkdir ucsc

##https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz20way/alignments/refGene.exonNuc.fa.gz

#rsync -avz --progress "rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz20way/alignments/refGene.exonNuc.fa.gz" ./ucsc

cd ucsc

#gunzip -k refGene.exonNuc.fa.gz # me to k an yparxei den tha to ksana unziparei

#grep -E -A1 "^>NM_001276351_.*_6_6" refGene.exonNuc.fa | sed 's/-//g' | sed 's/_6.*//' > NM_001276351.fa #Sed antikatastash

grep -E -A1 "^>${x}_.*${tr}" refGene.exonNuc.fa | sed 's/-//g' | sed 's/_6.*//' > ${x}${tr}.fa #Sed antikatastash

pwd

mafft ${x}${tr}.fa > ${x}${tr}.aln

raxml-ng --msa ./${x}${tr}.aln --model GTR+G

