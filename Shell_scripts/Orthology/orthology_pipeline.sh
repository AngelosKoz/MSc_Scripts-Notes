organisms=$1 #This is the file with organism names for the orthology
prot_threshold=${2:-0} #This will be the intended protein length.
first_aas=${3:-\.\*} #This one is used in case we want specific starting amino acids in the protein
#Simple comments
#_# Alternative code for the above approach

echo "=== Add the name of the organisms in a file. For humans use 'homo_sapiens', for mouse use 'mus_musculus' etc ==="
echo "=== All organisms can be found in ensembl_orgs.txt ==="
echo "=== Protein length can be stated as a second variable ==="
echo "=== The first aminoacids for the protein can be stated as a third variable when calling the script ==="
echo "======== Example of calling the script ./orthology_pipeline.sh myfile.txt 150 mmm msa iqtree (In order as seen: file with organism names, 150 is the minimum length of the protein, m would mean we want 3 methionine as the first amino acids, msa and iqtree are methods for tree inference. If skipped, dendroblast and fasttree will be applied. For anything else besides fasttree (iqtree, raxml, raxml-ng) use msa ========"


for check in $( cat ${organisms} ); do
    checking=$(grep -iwc ${check} ensembl_all_orgs.txt)
    if test ${checking} -lt 1
    then
	echo "!!!!=====${check} was not found in the list of organisms of ensembl. Check spelling or use ensembl_all_orgs.txt for verification.=====!!!!"
	break
    fi
done


rm downloaded_releases.list #Incase we run it a second time
touch downloaded_releases.list #Making a file to list the names of all of the downloaded files into

#We can use "grep -c ENS.*000 file" to check if inputs are the same as the original file. (For the first it is simpler with "grep -c \> file

for i in $( cat ${organisms} ); do
    if test ${checking} -lt 1
    then
        echo "!!!!=====${check} was not found in the list of organisms of ensembl. Check spelling or use ensembl_all_orgs.txt for verification.=====!!!!"
	echo "The file downloaded_releases.list can indicate where the issue occured."
	break
    fi
    wget https://ftp.ensembl.org/pub/current_fasta/${i}/pep/CHECKSUMS #Downloading the file with the ftp folder contents
    x=$( grep -iEo "${i}.*\.all\.fa\.gz" CHECKSUMS ) #using grep to find the wanted file and the output of that is saved as variable x
    rm CHECKSUMS
    echo "${x}" >> downloaded_releases.list #This can be used to check which release was downloaded.
    echo "Downloading fasta files for ${i}"
    wget https://ftp.ensembl.org/pub/current_fasta/${i}/pep/${x}  #Passing the grep result in the download link as the variable x
    echo "Unzipping and necessary transformations"
    gunzip -c ${x} > original_${i}.fa
    rm ${x} #Removing the zipped file to make up space
    #We then start making necessary transformations to have all of the aminoacids in the same line, while only keeping the protein name
    sed 's/ .*$/ /' original_${i}.fa | sed 's/>/ >/g' > temp_${i}.fa #The code below works even if we remove ">"
    tr -d '\n' < temp_${i}.fa | sed 's/ /\n/g' | awk 'NR>1' > fixed_${i}.fa #Since the first line will be a blank line, we remove it with awk NR>1
    rm temp_${i}.fa #Removing unecessary temporary file to make up space. We keep the original unzipped incase we want to check something.
    grep -iE -B1 "^${first_aas}[A-Z\*]{${prot_threshold},}[^0-9]$" fixed_${i}.fa | sed 's/--//' > cutoff_${i}.fa #An alternative would be adding [^\.0-9] with dot as an extra failsafe, but since pre-hand transformations may have been done, we keep it with numbers only. Another approach would be with the ENS at the beginning, but it is safer using numbers as ending. We use B1 to get the line before the intended protein (length can be adjusted by user by adding the variable prot_threshold)
    #We need the sed beacause -A and -B in grep leave behind -- in the file.
    #Since we kept ">" we can opt for a grep with [^\>] instead of the ending numbers
    sed -i '/^$/d' cutoff_${i}.fa #This will remove blank lines, -i means in-place
    check_1=$(grep -c ENS.*000 fixed_${i}.fa)
    check_2=$(grep -c \> original_${i}.fa)

    check_3=$(grep -c ENS.*000 cutoff_${i}.fa)
    echo "Original file proteins: ${check_2} | Pre-Cutoff file proteins: ${check_1} | After Cutoff of ${prot_threshold} file proteins: ${check_3}" >> downloaded_releases.list
done

#Another method would be using rsync and wildcard without the need for grep, like so, and then implementing the above transformations:

#_#for i in $( cat ${organisms} ); do
#_#    rsync -havP rsync://ftp.ensembl.org/ensembl/pub/current_fasta/${i}/pep/*\.all\.fa\.gz .
#_#done



mkdir original_fasta_files
mv original_*.fa original_fasta_files/

mkdir pre_cutoff_fasta_files
mv fixed_*.fa pre_cutoff_fasta_files/

mkdir cutoff_fasta_files
mv cutoff*.fa cutoff_fasta_files/

#This is alternative filtering provided by Orthofinder, it only keeps the longest variant of each gene.
#_#for f in *fa ; do python ~/orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done


#This will download orthofinder in case it is not installed
#_#wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz
#_#tar xzvf OrthoFinder.tar.gz #Unpacking


#We could use this approach to set variables 4,5 (method,tree) to something else, incase the user doesnt use the first_ass variable (the first amino acids of the protein). If for example it was ./script myfile threshold first_aas method tree and the user skipped one of those, we could replace it somehow (did not manage to completely fix it, instead i am passing my idea)

#_#nmbr='^[0-9]+$'
#_#if ! [[ ${var} =~ ${nmbr} ]] ; then
#_#    var=$sth_else
#_#else
#_#    va$same
#_#fi




method=${4:-dendroblast} #Options dendroblast and msa

tree=${5:-\ } #Option (iqtree, raxml-ng, fasttree, raxml). msa option is necessary to run something besides fasttree
#We use space here to try and automate it incase no inputs are used for method/tree


#With this i am trying to automate the script in case the user does not use the msa so it take default dendroblast and the fifth variable (tree) will be a space to complement the orthology analysis, since when using -T, the -M has to be msa. In case the user decides to use msa and an approach different than fasttree.
if [[ ${method} == msa ]] ; then
    echo "Using MSA with ${tree}"
    msa_tree=$( echo " -T ${tree}" )
    tree=${msa_tree}
else
    echo "Using the default dendroblast and fasttree inference"
fi

echo "Starting Orthology Analysis"

cd OrthoFinder/
./orthofinder -M ${method}${tree} -t 10 -f ../cutoff_fasta_files/ #-f is for files // -M is for method of inference. dendroblast is the default, we can set it to msa in order to then use -T to change inference tree. I chose to use iqtree since it was faster than raxml. I have this installed tho
#-y : Split paralogous clades below root of a HOG into separate HOGs
#-t : threads

echo " =========== Orthology analysis has finished ==========="
echo " To check results head on to cutoff_fasta_files/OrthoFinder/ folder"
