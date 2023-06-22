#`` #Toys ypoloipoys na mathw
# {}
# $() #To apotelesma ths entolhs pou einai mesa sth parenthesh tha paei san value sto i
#(()) epitrepoyn arithmitikes prakseis dhladh kanw a=1 | b=2 | d=$((a+b)) | echo $d | d =$((d+1)) | d=$((++d))

# a = $((a++)), prwta to a tha paei sto neo a kai meta tha aykshthei kata 1, alliws antistrofa



#echo -n > fileALL.txt #typwnei adeia grammh sto file all, an yparxei tha to katastrepsei kai tha typwsei kenh grammh // ftiaxnei arxeio xwris tipota kai vazei epanalhptika


for i in {1..100} #For, onoma metablhthw in, sth seira poy einai, prepei na metrhsw to teleytaio
do
    cat file${i}.txt
done > fileALL.txt


#Allos tropos

x=$1 ############## x=$@ --> aperioristes metablhtes

for i in {1..100} 
do
    cat file${i}.txt 
done > ${x} # Mpainoyn gia prostasia metablhths


-----------------------------


for i in {1..100..2} #to ..2 einai o bhmatismos, dhladh ana posa 
do
    echo $i #Typwnei 

done

----------------------------

#Thelw na kanw alignment file0.fa file2.fa ... file100.fa

for i in {0..100..2}
do
    mafft files${i}.fa > file${i}.aln
    raxml --msa file${i}.aln
done

#Seq 1 10 paragei tous arithmous 1 ews 10 me th seira

#`` The ektelestei entolh sto unix, kai to apotelesma kathe fora tha apothikeytei sthn allh. Mporw na xrhsimopoihsw to apotelesma
for i in `seq 1 2 10` #Edw o bhmatismos einai anamesa. dhladh tha typwsei 1,3,5,7,9
do
    echo file$i
done

----------------------------------

for i in $(seq 1 2 10);do #PROTEINOMENOS TROPOS
    echo example${i}.fa;
    grep ACACA example${i}.fa; done
-----------------------------------

for i in $( cat patterns.txt ); do
    grep $i example.txt; done #Tha vrei thn allhloyxia poy yparxei mesa sto patterns.txt sto arxeio example.txt


----------------------------------

#Thelw na kanw rename ola ta file*.txt se file*.txt.txt

for i in $(ls file*.txt) #Mporw na kanw kai file*
do
    #mv $i ${i}.txt #Den vazw to file prin, giati to ls epistrefei pisw kati pou exei mesa to file
    #ayto tha metonomasei
    y = $(echo $i | sed 's/\.txt/\.TXT/')
    mv $i $y
done

----------------------------------

for ((i=1; i<101; )); do
    echo $i ; done

    #i-- / --i tha dwsei aenah loupa
      #Orizw tis "synthikes" kai trexoyn mono an isxyei h deyterh

---------------------------------------
### If Conditional ###

if test 3 -lt 5
then
    echo '3 is less than 5'
fi #(me to fi kleinei)

--------------------------------------

if test -e file28.txt #Tha elegxei an yparxei file poy onomazetai etsi mesa sto fakelo
   #-f elegxei an yparxei to arxeio (file)
   #-d dir an yparxei o fakelos
   #-e to kanei se directory + arxeio
then
    mv file28.txt file28.TXT
fi

-------------------------------------

if test -d dir1
then
    if test -f file28.TXT; then
	mv file28.TXT dir1
    fi
else
    mkdir dir1
    if test -f file28.TXT; then
	mv file28.TXT dir1
    fi
fi


-------------------------------------

if [ -d dir1 ] #To idio me if test -d dir1... Opwsdhpote me keno sthn arxh/telos giati einai ENTOLH
then
    if test -f file28.TXT; then
	mv file28.TXT dir1
    fi
elif test -d dir2; then
    mv file28.TXT dir2
fi

-------------------------------------

if (( 1 < ))
then
    echo "A"
else
    echo "B"
fi

--------------------------------------

a = 1
b = 2
if (( $a+$b < 5 )) #Extended  version praksewn ///////// Antistoixa paei kai to [[ ]] 
   then
    echo "A"
else
    echo "B"
fi

-----------------------------------

a = 1
b = 2
if [ $a -lt 5 -a 1-lt 4 ] #To -a shmainei END
then
    echo "A"
else
    echo "B"
fi


