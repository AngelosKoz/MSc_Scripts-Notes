#To find psaxnei opoudhpote se upofakelous apo edw pou eimai ena arxeio
find DIRECTORY -name(flags) '*.txt'

./ -iname '*.txt' | xargs ls -s #-i (case insensitive)//-name (onoma) , psaxnei arxeia kai pernane diadoxika san teleytaio orisma sto ls //
#find, se poio fakelo, ti motivo.
#Xwris ta autakia tha doylepsei mono an den exw arxeio me .txt sto fakelo poy eimai

#xargs pairnei to apotelesma kai to vazei san orisma

for i in $(find ./ -iname '*.txt' ); do
    ls -s $i; done

for p in $(cat patterns.txt patterns2.txt); do grep -o $p example10.fa;done
#To idio me#
cat patterns.txt | xargs -I {} grep -o {} example10.fa | uniq -c


echo "example10.fa,example6.fa" | xargs -d, | xargs ls -s #shmainei oti diavazei san "," to diaxwristiko sta input (exa.fa , exa2.fa)

find -iname '*.txt' -print0 | xargs -0 ls -s #Se periptwsh poy exw mesa se ena arxeio keno sto onoma, ayto sto apofeygei

echo "ATC" | xargs.....


-----------------------------



#echo $PATH #mas dinei tous fakelous

#emacs -nw ~/.bashrc

#sto telos toy arxeioy vazw
#HOME : /home/student/

#PATH=$HOME/bin/:$PATH #Auto to directory KAI (:)auto to path otidhpote eixe mesa
#export PATH

#source ~/.bashrc #trexei tis metavlhtes ths bashrc
#to idio me
#. ~/.bashrch



ln -s $HOME/..... ./ #(bin)
ln -s /home/student/linuxclass2022/pavlos_pavlidis/lecture05/test/software/raxml-ng/bin/raxml-ng ./
#prepei prwta na paw ston bin (ekei poy mporw na treksw) kai apo ekei na kanw softlink
