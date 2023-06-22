H awk diavazei to arxeio grammh grammh apo to standard input. #tail
H awk pairnei to arxeio diavasmeno apo to standard input #xargs prin apo awk

awk '{print}' microarrays100.soft

awk '{print $2}' microarrays100.soft | sort | uniq -c | sort -k1 -g (general number) (k = sthlh)

awk '{print $2"\t"$3"\t"$4}' microarrays100.soft | sort | uniq -c | sort -k1 -g (general number) (k = sthlh) #"\t" gia na diaxwrisw

awk '{print $3/$4"\t"$4}' microarrays100.soft | sort | uniq -c | sort -k1 -g (general number) (k = sthlh) #Division by zero = leksh poy egine 0

tail +2 microarrays100.soft | awk '{ if($3 > 2) print $2"\t"$3}' |  sort -k2 -n/g -r #-n number? r = reverse

tail +2 microarrays100.soft | awk '{ if($3 < $4/2) print $2"\t"$3"\t"$4"\t"$3/$4}' |  sort -k2 -n/g -r #-n number? r = reverse

tail +2 microarrays100.soft | awk '{ if($3 < $4/2) print $0}' |  sort -k2 -n/g -r #To $0 typwnei olh th sthlh

print NF #number of fields (sth periptwsh mas einai column).Tha doyme dhladh posa "pragmata" yparxoyn anamesa sta tabs.
#Deixnei thn TELEUTAIA sthlh 

print $(NF-1) #me to dolario pairnw th timh ths NF. NF-1 = proteleutaia

awk '{ for(i=1; i<=3; i++) printf $i"\t"; print "$4"}' #H printf tha typwsei #apo th prwth ews thn 3h sthlh, kai sto telos tha typwsei thn sthlh 4 kai tha dwsei new-line, dhladh tha paei apo katw

awk '{ for(i=1; i<=3; i++) printf $i"\t"; print ""}' | grep "$2 " #ama den dhlwsw kati sto print, tha valei tab endiamesa giayto kanw grep me keno

awk '{print NR"\t"$2}' #NR = dinei ton arithmo ths ekastote grammhs. $NR = diavazei prwth grammh typwnei prwto pedio .... ara telika typwnei prwth sthlh

awk '/C4orf/ {print $3}' #psaxnei pantoy to gonidio C4orf kai to typwnei

awk -F,    // awk -F\t #Toy deixnw me ti na diavazei to arxeio (","  "\t")

x=2000; cat... | awk -v var=$x -F"," '{print $3+var}' #me to v ftiaxnw metablhth poy onomazetai var, kai ths dinw thn timh x

