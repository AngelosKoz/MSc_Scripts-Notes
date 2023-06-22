#Askhsh 1: mkdir ex3

#Askhsh 2 + 3:
#a = 0
#for i in {1..101..3}; do
#    touch f${i}.txt
#    echo (( ++a )) > f${i}.txt 
#done

#Tropos paylidh
a=0
for i in {1..100..3}; do
    a=$((++a)) #diaforetiko apo a++, tha dwsei 0
    #echo $((++a)) #dinei 2,4,6,8..
    echo ${a} > f${i}.txt
done

#Askhsh 4

for i in *.txt; do
    grep -iEl "1" ${i}; done

#Askhsh 5

for i in {22..37..3}; do
	 if [ -e f${i}.txt ]
	 then
	     echo "Yparxei"
	 else
	     echo "Den Yparxei"; fi
done




