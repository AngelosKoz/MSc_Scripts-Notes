1)grep -iEc "^hCoV-19*\-+A" brazil.FA
2)grep -n -o -E "^-*T" brazil.FA | wc
3)grep -o 'ACCAC[ACGT]' brazil.FA | sort | uniq -c
4)grep -Eo '^*[AGCT]' brazil.FA | grep Eo '[^]' | sort | uniq -c



2)grep -vc "\#" data.vcf    ///////    grep '#' data.vcf
3)grep -w A10000 "11815" data.vcf ////// grep -v 11815 | wc -l    ///////  grep -w 11815 data.vcf
///////  cat -n data.vcf | grep -w "11815"
4)grep '0/1' data.vcf /////// grep -E '1/0|0/1' ///// grep -o "./." data.vcf | sort | uniq -c
5)grep -v '#'  | head -100 | tail -1
6)grep '1/1' data.vcf | wc -l     -->    sed 's/1\/1/0\/0/g' | grep '0/0' | wc -l

grep -E -A1 "^>NM_001276351_.*_6_6" refGene.exonNuc.fa | sed 's/-//g' | sed 's/_6.*//' > NM_001276351.fa
grep -E -A1 "^>${x}_.*${tr}" refGene.exonNuc.fa | sed 's/-//g' | sed 's/_6.*//' > ${x}${tr}.fa 


for i in $( cat patterns.txt ); do
    grep $i example.txt; done #Tha vrei thn allhloyxia poy yparxei mesa sto patterns.txt sto arxeio example.txt








#Telestes

“grep” // Λειτουργει με το Control + Z, fg
- -color --> Χρωματιζει αυτο που θα ψαξω
-n --> δινει γραμμες ολοκληρες
-c --> μετραει το πληθος γραμμων που υπαρχει η λεξη, δινει αριθμο
-ο --> οι ευρεσεις του μοτιβου στο αρχειο, δινει το μοτιβο
-w --> βρισκει την λεξη ΑΚΡΙΒΩΣ οπως την εγραψα, και οχι στο περιπου
-v --> H γραμμη που δεν περιεχει το pattern
-l --> Δινει ονομα αρχειου
-E --> Extended version
-AN --> Γραμμη που περιεχει αυτο που ψαχνω +N γραμμες μετα (αρα Α1 --> γραμμη + επομενη)
"[]" -->  ψαχνει κάποιο απο τα γραμματα που μπορει να βρισκονται εδω π.χ grep -i “tgg[ACGT]gg”. [A-Z] απο A εως Z
“^” --> Βρισκει οτιδηποτε βρισκεται μετα απο αυτο
“[^]” -->  τοτε ΔΕΝ βρισκει αυτο που βρισκεται μεσα στο “[^]”
“\” --> escape characters
“A$” --> σημαινει να τελειωνει σε Α
‘^s.*A$’ --> Να ξεκιναει απο s να τελειωνει σε Α και οτιδηποτε αλλο υπαρχει ενδιαμεσα 
( + --> 1 η περισσοτερες φορες // * --> 0 ή περισσοτερες φορες // ? --> 0 ή 1 φορες)
{,8} --> το πολυ 8 χαρακτηρες, μπαινει μετα απο το γραμμα που θελω -->AG{,8}TTA (το πολυ 8 G)
{2,8} --> απο 2 εως 8 	// 	{8} --> ακριβως 8
"." --> οποιοδηποτε γραμμα ή αριθμος
"\<" --> Μπαινει στην αρχη, ψαχνει την λεξη που ξεκιναει απο αυτο που ψαχνω
"\>" --> Μπαινει στο τελος, ψαχνει τη λεξη που τελειωνει σε αυτο που ψαχνω   
“\|” προσθετει δευτερη εντολη (pipe): To output της προηγουμενης εντολης χρησιμποιειται ως input στην επομενη εντολη. Και στο grep