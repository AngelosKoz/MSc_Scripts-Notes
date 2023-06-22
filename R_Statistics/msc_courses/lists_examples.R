---
title: "Lists"
output:
  html_document:
    df_print: paged
  html_notebook: default
  pdf_document: default
  word_document: default
---
## Lists
Οι λιστες (list) ειναι μια από τις ισχυρότερες δομές δεδομένων στην R. 
- είναι γενικές, δηλαδή μπορεί να περιέχουν δεδομένα διαφορετικού τύπου (κάτι που δεν έχουν τα vectors και matrices)
- οι καταχωρήσεις μπορεί να έχουν όνομα (ισχύει και στα vectors αυτό)
- οι καταχωρήσεις μπορεί να έχουν διαφορετικό μήκος (κάτι που δεν ισχύει στα data.frames)

Άρα μπορούμε να πούμε ότι οι list ειναι δομές δεδομένων που επιτρέπουν την αποθήκευση διαφορετικού τύπου μεταβλητών, οι καταχωρήσεις ειναι ονοματισμένες και η κάθε καταχώρηση μπορεί να έχει διαφορετικό μήκος. 

# Κατασκευή
```{r}
a = list() # άδεια λίστα. H a ειναι μια λίστα που δεν περιέχει τίποτα

a = list(x=c(1,2,3)) # Η a ειναι μια λίστα που περιέχει ένα αντικείμενο, το x, που ειναι vector από 3 στοιχεία

a = list(x = c(1,2,3), y=c("A")) # Η a περιέχει 2 στοιχεία, το x και y, που το καθένα ειναι vector
## από αριθμούς και χαρακτήρες αντίστοιχα. 

a = list( x = list(x=c("A"), z=c(1,2,3)), y = "kalimera") # Η a είναι μια λίστα , που περιέχει μια λίστα (την x), και ένα character-string το y. Η λιστα x, περιέχει δύο vectors

a = list(c(1,2,3)) # λίστα που περιέχει ένα vector, μη ονοματισμένο (χωρίς όνομα)
```

Βλέπετε από τα παραπάνω παραδείγματα πόσο ισχυρές ειναι οι λίστες. 

# Πρόσβαση στα στοιχεία της λίστας
Στις λίστες έχουμε δύο operators (δύο τρόπους) να έχουμε πρόσβαση στα στοιχεία της. Το `[[ ]]` και το `[]`. 
Θεωρείστε το παρακάτω παράδειγμα.

```{r}
x = list(a=c(1,2,3), b=c("A", 'b'))
```

Tότε, το πρώτο στοιχείο της λίστας x, ειναι το a, και το δευτερο στοιχείο της λίστας x είναι το b. 

Έχουμε λοιπόν:
```{r}
## first element
## Πρόσβαση με το όνομα του στοιχείου
x$a ## το στοιχείο του x, που έχει το όνομα a

## Πρόσβαση με indexes

## Ώς vector
x[[1]] ## πρώτο στοιχείο ειναι το c(1,2,3) ΕΙΝΑΙ vector

## Ως λίστα
x[1]   ## πρώτο στοιχείο ΕΙΝΑΙ ΛΙΣΤΑ
```

ΠΡΟΣΟΧΗ! 
Όταν έχουμε πρόσβαση ως λιστα δηλαδή με το `[ ]` operator, τότε μπορούμε να κάνουμε slice, δηλαδή να πάρουμε και πανω από ένα στοιχεία. Για παράδειγμα

```{r}
a = list(x=c(1,2), y="a", z=c(5,4,5,3))

## μία λίστα με τα 2 πρώτα στοιχεία
a[1:2]

## μία λίστα με τα στοιχεία 1,3
a[c(1,3)]
```

ΠΡΟΣΟΧΗ!!! Το `a[[1:2]]` ίσως δεν βγάλει λάθος, όμως ειναι αποτέλεσμα χωρίς νόημα!


```{r}
a[[c(1,2)]]
```

Προσέξτε τις διαφορές:
```{r}
a = list(x=c(1,2,3), y=c("a", "b"))
## Το πρώτο στοιχείο (x)
a$x

## Επίσης το πρώτο στοιχείο
a[["x"]] ## θα ειναι vector

## Επίσης το πρώτο στοιχείο 
a["x"] # Θα ειναι λίστα
```

Παρατηρήστε την διαφορά

* `a$x`
* `a[["x"]]`

(χωρίς και με εισαγωγικά)

To γεγονός ότι έχουμε πρόσβαση με τα ονόματα, μας επιτρέπει να τις χρησιμοοποιήσουμε ως associative arrays (hashes)

```{r}
a=list(antonis=29, vaggelis=54, kostas=33)
# όπου σε κάθε όνομα έχουμε αποθηκευσει την ηλικία
# Ας υποθέσουμε ότι η λίστα έίχε εκατομμύρια καταχωρήσεις (αντί για μόνο 3)
```

Στην C δεν υπάρχει η δυαντότητα αυτή. Εκεί για να είχαμε τη σχέση μεταξύ ονόματος και ηλικίας θα έπρεπε:
```{r}
nam = c("antonis", "vaggelis", "kostas")
age = c(29, 54, 33)

# όμως για να βρούμε τι ηληκία έχει ο Κώστας, θα πρέπει να διατρέχουμε όλη τη λίστα
# Στην R απλά θα γράφαμε
a$kostas
```

Ένα δεύτερο συνηθισμένο παράδειγμα αφορά στην έκφραση γονιδίων. Σε έναν οργανισμό Α, που έχει 20,000 γονίδια, τα γονίδια a1, a2, a7 μας ενδιαφέρουν επειδή σχετίζονται με μία ασθένεια. Οι εκφράσεις γονιδίων στον Α, υπάρχουν σε μια λίστα `ΕΑ`. Για παράδειγμα, `EA[["a1"]]`, δίδει την έκφραση του a1, στον Α.  Αν θέλουμε να δούμε σε έναν οργανισμό Β, την έκφραση των τριων αυτών γονιδίων a1, a2, a7, θα πρέπει να διατρέξουμε όλη την λίστα των γονιδίων `ΕΒ` του B ώστε να βρούμε ποιες τιμές έκφρασης αντιστοιχούν στα γονίδια a1, a2, a7. Αν όμως έχουμε ονοματίσει την λίστα γονιδίων, τότε μπορούμε πολύ εύκολα να ζητήσουμε τις εκφράσεις ως εξής `ΕΒ[["a1"]]  ΕΒ[["a2"]]  ΕΒ[["a7"]]`. 

*Γενικά οι ονοματισμένες λίστες ειναι πολύ χρήσιμες στην βιοπληροφορική*

## Ονοματισμένα vectors
Και τα vectors μπορούν να φέρουν ονόματα στις τιμές που περιέχουν. Για παράδειγμα, 

```{r}
a = c(12,21,13) ## δεν φέρει ονόματα
names(a) = c("Gene1", "Gene2", "Gene3")
a

a["Gene1"] == a[1]
```


## Matrix
Η δομή matrix ειναι μια *δισδιάστατη* δομή δεδομένων. Μοιάζει με τα vectors. 
Έτσι:

* Η δομή matrix έχει πάντα τον ίδιο τυπο δεδομένων (πχ chr, numeric)

# Κατασκευή

*Με default τιμή*
```{r}
## όλα τα κελιά της m θα έχουν την τιμή 0
m = matrix(0, nrow=5, ncol=3)
m
```

*Με vector τιμών*
```{r}
m1 = matrix(1:10, nrow=5, ncol=2) ## δεκα τιμές σε δέκα κελιά
m1
m2 = matrix(1:10, nrow=5, ncol=2, byrow=TRUE)
m2
m3 = matrix(1:10, nrow=5, ncol=2, byrow=FALSE)
m3
# Ποια ειναι η default συμπεριφορά της R για το πώς να βάλει τις τιμές στα κελιά;
```

*κυκλικά με vector τιμών*
```{r}
m1 = matrix(1:2, nrow=2, ncol=5, byrow=TRUE)
m1
m2 = matrix(1:2, nrow=2, ncol=5, byrow=FALSE)
m2
```

*Τι γίνεται αν το μέγεθος του matrix δεν ειναι πολλαπλάσιο του μήκους του vector*
```{r}
m = matrix(1:3, nrow=2, ncol=5)
m
```

*Τι γίνεται αν το μέγεθος του matrix είναι μικρότερο από το μήκος του vector*
```{r}
m = matrix(1:15, nrow=2, ncol=5)
m
```

## Ονόματα γραμμών και στηλών
Στα matrix μπορούμε να έχουμε ονόματα γραμμών και στηλών
```{r}
m = matrix(1:10, nrow=2, ncol=5, byrow=TRUE)
colnames(m) = letters[1:5]
rownames(m) = c("a", "aa")
m
```

## Πρόσβαση στα στοιχεία του matrix
Για να έχουμε πρόσβαση στα στοιχεία του matrix χρησιμοποιούμε τον operator `[,]`. Πριν το κόμμα ',' αναφερόμαστε στις γραμμές και μετά το κόμμα στις στήλες. Για παράδειγμα:
```{r}
a = matrix(1:10, nrow=2)
a[1,2] ## στοιχείο 1ης γραμμής 2ης στήλης

a[1, 1:2] ## στοιχεία 1ης γραμμής, στήλης 1 και 2

a[1:2, c(1,2,4)] # στοιχεία γραμμών 1,2 και στηλών 1,2,4. Θα ειναι matrix 2x3. 
```

Αν υπάρχουν ονόματα γραμμών στηλών, τότε η πρόσβαση μπορεί να γίνει και με τα ονόματα
```{r}
a = matrix(1:10, nrow=2)
colnames(a) <- letters[1:5]
rownames(a) <- c("x", "y")
a[,"c"] ## στοιχείο 1ης γραμμής 2ης στήλης
a["x",]
```

*ΠΑΡΑΤΗΡΗΣΗ: Για να πάρουμε όλη την γραμμή ή την στήλη χρησιμοποιούμε 'άδεια' θέση στις γραμμές ή στις στήλες*. Δείτε το παράδειγμα παραπάνω. 


## Data frames
To data.frame μοιάζουν με matrix όμως...

* στις στήλες τους μπορούν να περιέχουν δεδομένα διαφορετικών τύπων. Για παράδειγμα, η πρώτη στήλη μπορεί να περιέχει αριθμούς, ενώ η δευτερη στήλη χαρακτήρες

Ο πιο συνηθισμένος τρόπος να κατασκευάσουμε data.frame είναι με την συνάρτηση read.table, όπου διαβάζουμε δεδομένα που βρίσκονται σε αρχειο και έχουν μορφή πίνακα

```{r}
##a = matrix(rnorm(20000, 1000, 10), nrow=200)
##colnames(a) = paste("Individual",1:100, sep="")
##rownames(a) = paste("Gene", 1:200, sep="")
##write.table(a, file="geneExpressionExample.csv", sep="\t", col.names=TRUE, row.names = F, quote = F)

## Διαβάζουμε το αρχείο με την συνάρτηση read.table
## Έχει επικεφαλίδα, άρα βάζουμε header=TRUE
a = read.table(file="geneExpressionExample.csv", header=TRUE)
dim(a) # οι Διαστάσεις του a
head(a) # οι πρώτες έξι γραμμές
```



#--------------------------------------------------------------#
---
title: "3_2020_R_lists_examples.Rmd"
author: "pavlos pavlidis"
date: "2/11/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## Examples with Lists

### Position weight matrices (PWM)

Let's first assume a long DNA sequence

```{r}
alphabet = c('A', 'C', 'G','T')
a <- sample(alphabet, 1000, replace = TRUE)
paste(a, collapse="")
```

Let's assume that this sequence is a regulatory sequence of a gene. This means that on it *at certain locations* specific *proteins are able to bind*. Such proteins are called transcription factors and they are responsible for the *regulation of gene expression for the gene located next to the region*. For an animation please see the video here:

[https://www.youtube.com/watch?v=vVDDV8jH1Wc](https://www.youtube.com/watch?v=vVDDV8jH1Wc)

We would like to detect the regions on the DNA, where a specific transcription factor can bind. Let's assume the following:

1. TF can bind to a region of specific length \( l = 13 \). 
2. It does not bind to a single specific sequence, e.g. ACCACTTTAAGTA, but it binds to a wide range of sequences with some constraints. For example, at the first position always an 'A' may be expected, while in the second position, C is preferred, but also a G is acceptable and so on. 

Such preferences in the binding site can be described by a position weight matrix (PWM). 
```{r}
library(Biostrings)
data(HNF4alpha)
pfm <- consensusMatrix(HNF4alpha)[1:4,]
pfm
###pwm1<- PWM(pfm)
pwm.tmp <- pfm / colSums(pfm) 
pwm.tmp
pwm.tmp <- pwm.tmp/0.25
pwm.tmp
pwm2 <- log2(pwm.tmp)
pwm2
```

Now, we can think of how to match the **pwm** to the string
```{r}
pwm2
a
diag(pwm2[a[1:13],])
score = sum(diag(pwm2[a[1:13],]))
print(score)
v <- c()
for(i in 1:988){
  score = sum(diag(pwm2[a[i:(i+12)],]))
  v <- c(v, score)
}
plot(v)
 ```
