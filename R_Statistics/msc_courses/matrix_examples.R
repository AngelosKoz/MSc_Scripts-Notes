---
title: "Examples with matrices and vectors"
author: "pavlos pavlidis"
date: "2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Examples with Matrices

### Matrices with one initial value
```{r}
m <- matrix(0, nrow=10, ncol=5)
```

### Matrices with one character value
```{r}
m <- matrix('a', nrow=10, ncol=5)
```

### Matrices initialized with a vector of elements with as many elements as the total number of matrix cells
```{r}
m <- matrix(1:50, nrow=10, ncol=5)
```

### Matrices initialized with a smaller vector than the number of total elements
```{r, echo=TRUE, eval=TRUE}
m1 <- matrix(1:10, nrow=10, ncol=5)
m2 <- matrix(1:10, nrow=10, ncol=5, byrow=TRUE)
```

### Matrices initialized with a vector that cannot be divided exactly with the total number of elements
```{r, echo=TRUE, eval=TRUE}
m1 <- matrix(1:9, nrow=10, ncol=5)
```

## Vectors

### Subsetting with negative indexes
```{r, echo=TRUE, eval=TRUE}
v <- c(1,2,3)
v[-c(1,3)]
```

### Subsetting with logical vectors
Here the logical vector should be the same length, or to be divided exactly with the vector that we want to subset
```{r, echo=TRUE, eval=TRUE}
v <- seq (1, 100, 2)
lv <- c(T,F)
v[lv]
```



### Named vectors

```{r, echo=TRUE, eval=TRUE}
v <- c(1,2,3)
names(v) <- c('a', 'b', 'c')

v[1]
v['b']
```
This means that vectors can be also used as hashes


#### Example with named vectors
##### Constructing the complementary DNA strand

```{r, echo=TRUE, eval=TRUE}
alphabet <- c('A', 'C', 'G', 'T')
names(alphabet) <- c('T', 'G', 'C', 'A')

dna <- sample(alphabet, 1000, TRUE)
compDNA <- alphabet[dna]

paste(dna, collapse='')
paste(compDNA, collapse='')
```

## A very useful link for subsetting
[http://adv-r.had.co.nz/Subsetting.html](http://adv-r.had.co.nz/Subsetting.html)
