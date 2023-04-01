#  `PyAGH` manual

## Introduction

`PyAGH` is a Python package developed to simply calculate the relationship matrix and handle the data required for the calculation (supports high marker density typing data; large pedigree data; microbiome data; additive, dominant and epistatic effects relationship matrix; relationship matrix across populations; pedigree tracing and visualization).

## Installation

It is recommended to use `Python 3.9`

### Using pypi

`PyAGH`can be installed by using `pip install`.

```shell
pip install PyAGH
```

## Contents

[loadEgPed()](#1)

[loadEgGeno()](#2)

[sortPed()](#3)

[selectPed()](#4)

[makeA()](#5)

[makeD()](#6)

[makeG()](#7)

[makeG_inter()](#8)

[makeM()](#9)

[makeH()](#10)

[coefKinship()](#11)

[coefInbreeding()](#12)

[cluster()](#13)

[pca()](#14)

[heat()](#15)

[gragh()](#16)

## Functions

---

`<span id="1">loadEgPed()`						Function to load pedigree example data

`<span id="2">loadEgGeno()`                      Function to load genotype example data

---

#### Description

Function `loadEgPed()` and `loadEgGeno()` can loads the pedigree and genotype example data sets built into `PyAGH`

#### Usage

```python
import PyAGH
ped = PyAGH.loadEgPed()
genofile = PyAGH.loadEgGeno() 
```

#### Arguments

None

#### Value

Return a pandas.DataFrame of example pedigree and a genotype file path.

---

`<span id="3">sortPed()`                         Function to sort the pedigree data

---

#### Description

Sort the pedigree  data according to the correct birth date of individuals and check for various errors in the pedigree like offspring born before its parents, same offspring have different parents, loop in pedigree and etc.

#### Usage

```python
ped_sorted = PyAGH.sortPed(ped)
```

#### Arguments

ped       a pandas.DataFrame with three columns  ID, SIRE, DAM. Null values in the data are replaced by 0.

#### Value

Return a pandas.DataFrame of sorted pedigree.

---

`<span id="4">selectPed()`                   Function to select pedigree

---

#### Description

`selectPed()` function can select pedigree based on specific individuals and generations.

#### Usage

```python
ped_selected = PyAGH.selectPed(data=ped,id=['9','10'],generation=3)
```

#### Arguments

data                          pandas.DataFrame of sorted pedigree data

id							   a list of individuals ID which want to be selected

generation			   a int type number represents the number of generations to be traced

#### Value

Return a pandas.DataFrame of selected pedigree.

---

`<span id="5">makeA()`				          Function to calculate the relationship matrix based on pedigree information

---

#### Description

Construct relationship matrix based on pedigree information for additive effect. Option to use sparse matrix for memory saving but slower.

#### Usage

```python
A = PyAGH.makeA(ped_sorted,Sparse=False)
```

#### Arguments

ped_sorted             pandas.DataFrame of sorted pedigree data

Sparse                      logic, if True uses sparse matrix

#### Value

Return a tuple with 2 elements. The first one  a numpy array or a sparse matrix. Second is a list of ID.

---

`<span id="6">makeD()`					   Function to calculate the relationship matrix based on pedigree information

---

#### Description

Construct relationship matrix based on pedigree information for dominant effect. Option to use  multithreading when there are multiple CPU.

#### Usage

```python
D = PyAGH.makeD(ped_sorted,multi=1)
```

#### Arguments

ped_sorted             pandas.DataFrame of sorted pedigree data

multi                        a int type number for multithreading.  Default value is 1. This function uses multi-threaded calculation by 								 default, if your computer has more than one cpu, you can set the value of multi equal to the number of cpu.

#### Value

Return a tuple with 2 elements. The first one  a numpy array . Second is a list of ID.

---

`<span id="7">makeG()`				      Function to calculate a relationship matrix from marker data, G matrix.

---

#### Description

Construct relationship matrix based on genotype data for additive effect. Option to use different methods. There are 0,1,2,3 four methods for G matrix construction. 0 and 1 methods are used in single population  developed by Yang et al. (2010) and VanRaden (2008). 2 and 3 methods are used in cross population developed by Wientjes et al. (2017) and Chen et al. (2013).

#### Usage

```python
G = PyAGH.makeG(File=genofile,method=1,File_list=False,n1=0,n2=0)
```

#### Arguments

File                           genotype file path. Input file type is traw of plink.

method                   int type number represents different methods.

File_list					logic, if True, File arguments need a file that writes the paths of genotype file for multiple  chromosomes. Then, 								 the function will calculate in chromosomes to use less memory.

n1, n2					  when the method is 2 or 3, the number of individuals in each population needs to be provided.

#### Value

Return a tuple with 2 elements. The first one  a numpy array . Second is a list of ID.

---

`<span id="8">makeG_inter()`		 Function to calculate a relationship matrix from marker data for epistatic effect

---

#### Description

Construct relationship based on genotype data for epistatic effect.  There are four epistatic effect algorithms represent as 'aa', 'dd', 'ad', 'da'. For the detailed formula, please refer to Xu (2013).

#### Usage

```python
G_inter = PyAGH.makeG_inter(geno,method,multi=1)
```

#### Arguments

geno           		 a numpy matrix of genotype which code as 0,1 and 2. Rows are individuals and columns are SNPs.

method               str value of 'dd','aa','ad' or 'da'.

multi 				   a int type number for multithreading.  Default value is 1. This function uses multi-threaded calculation by 							  default, if your computer has more than one cpu, you can set the value of multi equal to the number of cpu.

#### Value

Return a tuple with 2 elements. The first one  a numpy array . Second is a list of ID.

---

`<span id="9">makeM()`                   Function to calculate a relationship matrix based on microbiome data

---

#### Description

Construct relationship based on microbiome data.

#### Usage

```python
M = PyAGH.makeM(otu)
```

#### Arguments

otu                  OTU data in numpy ndarray type. Rows are individuals and columns are OTUs.

#### Value

Return a tuple with 2 elements. The first one  a numpy array . Second is a list of ID.

---

`<span id="10">makeH()`                Function to construct relationship based on pedigree and genotype data

---

#### Description

The `makeH()` function combines information of pedigree and genotypes to construct relationships for both genotyped and ungenotyped individuals used in single-step genomic best linear unbiased prediction (SSGBLUP) method.

#### Usage

```python
H = makeH(G,A,w=0.05)
```

#### Arguments

  G                a list with two elements, G[0] is relationship matrix based on genotype; G[1] is id series of genotyped individuals.

  A                 a list with two elements, A[0] is relationship matrix based on pedigree; A[1] is id series of all individuals.

  w                Default value is 0.05. The weights of the G-matrix and the A-matrix

#### Value

Return a tuple with 2 elements. The first one  a numpy array . Second is a list of ID.

---

`<span id="11">coefKinship()`                Calculates the relationship coefficients using relationship matrix

---

#### Description

Calculates the relationship coefficients using relationship matrix

#### Usage

```python
coef_kinship = PyAGH.coefKinship(A)
```

#### Arguments

A                  a list with two elements, A[0] is relationship matrix; A[1] is id series

#### Value

Return a pandas.DataFrame with relationship coefficients.

---

`<span id="12"> coefInbreeding()`           Calculates inbreeding coefficients for individuals

---

#### Description

Calculates the inbreeding coefficients using relationship matrix.

#### Usage

```python
coef_inbreeding = PyAGH.coefInbreeding(A)
```

#### Arguments

A                  a list with two elements, A[0] is relationship matrix; A[1] is id series

#### Value

Return a pandas.DataFrame with inbreeding coefficients.

---

`<span id="13">cluster()`                          Plot the cluster result of relationship matrix

---

#### Description

Cluster analysis of the relationship matrix and plot the result.

#### Usage

```python
cluster_example = PyAGH.cluster(A,color_threshold=0.9,above_threshold_color='gray')
plt.savefig('cluster_example.png', facecolor='w',dpi=300)
```

#### Arguments

A                  a list with two elements, A[0] is relationship matrix; A[1] is id series

color_threshold                 an optional parameter that can be used to set at which clustering distance the lines in the dendrogram plot should be changed to above_threshold_color.

#### Value

Figure

---

`<span id="14">pca()`                                 Principal component analysis of relationship matrix

---

#### Description

Principal component analysis of relationship matrix and plot the result.

#### Usage

```python
pca_example = PyAGH.pca(A,group,color)
plt.savefig('pca_example.png', facecolor='w',dpi=300)
```

#### Arguments

A                  a list with two elements, A[0] is relationship matrix; A[1] is id series

group          a list containing information about the group of individuals.

color            a list containing color of different group.

#### Value

Figure

---

`<span id="15">heat()`                          Plot the heatmap of relationship matrix

---

#### Description

 Plot the heatmap of a relationship matrix

#### Usage

```python
heatmap_example = PyAGH.heat(A)
plt.savefig('heatmap_example.png', facecolor='w',dpi=300)
```

#### Arguments

A                  a list with two elements, A[0] is relationship matrix; A[1] is id series

#### Value

Figure

---

`<span id="16">gragh()`                       Plot family tree in three generations of a individual

---

#### Description

 Plot the heatmap of a relationship matrix

#### Usage

```python
import graphviz
ped_selected = PyAGH.selectPed(data=ped,id=['17'],generation=3)
p = PyAGH.gragh(ped_selected,color_sire="blue",color_dam="red",fillcolor_sire="lightblue",fillcolor_dam="pink")
graphviz.Source(p)
```

#### Arguments

 data_ord                  the pedigree of one individual in three generations.

color_sire                  the color of sire.

color_dam                 the color of dam.

fillcolor_sire               the fillcolor of sire.

fillcolor_dam              the fillcolor of dam.

#### Value

Figure
