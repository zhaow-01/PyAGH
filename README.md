# PyAGH

![](https://raw.githubusercontent.com/zhaow-01/PyAGH/main/picture/logo.png)

## Introduction

`PyAGH` is a Python package developed for calculating kinship matrix using pedigree or genotype data as well as  processing, analysis and visualization for data. `PyAGH` provides fast and concise methods for calculating Amatrix (based on pedigree), Gmatrix (based on genotype) and Hmatrix (based on pedigree and genotype) used in breeding. `PyAGH`supports a variety of mainstream algorithms, you can easily choose different algorithms for calculation. With the obtained relationship matrix, you can use `PyAGH` for fast visualization, including PCA analysis, Cluster analysis and drawing Heatmaps.

To make it easier to use pedigree information, `PyAGH ` also provides targeted tools for specific needs, such as detecting pedigree errors,  selection individuals, pedigree visualization, sorting pedigree, calculating inbreeding coefficients and calculating  relationship coefficients . The latest stable release of the software can be installed conveniently via pip .

### Target Audience

The target audience of `PyAGH` includes:

- Students and researchers in the field of breeding and genetics, particularly those who want to perform genomic prediction.
- Government, enterprises, or other entities who expect efficient processing of pedigree and genomic information.

### Technical Features

- Provide a variety of  methods for calculating kinship matrix including combined reference population.
- Fast and support tens of millions of genotypes.

### Main Functions

Currently, `PyAGH` mainly provides the following methods:

- **Pedigree sorting:  **Provides methods to quickly obtain the correct pedigree in chronological order of birth clean multiple types of data error.
- **Pedigree selection:** Provides methods to select pedigree based on specific individuals and generations.
- **Calculation of coefficients:** Provides methods to easily obtain the inbreeding and relationship coefficients of specific individuals.
- **Kinship matrix: ** Provides different methods to calculating Amatix, Gmatrix and Hmatrix.
- **Visualization:** Provides methods to display results.

****

## Installation

It is recommended to use `Python 3.7, 3.8, 3.9`

### Using pypi

`PyAGH`can be installed by using `pip install`.Before installing `PyAGH`, make sure that you have support the   C++ compilation.

```shell
pip install PyAGH
```

## Examples

#### Load the example data

`loadEgPed()` and `loadEgGeno()` functions can load the example data.

```python
import PyAGH
ped = PyAGH.loadEgPed()
genofile = PyAGH.loadEgGeno() #is a filename
ped
```

| index |   id | sire |  dam |
| :---: | ---: | ---: | ---: |
|   0   |    9 |    1 |    2 |
|   1   |   10 |    3 |    4 |
|   2   |   11 |    5 |    6 |
|   3   |   12 |    7 |    8 |
|   4   |   13 |    9 |   10 |
|   5   |   14 |   11 |   12 |
|   6   |   15 |   11 |    4 |
|   7   |   16 |   13 |   15 |
|   8   |   17 |   13 |   14 |

#### Sort pedigree

`sortPed()` is a function that can sort the pedigree according to the correct birth date of individuals and check for various errors in the  pedigree.

```python
ped_sorted = PyAGH.sortPed(ped)
ped_sorted
```

| index |  id  | sire | dam  |
| :---: | :--: | :--: | :--: |
|   0   |  1   |  0   |  0   |
|   1   |  3   |  0   |  0   |
|   2   |  5   |  0   |  0   |
|   3   |  7   |  0   |  0   |
|   4   |  2   |  0   |  0   |
|   5   |  4   |  0   |  0   |
|   6   |  6   |  0   |  0   |
|   7   |  8   |  0   |  0   |
|   8   |  9   |  1   |  2   |
|   9   |  10  |  3   |  4   |
|  10   |  11  |  5   |  6   |
|  11   |  12  |  7   |  8   |
|  12   |  15  |  11  |  4   |
|  13   |  13  |  9   |  10  |
|  14   |  14  |  11  |  12  |
|  15   |  16  |  13  |  15  |
|  16   |  17  |  13  |  14  |

#### Select pedigree

`selectPed()` function can select pedigree based on specific individuals and generations.

```python
ped_selected = PyAGH.selectPed(data=ped,id=['9','10'],generation=3)
ped_selected
```

| index |  id  | sire | dam  |
| :---: | :--: | :--: | :--: |
|   0   |  1   |  0   |  0   |
|   1   |  3   |  0   |  0   |
|   2   |  2   |  0   |  0   |
|   3   |  4   |  0   |  0   |
|   4   |  9   |  1   |  2   |
|   5   |  10  |  3   |  4   |

#### Calculate kinship matrix

`makeA()`, `makeG()` and `makeH()` functions can calculate the kinship matrix using various methods based on pedigree , genotype and both together, respectively.

```python
A = PyAGH.makeA(ped_sorted)
G = PyAGH.makeG(File=genofile,method=1,File_list=False)
H = PyAGH.makeH(G,A,w=0.05)
####A,G and H are lists with 2 elements. The first one is kinship matrix in numpy.ndarray type and the other one is id labels in pandas.Series type. 
```

#### Calculate coefficients

```python
coef_kinship = PyAGH.coefKinship(A)
coef_kinship
```

| index |  ID  |    F    |
| :---: | :--: | :-----: |
|   0   |  1   | 0.0000  |
|   1   |  3   | 0.0000  |
|   2   |  5   | 0.0000  |
|   3   |  7   | 0.0000  |
|   5   |  4   | 0.0000  |
|  11   |  12  | -0.0000 |
|  ...  | ...  |   ...   |
|  12   |  15  | -0.0000 |
|  13   |  13  | -0.0000 |
|  14   |  14  | -0.0000 |
|  15   |  16  | 0.0625  |
|  16   |  17  | -0.0000 |

```python
coef_inbreeding = PyAGH.coefInbreeding(A)
coef_inbreeding
```

|      |  ID1 |  ID2 |        r |
| ---: | ---: | ---: | -------: |
|    0 |    1 |    1 | 1.000000 |
|    1 |    1 |    3 | 0.000000 |
|    2 |    1 |    5 | 0.000000 |
|    3 |    1 |    7 | 0.000000 |
|    4 |    1 |    2 | 0.000000 |
|  ... |  ... |  ... |      ... |
|  148 |   14 |   16 | 0.121268 |
|  149 |   14 |   17 | 0.500000 |
|  150 |   16 |   16 | 1.000000 |
|  151 |   16 |   17 | 0.333486 |
|  152 |   17 |   17 | 1.000000 |

153 rows × 3 columns

#### Visualization

```python
cluster_example = PyAGH.cluster(A)
```

![aaa](C:\Users\zhaow\Documents\aaa.png)

```python
group=['1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','2','2']
pca_example = PyAGH.pca(A,group=group)
```

![pca](C:\Users\zhaow\Documents\pca.png)

```python
heatmap_example = PyAGH.heat(A)
```

<img src="C:\Users\zhaow\Documents\heat.png" alt="heat" style="zoom:150%;" />

```python
import graphviz
ped_selected = PyAGH.selectPed(data=ped,id=['17'],generation=3)
p = PyAGH.gragh(ped_selected)
graphviz.Source(p)
```



![output](C:\Users\zhaow\Documents\output.svg)

## License

`PyAGH` is MIT licensed.