
"A. Finds genotype and methylation within a gene + creates a table with gene name and corresponding number of methylation sites and genotypes"

Run code: temp0.1.R

Input data:
1.	control.geno.file.txt
2.	case.geno.file.txt

3.	case.log.count.50perc.miss.txt
4.	control.log.count.nomiss.txt

5.	meth.case.txt
6.	meth.control.txt

7.	covariates_file.txt

Output files:

1.	genes.txt 
Gives Gene names. Only TSC1 gene is included in this example dataset

2.	mat.id.txt
Gives Subject names
Subject_ID_64
Subject_ID_65
Subject_ID_66
Subject_ID_67
Subject_ID_68
Subject_ID_69
…

3.	allgene.matrix.txt

Gives a table with columns GeneID, Gene.name, num.Genotype, num.Methylation.
First row: 1, TMEM165, 16, 10
Second row: 2, ZDHHC16, 5, 18
Third row: 3, PHC1, 0, 7
Fourth row: 4, RHD, 6, 3
Fifth row: 5, NIT1, 8, 2
and so on. 

4.	Generates data files such as, geneid_1_TMEM165.txt, geneid_2_ ZDHHC16.txt, geneid_3_ PHC1.txt, etc.

[generic name: geneid_GeneID_Gene.name.txt]



"B. Shell script for using TiMEG Pipeline and generating result files having p-values"

Run code: submit3.sh
Input data:
1.	geneid_GeneID_Gene.name.txt files

Output files:
Eg. fileid.1.TMEM165.1.2.txt

[generic name:
fileid.GeneID.Gene.name.methno.genotypeno.txt]


Run code: submit_R4.sh

Input code: 

1. submit_preR3.sh
2. TiMEG.R

Output files:

Eg. Rresult_TSC1.1.txt
[generic: Rresult_ Gene.name.j.txt, where, j=1,2,…, num.Genotype]



