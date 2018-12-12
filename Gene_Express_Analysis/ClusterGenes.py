import csv
import numpy as np
import matplotlib
import statistics
import matplotlib.pyplot as plt
import sys

line_cnt = 0
Samples = [None] * 4764
geneExpress = {}
geneExpress2 = {}

sampleName = [None] * 19
geneName = set([])
geneName2 = set([])

med = [None] * 19
med_ratio = [None] * 19
norm_med = [None] * 19
norm_med_ratio = [None] * 19

normGenes = [None] * 4764
normGenesHash = {}


def distance(gene1, gene2):
    # type: (bytearray, bytearray) -> float
    gene2 = normalize(gene1, gene2)
    return np.mean(np.absolute(np.subtract(gene1, gene2)))


def remove_from_set(entries, the_dict):
    temp = {}
    for key in entries:
        if key in the_dict:
            temp[key] = the_dict[key]
    return temp


def normalize(gene1, gene2):
    # type: (bytearray, bytearray) -> bytearray
    ratio = abs(np.mean(gene1) / np.mean(gene2))
    return np.multiply(gene2, ratio)


def convert_string_float(row):
    converted = [None] * len(row)
    for i in range(len(row)):
        converted[i] = float(row[i])
    return converted


### ---------------------------------- Main -------------------------------------------------------------------

input_file1 = '../19_Samples/19_samples.CPM_gt_1_in_19.gct'
with open(input_file1) as csvDataFile:
    print('   Reading file: ' + input_file1)
    csvReader = csv.reader(csvDataFile, delimiter='\t', quotechar='"')

    for row in csvReader:
        line_cnt += 1

        if line_cnt == 3:
            for n in range(2, 21):
                sampleName[n - 2] = row[n]

        if line_cnt > 3:
            geneExpress[row[1]] = convert_string_float(row[2:21])
            geneName.add(row[1])

    line_cnt = 0

input_file2 = '../19_Samples/rna_tissue.tsv'
with open(input_file2) as csvDataFile2:
    print('   Reading file: ' + input_file2)
    csvReader2 = csv.reader(csvDataFile2, delimiter='\t', quotechar='"')

    for row in csvReader2:
        line_cnt += 1
        if row[2] == "skin":
            geneName2.add(row[1])
            geneExpress2[row[1]] = row[3]

#
# Compute intersection of two gene lists
#
interGenes = set(geneName).intersection(geneName2)
geneList = sorted(list(interGenes))

print('   Number of genes in intersection: ' + str(len(interGenes)))

numGenes = len(interGenes)
geneExpress = remove_from_set(interGenes, geneExpress)
geneExp2ress2 = remove_from_set(interGenes, geneExpress2)

express_diff = [None] * numGenes

#
# Compute the pairwise difference between gene expression levels
#
print('   Compute differences between gene expression levels ...')
for i in range(numGenes):
    express_diff[i] = [None] * numGenes

    for j in range(numGenes):
        if j > i:
            express_diff[i][j] = distance(geneExpress[geneList[i]], geneExpress[geneList[j]])

max_diff = 7 # 10 # 9 # 8 # 7
gene_cluster = [None] * numGenes

#
# Compute clusters of genes such that the difference of expression levels between genes in
# a cluster is below the theshold max_diff.
#
print('   Computing clusters for max. difference: ' + str(max_diff))
for i in range(numGenes):
    for j in range(numGenes):
        if j > i and express_diff[i][j] < max_diff:
            if gene_cluster[i] is None:
                if gene_cluster[j] is None:
                    cluster = set([i, j])                               ### create new cluster
                    gene_cluster[i] = cluster
                    gene_cluster[j] = cluster
                else:
                    gene_cluster[j].add(i)                              ### add i to cluster of j
                    gene_cluster[i] = gene_cluster[j]
            else:
                if gene_cluster[j] is not None:
                    union = set(gene_cluster[i]).union(gene_cluster[j]) ### merge clusters 
                    for gene_idx in union:                              ### point each element of the cluster to the updated cluster
                        gene_cluster[ gene_idx ] = union
                else:
                    gene_cluster[i].add(j)                              ### add j to cluster of i
                    gene_cluster[j] = gene_cluster[i]

#
# Compute the list of unique clusters
#
cluster_list = []
for cluster in gene_cluster:
    if (    (not cluster in cluster_list)
        and (cluster != None)):
        cluster_list.append( cluster )

#
# Print the list of clusters:
#
cluster_no = 0
for cluster in cluster_list:
    cluster_no += 1
    print('   cluster ' + str(cluster_no) + ':')
    # print(cluster)
    if (cluster != None):
        sys.stdout.write('{')
        cnt = 0
        for gene_idx in cluster:
            if (cnt > 0):
                sys.stdout.write(', ')
            # sys.stdout.write(str(gene_idx))
            sys.stdout.write(geneList[gene_idx])
            cnt += 1
        print('}')
        

print('\nNunmber of clusters: ' + str(len(cluster_list)))
print('End of execution')
