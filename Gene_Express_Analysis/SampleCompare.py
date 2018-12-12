from scipy import stats
import csv
import numpy
import matplotlib
import statistics
import matplotlib.pyplot as plt
from operator import itemgetter
import copy


line_cnt = 0
Samples = [None] * 4764
geneExpress = {}
geneExpress2 = {}

sampleName = [None] * 19
geneName = set([])
geneName2 = set([])

z = [None] * 19
z2 = [None] * 30
avg = [None] * 19
avg2 = [None] * 30
norm_avg = [None] * 19
norm_avg2 = [None] * 30
norm_z = [None] * 19
norm_z2 = [None] * 30

normGenes = [None] * 4764
normGenesHash = {}


def convert_string_float(row):
    converted = [None] * len(row)
    for i in range(len(row)):
        converted[i] = float(row[i])
    return converted


def remove_from_set(entries, the_dict):
    temp = {}
    for key in entries:
        if key in the_dict:
            temp[key] = the_dict[key]
    return temp


def normalize(norm, the_dict):
    for key in the_dict:
        numpy.divide(the_dict[key], norm)

###
### Reads Dermal Papilla data and stores into dictionaries
###
def read_gene_express_data():
    line_cnt = 0
    global geneExpress
    global geneName
    global sampleName
    with open('../19_Samples/CPM_gt_1.DESeq_norm.FPKM.txt') as csvDataFile:
        csvReader = csv.reader(csvDataFile, delimiter='\t', quotechar='"')

        for row in csvReader:
            line_cnt += 1

            if (line_cnt == 3):
                for n in range(2, 21):
                    sampleName[n - 2] = row[n]

            if (line_cnt > 3):
                geneExpress[row[1]] = convert_string_float(row[2:21])
                geneName.add(row[1])

###
### Reads Fibroblast data and stores into dictionaries
###
def read_gene_express2_data():
    line_cnt = 0
    global geneExpress2
    global geneName2
    with open('../19_Samples/GSE63577_counts_rpkm_exvivo_jenage_data.csv') as csvDataFile2:
        csvReader2 = csv.reader(csvDataFile2, delimiter=',', quotechar='"')

        for row in csvReader2:
            line_cnt += 1
            if (line_cnt > 1):
                geneName2.add(row[1])
                geneExpress2[row[1]] = convert_string_float(row[4:34])


###
### Normalize expression set
###
def normalize_express():
    global tempSample
    global normExpress
    global normTempSample
    global norm_avg
    global norm_z
    for i in range(19):           ### iterate through 19 samples
        del tempSample[:]
        for key in geneExpress:   ### iterate through ~2000 genes
            tempSample.append(geneExpress[key][i])        ### 19 DP samples

        tempSample = numpy.log10(tempSample)
        tempSample = tempSample.tolist()

        avg[i] = statistics.mean(tempSample)
        z[i] = numpy.std(tempSample)

        for key in geneExpress:
            normExpress[key][i] = numpy.math.log10(normExpress[key][i])
            normExpress[key][i] = (normExpress[key][i] - avg[i]) / z[i]

        normTempSample = numpy.divide(numpy.subtract(tempSample, avg[i]), z[i])
        norm_avg[i] = statistics.mean(normTempSample)
        norm_z[i] = numpy.std(normTempSample)

###
### Normalize expression2
###
def normalize_express2():
    global tempSample2
    global normExpress2
    global normTempSample2
    global norm_avg2
    global norm_z2
    for i in range(30):           ### iterate through 30 Fibroblast samples
        del tempSample2[:]
        for key in geneExpress:   ### iterate through ~2000 genes
            if geneExpress2[key][i] == 0:
                tempSample2.append(10.0**(-10))
            else:
                tempSample2.append(geneExpress2[key][i])

        tempSample2 = numpy.log10(tempSample2)
        tempSample2 = tempSample2.tolist()

        avg2[i] = statistics.mean(tempSample2)
        z2[i] = numpy.std(tempSample2)

        for key in geneExpress:
            if normExpress2[key][i] == 0:
                normExpress2[key][i] = 10.0**(-10)
            normExpress2[key][i] = numpy.log10(normExpress2[key][i])
            normExpress2[key][i] = (normExpress2[key][i] - avg2[i]) / z2[i]

        normTempSample2 = numpy.divide(numpy.subtract(tempSample2, avg2[i]), z2[i])
        norm_avg2[i] = statistics.mean(normTempSample2)
        norm_z2[i] = numpy.std(normTempSample2)



### ---------------------------------- Main -------------------------------------------------------------------
read_gene_express_data()
read_gene_express2_data()

### Compute intersection of gene sets:
interGenes = set(geneName).intersection(geneName2)

numGenes = len(interGenes)
geneExpress = remove_from_set(interGenes, geneExpress)
geneExpress2 = remove_from_set(interGenes, geneExpress2)


print("\nSimilar genes before normalization")

preDiff = [None] * numGenes
i = 0

for key in geneExpress:
    set = (geneExpress[key][3], geneExpress[key][4], geneExpress[key][5],
           geneExpress[key][6], geneExpress[key][12], geneExpress[key][17])
    s1 = numpy.mean(geneExpress[key])
    s2 = numpy.mean(geneExpress2[key])
    preDiff[i] = [None] * 2
    ### Calculate t-stat
    preDiff[i][0] = abs((s1-s2)/numpy.sqrt(numpy.square(numpy.std(set))/6.0 +
                                        numpy.square(numpy.std(geneExpress2[key])/6.0)))
    ### Calculate p-value
    preDiff[i][1] = stats.t.sf(numpy.abs(preDiff[i][0]), 12-1)*2.0
    i += 1

num = 0
for i in range(numGenes):
    if preDiff[i][1] > 0.05:
        num += 1
print num


### Normalize:

tempSample2 = []
tempSample = []
normExpress = copy.deepcopy(geneExpress)
normExpress2 = copy.deepcopy(geneExpress2)

normalize_express()
normalize_express2()

### Print average and SD before and after normalization
print('\naverage pre-norm: ')
print('   19 samples:  ' + str(avg))
print('   skin sample: ' + str(avg2))
print('\nstandard deviation pre-norm: ')
print('   19 samples:  ' + str(z))
print('   skin sample: ' + str(z2))
print('\naverage post-norm: ')
print('   19 samples:  ' + str(norm_avg))
print('   skin sample: ' + str(norm_avg2))
print('\nstandard deviation pre-norm: ')
print('   19 samples:  ' + str(norm_z))
print('   skin sample: ' + str(norm_z2))

#xaxis = range(0, 19, 1)

#plt.figure(1)

#plt.plot(xaxis, med, label="pre-norm medians")
#plt.plot(xaxis, norm_med, label="normalized medians")
#matplotlib.pyplot.xlabel("19 Samples", fontdict=None, labelpad=None)
#matplotlib.pyplot.ylabel("Normalized Expression Values", fontdict=None, labelpad=None)
#matplotlib.pyplot.legend(loc=1, prop={'size': 6})


i = 0
diff = [None] * numGenes

### Calculate differential analysis components
for key in normExpress:
    diff[i] = [None] * 9
    diff[i][0] = key
    set = (normExpress[key][3], normExpress[key][4], normExpress[key][5],
           normExpress[key][6], normExpress[key][12], normExpress[key][17])
    s1 = numpy.mean(set)
    s2 = numpy.mean(normExpress2[key])

    ### Calculate Fold Change
    diff[i][1] = (10 ** s1 + .31) / (10 ** s2 + .31)

    ### Determine up or down regulation
    diff[i][2] = 'DP up'
    if s1 - s2 < 0:
        diff[i][2] = 'DP down'

    ### Calculate t-stat
    if s1 < 0 and s2 < 0:
        diff[i][3] = 0
    else:
        diff[i][3] = abs((s1-s2)/numpy.sqrt(numpy.square(numpy.std(set))/6.0 +
                                            numpy.square(numpy.std(normExpress2[key])/30.0)))
    ### Calculate p-value
    diff[i][4] = stats.t.sf(numpy.abs(diff[i][3]), 12-1)*2.0

    ### Calculate log z-score difference
    diff[i][5] = s1 - s2

    ### Calculate standard deviation
    diff[i][6] = numpy.sqrt(numpy.square(numpy.std(set))/3.0 + numpy.square(numpy.std(normExpress2[key]))/30.0)

    ### Mean of FPKM dermal papilla
    diff[i][7] = numpy.mean((geneExpress[key][3], geneExpress[key][4], geneExpress[key][5],
                             geneExpress[key][6], geneExpress[key][12], geneExpress[key][17]))
    ### Mean of RPKM fibroblast
    diff[i][8] = numpy.mean(geneExpress2[key])
    i = i + 1

sortedDiff = sorted(diff, key=itemgetter(3), reverse=True)
backSortedDiff = sorted(diff, key=itemgetter(1))

plt.figure(1)

num = 0
print("\nSimilar genes after normalization")
for i in range(numGenes):
    if backSortedDiff[i][4] > 0.05:
        num += 1
print num

with open('../Python_Programs/DPvsFibroblast2.csv', 'w') as csvOutFile:
    fieldnames = ['Gene', 'Fold Change', 'p-value', 't-stat', 'z-score Difference',
                  'Standard Deviation', 'DP FPKM', 'FB RPKM']
    writer = csv.DictWriter(csvOutFile, dialect='excel', fieldnames=fieldnames)
    writer.writeheader()


    for i in range(numGenes):
        writer.writerow({'Gene': sortedDiff[i][0], 'Fold Change': round(float(sortedDiff[i][1]), 4),
                         'p-value': round(float(sortedDiff[i][4]), 4),
                         't-stat' : round(float(sortedDiff[i][3]), 4),
                         'z-score Difference': round(float(sortedDiff[i][5]), 4),
                         'Standard Deviation': round(float(sortedDiff[i][6]), 4),
                         'DP FPKM': round(float(sortedDiff[i][7]), 4),
                         'FB RPKM': round(float(sortedDiff[i][8]))})

print("\nGenes of highest difference: ")
for i in range(20):
    print('   ' + '{:12s}'.format(sortedDiff[i][0]) + '	' + '{:7.1f}'.format(sortedDiff[i][1]) + '    ' + sortedDiff[i][2])
    plt.plot([i], [numpy.mean(normExpress[sortedDiff[i][0]])], marker='o', markersize=6, color="red",  label="19 Samples Expression")
    plt.plot([i], [numpy.mean(normExpress2[sortedDiff[i][0]])],            marker='o', markersize=6, color="blue", label="Skin Samples Expression")

matplotlib.pyplot.xlabel("Top 20 Genes of Highest Difference", fontdict=None, labelpad=None)
matplotlib.pyplot.ylabel("Normalized Expression Values",       fontdict=None, labelpad=None)

# temp1 = plt.subplot()
# temp2 = plt.subplot()
# legend_elements = [temp1.plot([0], [0], marker='o', markersize=3, color="red", label="19 Samples Expression"),
#                   temp2.plot([0], [0], marker='o', markersize=3, color="blue", label="Skin Samples Expression")]
# plt.legend(handles=legend_elements)
plt.show()
