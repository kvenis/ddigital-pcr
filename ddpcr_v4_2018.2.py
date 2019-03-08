# coding: utf-8

# In[76]:

import os
import sys
import pandas as pd
import numpy as np
import math
import csv
from collections import OrderedDict
from datetime import date
from datetime import datetime
from pandas import rolling_median
from io import StringIO


infile = sys.argv[1]
outfolder = "/BCR_ABL ddPCR pipeline output"
print(outfolder.replace(" ", "\ "))
print(infile)
data = pd.read_excel(sys.argv[1], sheet_name='Sheet')


print("Testing....")

#data.head()
#data.columns
data1 = data.loc[:, ['Well', 'Sample','Target', 'Conc(copies/µL)', 'Accepted Droplets', 'DyeName(s)']]
data1
data1.rename(columns={'Conc(copies/µL)': 'Conc', 'Accepted Droplets' : 'Droplets'}, inplace=True)
data1

wt = data1.loc[data1['DyeName(s)'] == "HEX"]
mut = data1.loc[data1['DyeName(s)'] == "FAM"]
wt = wt.reset_index(drop=True)
#wt
mut = mut.reset_index(drop=True)
#mut
wildtype = wt.drop('Target', 1)
wildtype
wildtype.rename(columns={'Conc':'ABL1'}, inplace=True)
mutant = mut.drop('Target', 1)
mutant
mutant.rename(columns={'Conc':'BCR-ABL1' }, inplace=True)
mutant
new = pd.concat([wildtype,mutant['BCR-ABL1']], axis =1)
#pritn(new)
new['Ratio'] = new['BCR-ABL1'] / new['ABL1']
#pritn(new)

new['Ratio'] = new['Ratio'].fillna(0)
#new['Ratio' == 'inf'] = 0
new['Ratio'] = new['Ratio'].replace(np.inf, 0)
#new['Ratio'] = new['Ratio'].('inf', 0)
new['%Ratio'] = new['Ratio']*100
new['%IS'] = new['%Ratio']*0.89
new

full = new.copy()
full
sample = full.groupby(['Sample'])
sample

##finding outliers
def s(one):
    one['%Ratio_Mean'] = one['%Ratio'].mean()
    one['2*std'] = 2*one['%Ratio'].std()
    one['notaOutlier'] = (one['%Ratio'] == 0) | ((one['%Ratio'] > (one['%Ratio'].mean() - 2*one['%Ratio'].std())) & (one['%Ratio'] < (one['%Ratio'].mean() + 2*one['%Ratio'].std())))
    return one

real = sample.apply(s)
real


def testing_outlier(row):
    if row['notaOutlier'] == True:
        return '.'
    else:
        return 'Yes'
real = real.assign(test_outlier = real.apply(testing_outlier, axis=1))
real

# selecting wells with outliers / seperating from table to get the outliers table
real_1_out = real.loc[real['notaOutlier'] == False]
real_1_out
real_1_out.rename(columns={ 'Well': 'WELL', 'Sample': 'SAMPLE','Ratio':'RATIO','%Ratio':'%RATIO','Droplets':'DROPLETS','test_outlier':'OUTLIERWELL'}, inplace=True)
real_1_out

# Outlier table
outlier_table = real_1_out[['WELL','SAMPLE','BCR-ABL1','ABL1','%RATIO','%IS','DROPLETS','OUTLIERWELL']]
###outlier_table_1 = outlier_table.round({'BCR_ABL1':6,'ABL1':6,'RATIO':6,'%Ratio':6})
#outlier_table_1['BCR_ABL1'] = outlier_table_1['BCR_ABL1'].round().astype('int')
#outlier_table_1['ABL1'] = outlier_table_1['ABL1'].round().astype('int')
###outlier_table_1

# selecting wells with no outliers, gropuing and getting the sum
sum_wells = real.loc[real['notaOutlier'] == True]
sum_wells
sum_wells_1 = sum_wells.groupby(['Sample'], sort=False, as_index=False).agg({'BCR-ABL1': np.sum, 'ABL1': np.sum, 'Droplets': np.sum })
sum_wells_1
sum_wells_1.rename(columns={'Sample': 'SAMPLE','Droplets':'DROPLETS'}, inplace=True)
sum_wells_1

#getting ratios and % ratios for the sum values.
sum_wells_1['RATIO'] = sum_wells_1['BCR-ABL1'] / sum_wells_1['ABL1']

sum_wells_1['%RATIO'] = sum_wells_1['RATIO'] *100
sum_wells_1['%IS'] = sum_wells_1['%RATIO'] *0.89
sum_wells_1

sum_wells_1['%IS'] = sum_wells_1['%IS'].replace(np.inf, 'ND')
sum_wells_1['RATIO'] = sum_wells_1['RATIO'].replace(np.inf, 'ND')
sum_wells_1['%RATIO'] = sum_wells_1['%RATIO'].replace(np.inf, 'ND')
#real['Sample'].value_counts()


# getting average of bcrabl for all samples
bcrabl_mean = full.groupby('Sample', sort=False, as_index=False)['BCR-ABL1'].mean()
bcrabl_mean.rename(columns={'Sample': 'SAMPLE', 'BCR-ABL1':'BCR_ABL1_average'}, inplace=True)
#bcrabl_mean
sum_wells_mean = pd.merge(bcrabl_mean,sum_wells_1, on = 'SAMPLE' )
#sum_wells_mean
## counting no. of wells per sample after outliers detected
count_wells = pd.DataFrame(sum_wells.Sample.value_counts().reset_index())
#count_wells
count_wells.columns = ['SAMPLE', 'count']
sum_wells_mean_counts = pd.merge(sum_wells_mean, count_wells, on='SAMPLE')
#print(counts)
sum_wells_mean_counts

## picking up right negative control value
negvalue = sum_wells_mean_counts.loc[sum_wells_mean_counts['SAMPLE'] == "Negative", 'BCR_ABL1_average']
negvalue

## function for correcting BCR-ABL1 copies
def calculate_correction(row):
    return row['BCR-ABL1'] - row['count'] * negvalue

sum_wells_mean_counts['corrected_BCR_ABL1_COPIES'] = sum_wells_mean_counts.apply(calculate_correction, axis=1)
sum_wells_mean_counts

##getting ratios for corrected bcr-abl1 copies
sum_wells_mean_counts['RATIO'] = sum_wells_mean_counts['corrected_BCR_ABL1_COPIES'] / sum_wells_1['ABL1']
sum_wells_mean_counts['%RATIO'] = sum_wells_mean_counts['RATIO'] *100
sum_wells_mean_counts['%IS'] = sum_wells_mean_counts['%RATIO'] *0.89

#sum_wells_mean_counts['%IS'] = sum_wells_mean_counts['%IS'].replace(np.inf, 'ND')
#sum_wells_mean_counts['RATIO'] = sum_wells_mean_counts['RATIO'].replace(np.inf, 'ND')
#sum_wells_mean_counts['%RATIO'] = sum_wells_mean_counts['%RATIO'].replace(np.inf, 'ND')
#sum_wells_mean_counts

## selecting and renaming required columns
sum_wells_mean_counts_1 = sum_wells_mean_counts.loc[:, ('SAMPLE','BCR-ABL1','corrected_BCR_ABL1_COPIES','ABL1','%RATIO','%IS','DROPLETS')]
sum_wells_mean_counts_1
sum_wells_mean_counts_1.rename(columns={'corrected_BCR_ABL1_COPIES': 'corr.BCR-ABL1','%RATIO': 'corr.%RATIO','%IS': 'corr.%IS'}, inplace=True)
sum_wells_mean_counts_1
## setting negative value to zero for bcrabl1copies, ratios to ND
sum_wells_mean_counts_1.loc[~(sum_wells_mean_counts_1['corr.BCR-ABL1'] > 0), 'corr.BCR-ABL1'] = 0
sum_wells_mean_counts_1.loc[~(sum_wells_mean_counts_1['corr.%RATIO'] > 0), 'corr.%RATIO'] = 0
sum_wells_mean_counts_1.loc[~(sum_wells_mean_counts_1['corr.%IS'] > 0), 'corr.%IS'] = 0
sum_wells_mean_counts_1['corr.%RATIO'] = sum_wells_mean_counts_1['corr.%RATIO'].fillna('ND')
sum_wells_mean_counts_1['corr.%IS'] = sum_wells_mean_counts_1['corr.%IS'].fillna('ND')
## replacing inf to ND
sum_wells_mean_counts_1['corr.%IS'] = sum_wells_mean_counts_1['corr.%IS'].replace(np.inf, 'ND')
#sum_wells_mean_counts['RATIO'] = sum_wells_mean_counts['RATIO'].replace(np.inf, 'ND')
sum_wells_mean_counts_1['corr.%RATIO'] = sum_wells_mean_counts_1['corr.%RATIO'].replace(np.inf, 'ND')

backgroundvalue = sum_wells_mean_counts_1.loc[sum_wells_mean_counts_1['SAMPLE'] == "Negative", 'BCR-ABL1']
backgroundvalue
# get the wells with low droplet count
def low_droplets(x):
    if x < 12000:
        return x
    else:
        return '.'

real['DROPLETS<12000'] = real['Droplets'].apply(low_droplets)
real
# get the raw individual wells data
individual_wells = real.copy()
individual_wells.rename(columns={'Well':'WELL', 'Sample':'SAMPLE','Droplets':'DROPLETS','Ratio':'RATIO','%Ratio':'%RATIO','test_outlier':'OUTLIERWELL'},inplace=True)
individual_wells_1 = individual_wells[['WELL','SAMPLE', 'BCR-ABL1', 'ABL1','%RATIO','%IS','DROPLETS','OUTLIERWELL', 'DROPLETS<12000']]
#print(individual_wells_1)

###individual_wells_2 = individual_wells_1.round({'BCR_ABL1':6,'ABL1':6,'RATIO':6,'%Ratio':6})
###individual_wells_2

###sum_wells_3 = sum_wells_1.round({'BCR_ABL1_COPIES':6,'ABL1_COPIES':6,'RATIO':6, '%Ratio':6, '%IS':6})
###sum_wells_3

#sum_wells_3['BCR_ABL1_COPIES'] = sum_wells_3['BCR_ABL1_COPIES'].astype('int')
#sum_wells_3['ABL1_COPIES'] = sum_wells_3['ABL1_COPIES'].astype('int')
#sum_wells_3
#print(sum_wells_3)
#sum_wells_4 = sum_wells_3.reset_index(True)

#sum_wells_1['SAMPLE_NO'] = sum_wells_1.index+1
#sum_wells_1.assign(C=" ")
#sum_wells_fin = sum_wells_1.loc[:, ('SAMPLE_NO','SAMPLE', 'BCR-ABL1', 'ABL1', '%RATIO', '%IS', 'DROPLETS','C')]
###sum_wells_fin = sum_wells_1.loc[:, ('SAMPLE', 'BCR-ABL1', 'ABL1', '%RATIO', '%IS', 'DROPLETS')]
#print(sum_wells_fin)
##merging both frames to get out output
###sum_wells_final = pd.merge(sum_wells_fin,sum_wells_mean_counts_1, on = 'SAMPLE')

time_stamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S%z')
list_4 = ['Date and Time:']
list_5 = ['ddPCR pipeline version V.04.2018.02']
list_6 = ['%IS correction Parameter:0.89 ']
list_7 = ['Negative control BCR-ABL1 background value:'+ str(backgroundvalue.values).strip('[]')]

list_1 = ['Summary Results:', '', '',str(backgroundvalue.values).strip('[]'), '', '',str(0.89)]
list_2 = ['Detected Outlier Values + or -2 std:']
list_3 = ['Individual wells data:']

sum_wells_mean_counts_1.insert(loc=0, column='', value="")

#print(''+sum_wells_fin.columns)

#print(sum_wells_fin.values)

# Writing the data to output csv file
#with open('v3_Ratios_outlier_nocount.csv', 'w') as f:
with open(infile + "_out.csv", 'w') as f:
    linewriter = csv.writer(f, lineterminator='\n')
    linewriter.writerow(list_1)
    linewriter.writerow(sum_wells_mean_counts_1.columns)
    linewriter.writerows(sum_wells_mean_counts_1.values)
    linewriter.writerow(" ")
    linewriter.writerow(" ")
    linewriter.writerow(list_2)
    linewriter.writerow(outlier_table.columns)
    linewriter.writerows(outlier_table.values)
    linewriter.writerow(" ")
    linewriter.writerow(" ")
    linewriter.writerow(list_3)
    linewriter.writerow(individual_wells_1.columns)
    linewriter.writerows(individual_wells_1.values)
    linewriter.writerow(" ")
    linewriter.writerow(" ")
    linewriter.writerow(list_5)
    linewriter.writerow(list_4 + [time_stamp])
    linewriter.writerow(list_6)
    linewriter.writerow(list_7)

os.system("mv "+infile.replace(" ", "\ ")+"_out.csv "+outfolder.replace(" ", "\ "))
print(outfolder.replace(" ", "\ "))

print("*Analysis is completed*")

#>                     ddpcr plate
#>                    -------------
#>             Dataset name : small
#>             Data summary : 5 wells; 72,727 drops
#>               Plate type : ddpcr_plate
#> analysis steps : REMOVE_FAILURES, REMOVE_OUTLIERS, REMOVE_EMPTY
