# coding: utf-8

# In[76]:

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from collections import OrderedDict
from datetime import date
from ipykernel import kernelapp as app
from datetime import datetime
from pandas import rolling_median
from io import StringIO
from matplotlib import style

infile = sys.argv[1]
outfolder = "/ddPCR BCR-ABL/Experiments and Raw Data/Analyzed_results_ddpcr/ddPCR_Results_version3.1"
print(outfolder.replace(" ", "\ "))
print(infile)
data = pd.read_excel(sys.argv[1], sheet_name='Sheet')


print("Testing....")

#data.head()
#data.columns
data1 = data.loc[:, ['Well', 'Sample','Target', 'Conc(copies/µL)', 'Accepted Droplets']]
data1
data1.rename(columns={'Conc(copies/µL)': 'Conc', 'Accepted Droplets' : 'Droplets'}, inplace=True)
data1


# In[77]:

wt = data1.loc[data1['Target'] == "HEX"]
mut = data1.loc[data1['Target'] == "FAM"]
wt = wt.reset_index(drop=True)
#wt
mut = mut.reset_index(drop=True)
#mut
wildtype = wt.drop('Target', 1)
wildtype
wildtype.rename(columns={'Conc':'ABL1'}, inplace=True)
mutant = mut.drop('Target', 1)
mutant
mutant.rename(columns={'Conc':'BCR_ABL1' }, inplace=True)
mutant
new = pd.concat([wildtype,mutant['BCR_ABL1']], axis =1)
new
new['Ratio'] = new['BCR_ABL1'] / new['ABL1']
new



# In[78]:

new['Ratio'] = new['Ratio'].fillna(0)
new['Ratio'] = new['Ratio'].replace('inf', 0)
new['%_Ratio'] = new['Ratio']*100
new


# In[79]:

full = new.copy()
full
sample = full.groupby(['Sample'])
sample


# In[80]:

def s(one):
    one['%_Ratio_Mean'] = one['%_Ratio'].mean()
    one['2*std'] = 2*one['%_Ratio'].std()
    one['notaOutlier'] = (one['%_Ratio'] == 0) | ((one['%_Ratio'] > (one['%_Ratio'].mean() - 2*one['%_Ratio'].std())) & (one['%_Ratio'] < (one['%_Ratio'].mean() + 2*one['%_Ratio'].std())))
    return one

real = sample.apply(s)
real


# In[81]:

def testing_outlier(row):
    if row['notaOutlier'] == True:
        return '.'
    else:
        return 'Yes'
real = real.assign(test_outlier = real.apply(testing_outlier, axis=1))
real


# In[82]:

# selecting wells with outliers / seperating from table to get the outliers table
real_1_out = real.loc[real['notaOutlier'] == False]
real_1_out
real_1_out.rename(columns={ 'Well': 'WELL', 'Sample': 'SAMPLE_ID','Ratio':'RATIO','%_Ratio':'%_RATIO','Droplets':'DROPLETS','test_outlier':'OUTLIER_WELL'}, inplace=True)
real_1_out



# In[83]:

# Outlier table
outlier_table = real_1_out[['WELL','SAMPLE_ID','BCR_ABL1','ABL1','RATIO','%_RATIO','DROPLETS','OUTLIER_WELL']]
outlier_table_1 = outlier_table.round({'RATIO':4,'%_RATIO':4})
outlier_table_1['BCR_ABL1'] = outlier_table_1['BCR_ABL1'].round().astype('int')
outlier_table_1['ABL1'] = outlier_table_1['ABL1'].round().astype('int')
outlier_table_1


# In[85]:

# selecting wells with no outliers, gropuing and getting the sum
sum_wells = real.loc[real['notaOutlier'] == True]
sum_wells
sum_wells_1 = sum_wells.groupby(['Sample'], sort=False, as_index=False).agg({'BCR_ABL1': np.sum, 'ABL1': np.sum, 'Droplets': np.sum })
sum_wells_1
sum_wells_1.rename(columns={'Sample': 'SAMPLE_ID','BCR_ABL1':'BCR_ABL1_COPIES', 'ABL1':'ABL1_COPIES','Droplets':'DROPLETS_SUM'}, inplace=True)
sum_wells_1


# In[50]:

#getting ratios and % ratios for the sum values.
sum_wells_1['RATIO'] = sum_wells_1['BCR_ABL1_COPIES'] / sum_wells_1['ABL1_COPIES']
sum_wells_1['%_RATIO'] = sum_wells_1['RATIO'] *100
sum_wells_1


#real['Sample'].value_counts()


# In[86]:

def low_droplets(x):
    if x < 12000:
        return x
    else:
        return '.'

real['DROPLETS<12000'] = real['Droplets'].apply(low_droplets)
real




# In[87]:

individual_wells = real.copy()
individual_wells.rename(columns={'Well':'WELL', 'Sample':'SAMPLE_ID','Droplets':'DROPLETS','Ratio':'RATIO','%_Ratio':'%_RATIO','test_outlier':'OUTLIER_WELL'},inplace=True)
individual_wells_1 = individual_wells[['WELL','SAMPLE_ID', 'BCR_ABL1', 'ABL1','RATIO','%_RATIO','DROPLETS','OUTLIER_WELL', 'DROPLETS<12000']]
#individual_wells_1


# In[88]:

individual_wells_2 = individual_wells_1.round({'RATIO':4,'%_RATIO':4})
individual_wells_2


# In[89]:

individual_wells_2['BCR_ABL1'] = individual_wells_2['BCR_ABL1'].round().astype('int')
individual_wells_2['ABL1'] = individual_wells_2['ABL1'].round().astype('int')

sum_wells_3 = sum_wells_1.round({'RATIO':4, '%_RATIO':4})
sum_wells_3


# In[93]:

sum_wells_3['BCR_ABL1_COPIES'] = sum_wells_3['BCR_ABL1_COPIES'].astype('int')
sum_wells_3['ABL1_COPIES'] = sum_wells_3['ABL1_COPIES'].astype('int')
sum_wells_3


# In[94]:
#print(sum_wells_3)
#sum_wells_4 = sum_wells_3.reset_index(True)

sum_wells_3['SAMPLE_NO'] = sum_wells_3.index+1


# In[95]:

#sum_wells_5 = sum_wells_3.drop('index', 1)
#sum_wells_5

sum_wells_fin = sum_wells_3.loc[:, ('SAMPLE_NO','SAMPLE_ID', 'BCR_ABL1_COPIES', 'ABL1_COPIES', 'RATIO', '%_RATIO', 'DROPLETS_SUM')]
sum_wells_fin


# In[96]:

time_stamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S%z')
list_4 = ['Date and Time:']
list_5 = ['version: ddpcr v.3.1']


# In[98]:

list_1 = ['Summary Results:']
list_2 = ['Detected Outlier Values + or -2 std:']
list_3 = ['Individual wells data:']

#with open('v3_Ratios_outlier_nocount.csv', 'w') as f:
with open(infile + "_out.csv", 'w') as f:
    linewriter = csv.writer(f, lineterminator='\n')
    linewriter.writerow(list_1)
    linewriter.writerow(sum_wells_fin.columns)
    linewriter.writerows(sum_wells_fin.values)
    linewriter.writerow(" ")
    linewriter.writerow(" ")
    linewriter.writerow(list_2)
    linewriter.writerow(outlier_table_1.columns)
    linewriter.writerows(outlier_table_1.values)
    linewriter.writerow(" ")
    linewriter.writerow(" ")
    linewriter.writerow(list_3)
    linewriter.writerow(individual_wells_2.columns)
    linewriter.writerows(individual_wells_2.values)
    linewriter.writerow(" ")
    linewriter.writerow(" ")
    linewriter.writerow(list_4)
    linewriter.writerow(list_5 + [time_stamp])

os.system("mv "+infile.replace(" ", "\ ")+"_out.csv "+outfolder.replace(" ", "\ "))
print(outfolder.replace(" ", "\ "))

print("*Analysis is completed*")

#>                     ddpcr plate
#>                    -------------
#>             Dataset name : small
#>             Data summary : 5 wells; 72,727 drops
#>               Plate type : ddpcr_plate
#> analysis steps : REMOVE_FAILURES, REMOVE_OUTLIERS, REMOVE_EMPTY


