#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:00:30 2020

@author: fatemehzahed
"""
import DNA_DECODING as DNA_fun
from time import time


t0 = time()

# Create Object from DNA decoding 
Object = DNA_fun.DNA_DECODING('data_scientist_exercise_file.fq',range(6))

# Get a list of all DNAs
data = Object.get_data()

# Create a dataframe of pairs with four columns:
#    (Barcode1, NumberB1, Barcode2, NumberB2)
df = Object.count_barcodes(data)

# Save the dataframe as .csv
Object.save_pairs_to_csv(df)

# Create another dataframe of non-pairs and save as .csv
df2 = Object.save_unpairs_to_csv(df)

print('Head of dataframe of pairs with counts: \n')
print(df.head())
print('\n')
print('Verify the dataframe with the document: \n')
print(df.loc[df.Barcode1 == 'TATACA'])
print('\n')
print('Head of dataframe of unpairs with counts: \n')
print(df2.head())

t1 = time()
print('\n execution time = ',t1-t0)
