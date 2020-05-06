#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:02:25 2020

@author: fatemehzahed

This module implements an algorithm to extract pairs of reads from a FASTQ 
    file containing sequencing results.
    
Input: file_name, base(range of base barcode)
Outputs: 'pairs.csv' and 'unpaired.csv'

"""


class DNA_DECODING():
    
    
    def __init__(self, file_name, base):
        self.file_name = file_name
        self.base = base    # Input range
        self.offset = 0     # Number of mismatches
    
    
    def get_data(self):
        """ Read the DNA string from FASTQ file """
        file_name = self.file_name
        with open(file_name) as f:
            lines=f.readlines()
        data=[item[:-1] for item in lines[1::4]]
        return data
    
    
    def find_reverse(self, DNA):
        """ Get DNA and return barcode and reverse complement """
        base = self.base
        dic = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        barcode1 = ''
        for i in base:
            barcode1 += DNA[i]
        barcode2 = ''
        for i in barcode1:
            barcode2 += dic[i]
        return barcode1, barcode2[::-1]
    
    
    def count_barcodes(self, data):
        """ Count the number of barcode accurance and return a dataframe """
        import pandas as pd
        
        offset = self.offset

        columns = {'Barcode1', 'NumberB1', 'Barcode2', 'NumberB2'}
        df = pd.DataFrame(index=range(len(data)), columns=columns)
        df.NumberB1 = df.fillna(0)
        df.NumberB2 = df.fillna(0)

        counter = 0

        for DNA in data:
            bar1, bar2 = self.find_reverse(DNA)
            
            if not offset:  # No mismatch

                if sum(df.Barcode1 == bar1):  # If bar1 has been called before
                    df.loc[(df.Barcode1 == bar1),'NumberB1'] += 1
                elif sum(df.Barcode2 == bar1): # If bar1 has been called before
                    df.loc[(df.Barcode2 == bar1),'NumberB2'] += 1
                else: # New barcode
                    df.loc[counter, 'Barcode1'] = bar1
                    df.loc[counter,'Barcode2'] = bar2
                    df.loc[counter,'NumberB1'] += 1
            else:   # With mismatch
                
                if sum(self.barcode_comparison(bar1,df.Barcode1)):
                    df.loc[(self.barcode_comparison(bar1,df.Barcode1)),
                           'NumberB1'] += 1
                elif sum(self.barcode_comparison(bar1,df.Barcode2)): 
                    df.loc[(self.barcode_comparison(bar1,df.Barcode2)),
                           'NumberB2'] += 1
                else:
                    df.loc[counter, 'Barcode1'] = bar1
                    df.loc[counter,'Barcode2'] = bar2
                    df.loc[counter,'NumberB1'] += 1

            counter +=1
        
        df = df.dropna()
        df = self.reorder(df,data)
        
        return df
    
    
    def reorder(self, df, data):
        """ Order by minimum coverage (NumberB1<NumberB2) """
        for i in range(len(data)):
            try:
                if df.NumberB1[i] > df.NumberB2[i]:
                    df.loc[i,'NumberB2'], df.loc[i,'NumberB1'] = (
                        df.loc[i,'NumberB1'], df.loc[i,'NumberB2'])
                    df.loc[i,'Barcode2'], df.loc[i,'Barcode1'] = (
                        df.loc[i,'Barcode1'], df.loc[i,'Barcode2'])
            except:
                continue
        return df
    
    
    def barcode_comparison(self, bar1, bars):
        """ 
        This function counts the number of mismatches between barcodes, 
        returns a string of booleans
        """
        offset = self.offset
        difference = []
        for b in bars:
            if not (b != b): # b is not NAN
                temp = [bar1[i] for i in range(len(bar1)) if bar1[i] != b[i]]
                if len(temp)>offset:
                    difference.append(False)
                else:
                    difference.append(True)
            else:
                difference.append(False)
        return difference
    
    
    def save_pairs_to_csv(self, df):
        """ Save to .csv """
        df = df[['Barcode1', 'NumberB1', 'Barcode2', 'NumberB2']]
        df.to_csv('pairs.csv', index=False)
      
        
    def save_unpairs_to_csv(self, df):
        """ Create a dataframe for unpairs and save as .csv """
        df2 = df[df.NumberB1 == 0]
        df2 = df2.sort_values(by=['NumberB2'], ascending=False)
        df2 = df2.drop(['Barcode1', 'NumberB1'], axis = 1)
        df2.rename(columns={'Barcode2': 'barcode', 'NumberB2': 'numbers'},
                   inplace=True)
        df2 = df2[['barcode', 'numbers']]
        df2.to_csv('unpaired.csv', index=False)
        return df2

