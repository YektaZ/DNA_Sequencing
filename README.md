## DNA_Decoding_Function

### Problem Statement
Sequencing often involves decoding a DNA fragment in two directions: forward (+) and reverse (-). 
Every individual fragment is identified by a unique “barcode”, 
which is used to pair forward and reverse reads coming from the same fragment. 
The task is to implement an algorithm to extract pairs of reads from a FASTQ file containing sequencing results.

### Problem Specifications
1 - Open the files "OpenMe.py" and "DNA_DECODING.py" to find the code for the algorithm. The algorithm is written in Python3 and calls "pandas" library. 

2 - Open the files "paires.csv" and "unpaired.csv" to find the output of the algorithm.

3 - The time complexity of the algorithm is O(NlogN) when offset = 0 and O(N^3) when offset != 0. 

  The space complexity of the algorithm is O(N+ N/4 + 5N/4 + 2N/4 + etc.) = O(3N + etc.).
  Therefore, the space complexity of the algorithm changes linearly with the size of the input. 
  
  Usually there is a trad-off between time and space complexity. In case of larger data sets (hundreds of millions of DNA) maybe we could alter the algorithm to less time complex but higher complexity in space. For example by using hash tables. 

4 - A variable "offset" is defined inside the function. This value determines how many letter differences can exist in the barcode search.
  For example offset = 0 means that the search for barcode "XXXXXX" finds the exact same barcode. 
  However, when offset = 1, the algorithm will detect "XXXXXX" and any other barcode with only 1 letter mismatch. 
  If offset is greater than zero, the time complexity of the algorithm increases to O(N^3). 
  
5 - The is a special case of barcode: "AAATTT" where the reverse barcode is similar to the original barcode. 
This special case has been categorized in the unpaired category, however, further consideration can be taken into account for this special case.
