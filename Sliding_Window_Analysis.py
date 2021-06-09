########################################################################################################################################################
## July 2019 Michelle Amaral
## This script performs burden analysis using a sliding window by counting the number of qualifying variants within specifically defined window lengths:
## 500, 1000, 5000, 10000, 25000, 50000, 1000000
## 
## To run: python Sliding_Window_Analysis.py filtered_Variant_Counts_for_Chr_22_CADD_15.0_MAF_0.005 chrNum  sizeOfSlide outputDirectory logFile
## 
## Example for chromosome 22, 50000 bp windows, sliding 500 bp:
## python testing_Sliding_Window_Analysis.py filtered_Variant_Counts_for_Chr_22_CADD_15.0_MAF_0.005 22 500 . log_output_Chr_22_CADD_15.0_MAF_0.005.log
## 
## output file contains the name of the previously filtered input file, for example: filtered_Variant_Counts_for_Chr_22_CADD_15.0_MAF_0.005 22 500
## 
########################################################################################################################################################

import sys,re,os
import pandas as pd
import numpy as np
#import scipy.stats as stats

inputFileName = sys.argv[1]
splitFileName = inputFileName.strip().split('_')
CADD = splitFileName[7] # from the input file name, parse CADD used for filtering  
MAF = splitFileName[9] # from the input file name, parse MAF used for filtering 

chromosomeNumber = sys.argv[2]
#windowSize = int(sys.argv[3])     ********************
slideSize = int(sys.argv[3])

OUTPUT_directory = sys.argv[4]
#main_results_output = 'Sliding_Window_Analysis_for_Chr_' + chromosomeNumber + '_Window_' + str(windowSize) + 'bp_SlideSize_' + str(slideSize) ******
#print(main_results_output) *****

# this argument is the file name of the log file from filtering script "variant_dictionary.py"
logOutputFile = sys.argv[5]


###########################################################################################################
## function to obtain length of chromosome from Chromosome_Lengths.txt
## copied & pasted from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
###########################################################################################################

def locate_chromosome_length(chromosomeNumber):
  with open('/gpfs/gpfs1/home/mamaral/Pritzker_Whole_Genome/2019_batch_call/Chromosome_Lengths.txt', "r") as chrLens:
    for line in chrLens:
      split_chr_file = line.strip().split('\t')
      obtainChrNum = split_chr_file[0].split('chr')
      chrNumTemp = obtainChrNum[1]
      if chrNumTemp == chromosomeNumber:
        thisChrLength = split_chr_file[1]
        pass
      else:
        continue

    return(thisChrLength)


###########################################################################################################
## function to generate chromosome ranges of n base pair window size along a chromosome 
## of a certain length
## {(lower end, upper end): {'Control': [0], 'Case': [0]}}
###########################################################################################################

def generate_chrom_ranges(this_chr_length, window_size, slide_size): 
  chr_ranges = []
  for i in range(1, int(this_chr_length)):    
    if ( i % slide_size == 0 and (i - (window_size - 1)) > 0 ): 
      chr_ranges.append((i - (window_size - 1), i))
    else:
      continue

  chrRange_count_dictionary = {}
  for k in chr_ranges:
    chrRange_count_dictionary[k] = {}
    chrRange_count_dictionary[k]['Case'] = [0]
    chrRange_count_dictionary[k]['Control'] = [0]
    
  return(chrRange_count_dictionary)


###########################################################################################################
## function sums all control variant counts and all case variant counts for each chromosome range in the 
## dictionary
###########################################################################################################

def add_all_counts(results_dictionary):
  for pos_range in results_dictionary:
    finalTallyControls = (sum(results_dictionary[pos_range]['Control']))
    results_dictionary[pos_range]['Control'] = finalTallyControls

    finalTallyCases = (sum(results_dictionary[pos_range]['Case']))
    results_dictionary[pos_range]['Case'] = finalTallyCases

  return(results_dictionary)

###########################################################################
## function to obtain total number of control subjects
## insert path to subject key
###########################################################################

def get_number_controls():
  controlCounter = 0
  with open('', "r") as subjectKey:
    for line in subjectKey:
      split_subjectKey = line.strip().split('\t')
      if split_subjectKey[0].startswith('Use'):
        continue
      else:
        status = split_subjectKey[9]
        if status == 'Control':
          controlCounter = controlCounter + 1
    #print(controlCounter)
    return (controlCounter)


###########################################################################
## function to obtain total number of case (affected) subjects
## insert path to subject key
###########################################################################

def get_number_cases():
  caseCounter = 0
  with open(' ', "r") as subjectKey:
    for line in subjectKey:
      split_subjectKey = line.strip().split('\t')
      if split_subjectKey[0].startswith('Use'):
        continue
      else:
        status = split_subjectKey[9]
        if status == 'Case':
          caseCounter = caseCounter + 1
    #print(caseCounter)
    return (caseCounter)


###########################################################################################################
## BODY
## Insert path to log file and the input file
###########################################################################################################
# want the variantCount variable to be global
variantCount = 0

# create global variant for total control subjects
totalControls = get_number_controls()

# create global variant for total case subjects
totalCases = get_number_cases()

# open the log file, which is the output from the filtering script "variant_dictionary.py"
# extract the count of total variants from the file
with open (' ' + logOutputFile, 'r') as logOutput:
  for line in logOutput:
    split_the_lines = line.strip().split(' ')
    if split_the_lines[0].startswith('#After'):
      variantCount = split_the_lines[4]
    else:
      continue
  
chromosomeLength = locate_chromosome_length( chromosomeNumber )
theWindowSizes = [500, 1000, 5000, 10000, 25000, 50000, 100000]

for windowSize in theWindowSizes:
  with open (' ', 'r') as countTable:
    countLinesProcessed = 0
    numberOfTests = 0
    main_results_output = 'data_Sliding_Window_Analysis_for_Chr_' + chromosomeNumber + '_Window_' + str(windowSize) + 'bp_SlideSize_' + str(slideSize) + '_CADD_' + CADD + '_MAF_' + MAF
    print(main_results_output)

    rangeDict = generate_chrom_ranges(chromosomeLength, windowSize, slideSize)

    for line in countTable:

      if countLinesProcessed % 1000 == 0:
          print('Still working on ', windowSize, ' bp window size. I am at line: ', countLinesProcessed, ' of ', variantCount)
          print('This is ', (int(countLinesProcessed)/int(variantCount))*100, '% ', 'complete' )

      split_line = line.strip().split('\t')
      # skip header
      if split_line[0].startswith('POSITION'):
        continue
      else:
        countLinesProcessed = countLinesProcessed + 1
        #print(countLinesProcessed)
        # the list below is created for every data line that is processed
        # it will hold a list of the ranges that each chromosome position fits into
        rangesItBelongsTo = []
        
        # extended position format: chrNumber_chrPosition_Ref_Alt
        extendedPosition = split_line[0]
        deconstructPosition = extendedPosition.strip().split('_')
        chromNumExtendPos = deconstructPosition[0]
        chrPosition = int(deconstructPosition[1])

        # QC check
        if chromNumExtendPos != chromosomeNumber:
          print('Chromosome number does not match')

        # control counts are in the first column after the extended position
        controlCounts = int(split_line[1])
        # case counts are in the second column after the extended position
        caseCounts = int(split_line[2])

        for this_range in rangeDict:
          if chrPosition >= this_range[0] and chrPosition <= this_range[1]:
            rangesItBelongsTo.append(this_range)
          #else:
            #continue
            #doesNotBelongHere = doesNotBelongHere + 1
            #print('does not belong here')
        for k in rangesItBelongsTo:
          rangeDict[k]['Control'].append(controlCounts)
          rangeDict[k]['Case'].append(caseCounts) 

    # after all lines of the file have been processed add together the total number of variant counts per range
    totalsDict = add_all_counts(rangeDict)
    #df = pd.DataFrame.from_dict(totalsDict,orient='index')
    #df.to_csv()

    with open (main_results_output, 'x') as mainOutput:
      mainOutput.write('#This is an analysis of the data from: ' + inputFileName + '\n')
      mainOutput.write('#START' + '\t' + 'STOP' + '\t' + 'CONTROL' + '\t' + 'CASE' + '\n')
      for posRange in totalsDict:
        if (totalsDict[posRange]['Control'] != 0 and totalsDict[posRange]['Case'] != 0):
          rangeStartStop = ''.join(str(posRange))
          mainOutput.write(rangeStartStop + '\t' + str(totalsDict[posRange]['Control']) + '\t' + str(totalsDict[posRange]['Case']) + '\n')









######################################## Exercise using code below

#    with open (main_results_output, 'x') as mainOutput:
#      mainOutput.write('#This is an analysis of the data from: ' + inputFileName + '\n')
#      mainOutput.write('#START' + '\t' + 'STOP' + '\t' + 'CONTROL' + '\t' + 'CASE' + '\t' + 'PVAL' + '\n')
#      for posRange in totalsDict:
#        if (totalsDict[posRange]['Control'] != 0 and totalsDict[posRange]['Case'] != 0):
#          numberOfTests = numberOfTests + 1
#          obs = np.array([[totalsDict[posRange]['Control'], totalsDict[posRange]['Case']], [totalControls - totalsDict[posRange]['Control'], totalCases - totalsDict[posRange]['Case']]])
#          fisher_result = stats.fisher_exact(obs)
#          pval = fisher_result[1]
          #may want to add something here to distinguish places in which 
          #case counts are greater than controls? would decrease the number
          #of tests performed
#          if pval > 0.05:
#            continue
#          else: 
#            rangeStartStop = ''.join(str(posRange))
#            mainOutput.write(rangeStartStop + '\t' + str(totalsDict[posRange]['Control']) + '\t' + str(totalsDict[posRange]['Case']) + '\t' + str(pval) + '\n')
#      mainOutput.write('#The number of fisher exact tests conducted was: ' + str(numberOfTests))

  #for thing in totalsDict:
  #  if (totalsDict[thing]['Control'] != 0 and totalsDict[thing]['Case'] != 0):
  #    print(thing, totalsDict[thing])

  #for thing in totalsDict:
  #  if (totalsDict[thing]['Control'] != 0 and totalsDict[thing]['Case'] != 0):
      #print(thing, totalsDict[thing])
  #    oddsratio, pvalue = stats.fisher_exact(table)


#obs = np.array([[totalsDict[posRange]['Control'], totalsDict[posRange]['Case']], [totalControls - totalsDict[posRange]['Control'], totalCases - totalsDict[posRange]['Case']]])
#fisher_result = stats.fisher_exact(obs)
#pval = fisher_result[1]
#if pval >

#[[       totalsDict[posRange]['Control'], totalsDict[posRange]['Case']         ]
# [       totalControls - totalsDict[posRange]['Control'], totalCases - totalsDict[posRange]['Case']      ]]
















