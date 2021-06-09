########################################################################################################################################################
## August 2019 Michelle Amaral
## This script performs gene burden analysis. Each subject is restricted to two qualifying variants per window
## 
## To run: python protein_altering_gene_burden_2_variants_per_interval.py pickledDictionary chrNum outputDirectory
## 
## Example for chromosome 22:
## python protein_altering_gene_burden_2_variants_per_interval.py big_variant_dict_for_Chr_22_CADD_15.0_MAF_0.005.p 22 . 
##
##
########################################################################################################################################################

import sys,re,os
import pandas as pd
import numpy as np
import pickle

inputFileName = sys.argv[1]
splitFileName = inputFileName.strip().split('_')
CADD = splitFileName[7] # from the input file name, parse CADD used for filtering  
splitMAF = splitFileName[9].strip().split('.p')
MAF = splitMAF[0] # from the input file name, parse MAF used for filtering 

chromosomeNumber = sys.argv[2]

OUTPUT_directory = sys.argv[3]

def generate_gene_dict(geneRanges, nameOfGenes):

  geneRange_count_dictionary = {}
  for k in list(range(0, len(geneRanges), 1)):
    empty_ids = get_case_control_ids()
    geneRange_count_dictionary[geneRanges[k]] = empty_ids
    geneRange_count_dictionary[geneRanges[k]]['Case'] = [0]
    geneRange_count_dictionary[geneRanges[k]]['Control'] = [0]
    geneRange_count_dictionary[geneRanges[k]]['Gene_Name'] = nameOfGenes[k]
  
  return(geneRange_count_dictionary)

#############################################################################################################################################
## Other info about the UCSC table browser output: "Genes and Gene Prediction", Track: "NCBI RefSeq", Table: "RefSeq Curated", region:
## "genome"
#############################################################################################################################################

def process_gene_coordinates(chromosomeNumber):
  with open('Gene_coordinates.txt', "r") as annoFile:
    chr_ranges = []
    gene_names = []
    gene_lengths = []
    for line in annoFile:
      if line.startswith('#'):
        continue
      else:
        chrNumList = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
        annoChrNum = line.strip().split('\t')[0]
        chrNum = annoChrNum.strip().split('chr')[1]
        rangeStart = int(line.strip().split('\t')[1]) - 50000
        rangeStop = int(line.strip().split('\t')[2]) + 50000
        geneName = line.strip().split('\t')[3]
        geneLength = rangeStop - rangeStart
     
        if not chrNum in chrNumList:
          continue
        if int(chrNum) != int(chromosomeNumber):
          continue
        else:
          chr_ranges.append((rangeStart, rangeStop))
          gene_names.append(geneName)
          gene_lengths.append(geneLength)
    
    print('length of chr ranges: ', len(chr_ranges))
    print('length of gene_names: ', len(gene_names))

    final_chr_ranges = []
    final_gene_names = []
    final_gene_lengths = []

    i = 0

    while i < (len(gene_names) - 1):

      if (i+1) == len(gene_names):
        final_chr_ranges.append(final_chr_range)
        final_gene_names.append(final_gene_name)

      else:
        final_gene_length = gene_lengths[i]
        final_chr_range = chr_ranges[i]
        final_gene_name = gene_names[i]
        #print('This is chr_ranges[i]: ', chr_ranges[i])
        #print('Before inner while loop. Final gene length is: ', gene_lengths[i])
        #print('Before inner while loop. Final chr range is: ', final_chr_range)
        #print('Before inner while loop. Final gene name is: ', gene_names[i])
        while gene_names[i] == gene_names[(i+1)]:
          #print('Yes they equal. Gene_names[i] is: ', gene_names[i], 'Gene_names[i+1] is: ', gene_names[i+1])
          if gene_lengths[i] < gene_lengths[(i+1)]:
            #print('They equal and gene_lengths[i] ', gene_lengths[i], ' is less than ', gene_lengths[(i+1)])
            final_gene_length = gene_lengths[(i+1)]
            final_chr_range = chr_ranges[(i+1)]
          i+=1
          #print('i after increment is: ', i)
          if (i+1) == len(gene_names):
            break

        final_chr_ranges.append(final_chr_range)
        final_gene_names.append(final_gene_name)
        i+=1

    emptyDict = generate_gene_dict(final_chr_ranges, final_gene_names)
    return(emptyDict)


################################################################################
## output is a dictionary containing subject ID, ethnicity, and disease status
################################################################################

def get_case_control_counts():
  subjectDict = {}
  with open('master_key.txt', "r") as subjectKey:
    for line in subjectKey:
      split_subjectKey = line.strip().split('\t')
      if split_subjectKey[0].startswith('Use'):
        continue
      else:
        subjID = split_subjectKey[0]
        ethnicity = split_subjectKey[6] 
        status = split_subjectKey[9]

        holdingDict = {}
        holdingDict['status'] = status
        holdingDict['ethnicity'] = ethnicity
        
        subjectDict[subjID] = holdingDict

    return (subjectDict)


###########################################################################
## function to obtain case / control ids
###########################################################################

def get_case_control_ids():
  subjectDict = {}
  with open('master_key_Jul_18.txt', "r") as subjectKey:
    for line in subjectKey:
      split_subjectKey = line.strip().split('\t')
      if split_subjectKey[0].startswith('Use'):
        continue
      else:
        subjID = split_subjectKey[0]        
        subjectDict[subjID] = [0]

    return (subjectDict)

###########################################################################################################
## function sums all control variant counts and all case variant counts for each chromosome range in the 
## dictionary
###########################################################################################################

def add_all_counts(results_dictionary):
  
  for pos_range in results_dictionary:
    for subject in results_dictionary[pos_range]:
      if len(results_dictionary[pos_range][subject]) > 3 or len(results_dictionary[pos_range][subject]) == 1:
        continue
      if subject == 'Case' or subject == 'Control' or subject == 'Gene_Name':
        continue      
      if subjectDict[subject]['status'] == 'Case':
        case_count = sum(results_dictionary[pos_range][subject])
        results_dictionary[pos_range]['Case'].append(case_count)
      else:
        if subjectDict[subject]['status'] == 'Control':
          control_count = sum(results_dictionary[pos_range][subject])
          results_dictionary[pos_range]['Control'].append(control_count)

    finalTallyCases = (sum(results_dictionary[pos_range]['Case']))
    results_dictionary[pos_range]['Case'] = finalTallyCases
    if finalTallyCases != 0:
      print('this is the final Tally Cases for this range: ', results_dictionary[pos_range]['Case'])

    finalTallyControls = (sum(results_dictionary[pos_range]['Control']))
    results_dictionary[pos_range]['Control'] = finalTallyControls
    if finalTallyCases != 0:
      print('this is the final Tally Controls for this range: ', results_dictionary[pos_range]['Control'])
  return(results_dictionary)


###########################################################################################################
## BODY
## 
###########################################################################################################

# generate a dictionary structured key:value with subjectID:case/control status
subjectDict = get_case_control_counts()

pickleDict = pickle.load(open(inputFileName, 'rb')) 
main_results_output = 'data_Protein_Altering_Gene_Burden_Analysis_for_Chr_' + chromosomeNumber + '_CADD_' + CADD + '_MAF_' + MAF + '_2_variants_per_interval'
print(main_results_output)

variantCountDict = process_gene_coordinates(chromosomeNumber)

# list of variant consequences
consequencesDesired = ['splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'coding_sequence_variant']

for thing in pickleDict:

  rangesItBelongsTo = []

  # extended position format: chrNumber_chrPosition_Ref_Alt
  variantPosition = int(thing.strip().split('_')[1])
  # results of consequence, sift, and polyphen will be lists
  consequence = pickleDict[thing]['Consequence']
  sift = pickleDict[thing]['SIFT']
  polyphen = pickleDict[thing]['PolyPhen']

  for this_range in variantCountDict:
    if not (variantPosition >= this_range[0] and variantPosition <= this_range[1]): 
      continue

    for j in list(range(0, len(consequence), 1)):
      if consequence[j] in consequencesDesired:
        if (consequence[j] == 'missense_variant') and (sift[j].strip().split('(')[0] == 'benign' or sift[j].strip().split('(')[0] == 'tolerated' or polyphen[j].strip().split('(')[0] == 'benign'):
          continue
        # calculate case control counts for the position
        rangesItBelongsTo.append(this_range)
        break

      else:
        continue

  for k in rangesItBelongsTo:
    badIDs = []
    for subject in pickleDict[thing]['subject_ids']:
      if subject in subjectDict:

        if pickleDict[thing]['subject_ids'][subject] == 'NA':
          continue
        if pickleDict[thing]['subject_ids'][subject] == 0:
          #print(type(pickleDict[thing]['subject_ids'][subject]))
          continue
        variantCountDict[k][subject].append(pickleDict[thing]['subject_ids'][subject])
      else:
        if subject not in badIDs:
          badIDs.append(subject)
        else:
          continue

# after all lines of the file have been processed add together the total number of variant counts per range
totalsDict = add_all_counts(variantCountDict)

with open (main_results_output, 'x') as mainOutput:
  mainOutput.write('#This is an analysis of the data from: ' + inputFileName + '\n')
  mainOutput.write('#START' + '\t' + 'STOP' + '\t' + 'CONTROL' + '\t' + 'CASE' + '\t' + 'GENE' + '\n')
  for posRange in totalsDict:
    if (totalsDict[posRange]['Control'] == 0 and totalsDict[posRange]['Case'] == 0):
      continue
    else:
      rangeStartStop = ''.join(str(posRange))
      mainOutput.write(rangeStartStop + '\t' + str(totalsDict[posRange]['Control']) + '\t' + str(totalsDict[posRange]['Case']) + '\t' + str(totalsDict[posRange]['Gene_Name']) + '\n')

