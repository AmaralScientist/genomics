#################################################################################################################################################
## Michelle Amaral
## February 2020
## The output of this script is information about variants that contributed to the counts reported from a protein-altering gene burden test. 
##
## Script needs to be run in the current directory you are working from.
##
## Input to run script is as follows:
## python Interrogate_Specific_Variants.py <txt file containing output from Fisher tests R script> <number of lines from R script output file>
## For example:
## python Interrogate_Specific_Variants.py pvals_Protein_Altering_Gene_Burden_CADD_15.0_MAF_0.005_2_variants_per_interval_ethnicity_AA.txt 10
##
#################################################################################################################################################

import sys,re,os
import pickle
import numpy as np

INPUT_pVals_file = sys.argv[1]
split_file_name = INPUT_pVals_file.strip().split('_')

if "AA.txt" in split_file_name:
  input_ethnicity = "AA"
if "EA.txt" in split_file_name:
  input_ethnicity = "EA"
if "all" in split_file_name:
  input_ethnicity = "All"

print("input ethnicity: ", input_ethnicity)

preWindowSize = split_file_name[2]
windowSize = preWindowSize.strip().split('bp')[0]

CADD = split_file_name[6]
toSplitMAF = split_file_name[8].strip().split('.')
MAF = toSplitMAF[0] + '.' + toSplitMAF[1]

# number of lines from the input pvals file is specified on the command line
number_of_lines = int(sys.argv[2])


###########################################################################
## function to obtain case / control variant counts at given position
###########################################################################

def get_case_control_counts():
  subjectDict = {}
  with open('Master_key.txt', "r") as subjectKey:
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
## body
###########################################################################

with open(INPUT_pVals_file, "r") as inputFile:
  subjectDict = get_case_control_counts()

  lineCount = 0

  for line in inputFile:
    if lineCount <= number_of_lines:
      measureMultipleVariants = []
      measureMultiplePositions = []

      if line.startswith('Chromosome'):
        continue

      split_input = line.strip().split('\t')

      controlCounts = int( split_input[3] )
      caseCounts = int( split_input[4] )

      chromosomeNumber = split_input[1]
      chromosomeRange = split_input[2]
      splitChromosomeRange = chromosomeRange.strip().split(',')

      preRangeStart = splitChromosomeRange[0]
      preRangeStop = splitChromosomeRange[1]

      rangeStart = int(preRangeStart.strip().split('(')[1])
      print("Range start: ", rangeStart)
      rangeStop = int(preRangeStop.strip().split(')')[0])   
      print("Range stop: ", rangeStop)

      if lineCount != 0:
        if rangeStart < ( store_rangeStart + int( windowSize ) ) and ( store_controlCounts == controlCounts and store_caseCounts == caseCounts ):
          continue

      main_results_output = 'parsed_sliding_window_variants_for_Chr_' + chromosomeNumber + '_Range_' + preRangeStart.strip().split('(')[1] + '_through_' + 
          preRangeStop.strip().split(')')[0] + '_' + windowSize + 'bp_Window_CADD_' + CADD + '_MAF_' + MAF + '_ethnicity_' + input_ethnicity
      print(main_results_output)

      fileName = 'big_variant_dict_for_Chr_' + chromosomeNumber + '_CADD_' + CADD + '_MAF_' + MAF + '.p'
      print(fileName)

      pickleDict = pickle.load(open(fileName, 'rb'))
      

      with open( main_results_output, 'w' ) as outputFile:
        outputFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%("Chr", "Position", "Ref", "Alt", "BatchAC", "GnomadFreq", "Gnomad3Freq", 
            "GnomadFreq_AFR", "Consequence", "AAchange", "CADD", "SubjectID", "Status", "VariantCount", "Ethnicity"))
        for thing in pickleDict:
          if int(thing.strip().split('_')[1]) <= rangeStop and int(thing.strip().split('_')[1]) >= rangeStart:
            print('CADD: ', pickleDict[thing]['CADD'])
            print('gnomAD3: ', pickleDict[thing]['gnomAD3']['gnomad3_AF_afr'][0])
            print('Batch allele count: ', pickleDict[thing]['batch_AC'])
            print(thing, pickleDict[thing]['gnomadGe_AF'])
            print(thing, pickleDict[thing]['Consequence'])
            print(thing, pickleDict[thing]['SYMBOL'])
            print('Variant Class: ', pickleDict[thing]['VARIANT_CLASS'])
            print('Feature type: ', pickleDict[thing]['Feature_type'])
            print('Exon: ', pickleDict[thing]['EXON'])
            print('Amino acids: ', pickleDict[thing]['Amino_acids'])
            print('Motif name: ', pickleDict[thing]['MOTIF_NAME'])
            print('SIFT: ', pickleDict[thing]['SIFT'])
            print('PolyPhen: ', pickleDict[thing]['PolyPhen'])

            for subject in pickleDict[thing]['subject_ids']:
              if pickleDict[thing]['subject_ids'][subject] == 'NA':
                continue
              if input_ethnicity == "AA":
                if pickleDict[thing]['subject_ids'][subject] != 0 and subject in subjectDict and subjectDict[subject]['ethnicity'] == "Black":
                  print('Subject ID: ', subject, ', Disease status: ', subjectDict[subject], ', Variant count: ', pickleDict[thing]['subject_ids'][subject])
                  outputFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chromosomeNumber, thing.strip().split('_')[1], 
                      thing.strip().split('_')[2], thing.strip().split('_')[3], pickleDict[thing]['batch_AC'], pickleDict[thing]['gnomadGe_AF'][0], 
                      pickleDict[thing]['gnomAD3']['gnomad3_AF'][0], pickleDict[thing]['gnomAD3']['gnomad3_AF_afr'][0], pickleDict[thing]['Consequence'][0], 
                      pickleDict[thing]['Amino_acids'][0], pickleDict[thing]['CADD'], subject, subjectDict[subject]['status'], 
                      pickleDict[thing]['subject_ids'][subject], subjectDict[subject]['ethnicity'] ))
                  measureMultipleVariants.append(subject)
                  measureMultiplePositions.append(thing.strip().split('_')[1])
              if input_ethnicity == "EA":
                if pickleDict[thing]['subject_ids'][subject] != 0 and subject in subjectDict and subjectDict[subject]['ethnicity'] == "White":
                  print('Subject ID: ', subject, ', Disease status: ', subjectDict[subject], ', Variant count: ', pickleDict[thing]['subject_ids'][subject])
                  outputFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chromosomeNumber, thing.strip().split('_')[1], 
                      thing.strip().split('_')[2], thing.strip().split('_')[3], pickleDict[thing]['batch_AC'], pickleDict[thing]['gnomadGe_AF'][0], 
                      pickleDict[thing]['gnomAD3']['gnomad3_AF'][0], pickleDict[thing]['gnomAD3']['gnomad3_AF_afr'][0], pickleDict[thing]['Consequence'][0], 
                      pickleDict[thing]['Amino_acids'][0], pickleDict[thing]['CADD'], subject, subjectDict[subject]['status'], 
                      pickleDict[thing]['subject_ids'][subject], subjectDict[subject]['ethnicity'] ))
                  measureMultipleVariants.append(subject)
                  measureMultiplePositions.append(thing.strip().split('_')[1])
              if input_ethnicity == "All":
                if pickleDict[thing]['subject_ids'][subject] != 0 and subject in subjectDict:
                  print('Subject ID: ', subject, ', Disease status: ', subjectDict[subject], ', Variant count: ', pickleDict[thing]['subject_ids'][subject])
                  outputFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chromosomeNumber, 
                      thing.strip().split('_')[1], thing.strip().split('_')[2], thing.strip().split('_')[3], pickleDict[thing]['batch_AC'], 
                      pickleDict[thing]['gnomadGe_AF'][0], pickleDict[thing]['gnomAD3']['gnomad3_AF'][0], pickleDict[thing]['gnomAD3']['gnomad3_AF_afr'][0], 
                      pickleDict[thing]['Consequence'][0], pickleDict[thing]['Amino_acids'][0], pickleDict[thing]['CADD'], subject, subjectDict[subject]['status'], 
                      pickleDict[thing]['subject_ids'][subject], subjectDict[subject]['ethnicity'] ))
                  measureMultipleVariants.append(subject)
                  measureMultiplePositions.append(thing.strip().split('_')[1])

            print('\n')
        print(len(measureMultipleVariants),' variants in this region.') 
        outputFile.write('\n#There are ' + str(len(measureMultipleVariants)) + ' variants in this region.')
        outputFile.write('\n#There are ' + str(len(np.unique(measureMultiplePositions))) + ' unique variant positions in this region.')
        outputFile.write('\n#There are ' + str(len(np.unique(measureMultipleVariants))) + ' unique subjects with variant in this region.')
        outputFile.write('\n#These variants are taken from the file: ' + INPUT_pVals_file + '\n')
        print('The subjects are: ')
        [print(measureMultipleVariants[g]) for g in range(len(measureMultipleVariants))]
        print(len(np.unique(measureMultipleVariants)), ' unique subjects with variant in this region: ', np.unique(measureMultipleVariants))
        lineCount = lineCount + 1

      store_rangeStart = rangeStart
      store_rangeStop = rangeStop
      store_controlCounts = controlCounts
      store_caseCounts = caseCounts

    else:
      break
  sys.exit("Completed")

