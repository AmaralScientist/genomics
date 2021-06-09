
###############################################################################################################################
## November 2019 Michelle Amaral
## To run: python Variant_dictionary.py CADDforFiltering MAFlowerBound MAFupperBound OUTPUT_directory
## python CBD_variant_dictionary.py 0 0.04 0.75 .
## The above example is for CADD greater than 0, and MAF greater than 0.04 but less than 0.75
## 
## Example: python Variant_dictionary.py 0 0.04 0.75 .
##
## first the script creates a dictionary of variants (big_variant_dict) filtered on CADD and MAF; more filtering
## criteria can be easily added or changed
##
## quality control checks are built in that ensure the same number of variants that enter are filtered out or retained.
## 
## There are 3x output files: (1) the main result output begins with prefix 'filtered_Variant_Counts_for_Chr_' and consists of 
## variant position, number of counts in the control group and number of counts in the case group. This file is input to 
## Sliding_Window_Analysis.py (2) the big variant dictionary is saved for later use by invoking pickle. The file prefix is
## 'big_variant_dict_for_Chr_' and is a .p file. This can be used to obtain information on the variants of interest. (3) a log
## file is output with prefix 'log_output_Chr_'. It is also used as input to Sliding_Window_Analysis.py to extract the total
## number of variants remaining after filtering.
###############################################################################################################################



import sys,re,os
import gzip
from itertools import islice
import pickle

#INPUT_file = sys.argv[1]
#INPUT_Chromosome = sys.argv[2]
CADD_to_filter = float(sys.argv[1])
lower_MAF_to_filter = float(sys.argv[2])
upper_MAF_to_filter = float(sys.argv[3])
OUTPUT_directory = sys.argv[4]

#######################################
## functions to wrangle vep annotations
#######################################

def multiples_pull_annotations(header, descriptors, annotationList, position, how_many_alts):
  # handles parsing and key, values of variant with multiple alternate alleles and multiple vep annotations per allele
  for line in header: 
    line_split = line.split()
    if (line_split[0].startswith('##')) and ("CSQ" in line):
      sColon = re.search(':(.*?)>',line)
      colon = sColon.group(1)

      descriptors_split = colon.strip().split("|")
      for word in descriptors_split:
        descriptors.append(word)

  # determine how many vep annotations there are
  vepSplitAnnotationOnComma = annotationList.strip().split(',')
  numVepAnnotations = len(vepSplitAnnotationOnComma)

  # pull out the first annotation for the respective alternative allele, note use of position as index
  # position is k, the index from main body that loops through each of the alternate alleles 
  firstAnnotation = vepSplitAnnotationOnComma[position].strip().split('|')
  # declare dictionary and initialize it with the values from the first annotation for the respective alternative allele
  annotation_dictionary = {k:[v] for k,v in zip(descriptors_split, firstAnnotation)}
  for key in annotation_dictionary:
    if annotation_dictionary[key] == ['']:
      annotation_dictionary[key] = ['NA']

  # since there are multiple vep annotations per alternate allele, increment the position counter
  # by the number of alternate alleles. This 'secondPosition' will be the second vep annotation added to the dictionary 
  secondPosition = position + how_many_alts

  # to add the remaining vep annotations to the dictionary, loop through the list of full vep annotations
  # (vepSplitAnnotationOnComma) and pull out those that will correspond to the respective allele.
  # to do this, begin with vepSplitAnnotationOnComma[secondPosition] and increment the index by the number of
  # alternative alleles 
  for m in list(range(secondPosition, numVepAnnotations, how_many_alts)):
    vepSplit = vepSplitAnnotationOnComma[m].strip().split('|')        
    for y in list(range(0, len(vepSplit), 1)):
      if vepSplit[y] == '':
        vepSplit[y] = 'NA'
        annotation_dictionary[descriptors[y]].append(vepSplit[y])
      else:
        annotation_dictionary[descriptors[y]].append(vepSplit[y])

  return(annotation_dictionary)

def single_alt_pull_annotations(header, descriptors, annotationList):

  for line in header: 
    line_split = line.split()
    if (line_split[0].startswith('##')) and ("CSQ" in line):
      sColon = re.search(':(.*?)>',line)
      colon = sColon.group(1)

      descriptors_split = colon.strip().split("|")
      for word in descriptors_split:
        descriptors.append(word)
  # one vep annotation contains 232 different features 
  if len(annotationList) == 232:
    annotation_dictionary = {k:[v] for k,v in zip(descriptors_split, annotationList)}
    for key in annotation_dictionary:
      if annotation_dictionary[key] == ['']:
        annotation_dictionary[key] = ['NA']

  # if there is more than one vep annotation but only one alternate allele
  else:      
    # declare dictionary and initialize with first vep annotation 
    firstAnnotation = annotationList[0].strip().split('|')
    annotation_dictionary = {k:[v] for k,v in zip(descriptors_split, firstAnnotation)}
    for key in annotation_dictionary:
      if annotation_dictionary[key] == ['']:
        annotation_dictionary[key] = ['NA']

    # split remaining vep annotations and append values to the existing keys
    for m in list(range(1, len(annotationList), 1)):
      vepSplit = annotationList[m].strip().split('|')

      for y in list(range(0, len(vepSplit), 1)):
        if vepSplit[y] == '':
          vepSplit[y] = 'NA'
          annotation_dictionary[descriptors[y]].append(vepSplit[y])
        else:
          annotation_dictionary[descriptors[y]].append(vepSplit[y])

  return(annotation_dictionary)


def check_MAF(annotation_dictionary, num_of_veps, extended_position):  #might not need num of veps here
  try:
    for b in list(range(0, len(num_of_veps), 1)):
      if (annotation_dictionary['gnomadGe_AF'][0] == '.' or annotation_dictionary['gnomadGe_AF'][0] == 'NA' or annotation_dictionary['gnomAD_AF'][0] == '.' or annotation_dictionary['gnomAD_AF'][0] == 'NA' or annotation_dictionary['TOPMed_AF'][0] == '.' or annotation_dictionary['TOPMed_AF'][0] == 'NA'):
        return(1) # filter out
      # is the gnomAD frequency less than the filtering threshold (MAF_to_filter)?
      elif ( lower_MAF_to_filter < float(annotation_dictionary['gnomadGe_AF'][0]) and float(annotation_dictionary['gnomadGe_AF'][0]) < upper_MAF_to_filter ) and ( lower_MAF_to_filter < float(annotation_dictionary['gnomAD_AF'][0]) and float(annotation_dictionary['gnomAD_AF'][0]) < upper_MAF_to_filter ):
        # if gnomAD frequency is less than threshold, keep the variant if the freq in the Bravo database is 'NA'
        if annotation_dictionary['TOPMed_AF'][0] == 'NA':
          return(1) # filter out
        # if gnomAD frequency is less than threshold and the freq in the Bravo database is not 'NA'
        # check to see if the Bravo annotation contains a '&' and then check whether the frequency in 
        # Bravo database is less than the threshold
        else:
          bravoSplit = annotation_dictionary['TOPMed_AF'][0].split('&')
          numberOfFrequencies = len(bravoSplit)
          # There will only be 1 frequency value if there is no '&'
          # if that value is less than the threshold, keep the variant
          if numberOfFrequencies == 1 and ( lower_MAF_to_filter < float(bravoSplit[0]) ) and ( float(bravoSplit[0]) < upper_MAF_to_filter ):
            return(0) # keep
          # if that value is greater than the threshold, filter out the variant
          elif ( numberOfFrequencies == 1 and (lower_MAF_to_filter > float(bravoSplit[0])) ) or ( numberOfFrequencies == 1 and (upper_MAF_to_filter < float(bravoSplit[0])) ):
            return(1) # filter out
          # if there are '&' in the Bravo annotation
          elif numberOfFrequencies > 1:
            # counter for number of Bravo annotations less than threshold
            correctFreqInterval = 0
            # counter for number of Bravo annotations greater than threshold
            outsideFreqInterval = 0
            # iterate through the number of annotations that were 
            # separated by '&'
            for i in list(range(0, numberOfFrequencies, 1)):
              if ( lower_MAF_to_filter < float(bravoSplit[i]) ) and ( float(bravoSplit[i]) < upper_MAF_to_filter ) :
                # freq is less than threshold, so increment counter
                correctFreqInterval+=1
              else:
                # freq is greater than threshold, so increment counter
                outsideFreqInterval+=1
            # if at least one of those values are less than threshold, keep variant
            if correctFreqInterval > 0:
              return(0) # keep
            # otherwise all were above the threshold, so filter out
            else: 
              return(1) # filter out
      # The gnomAD frequency is greater than the filtering threshold so filter it out
      else:
        return(1) # filter out
  except ValueError:
    print('Value error at: ', extended_position)
    print('Bravo annotation: ', annotation_dictionary['TOPMed_AF'][0])
    return (0)



###############################################
## functions to obtain subject variant counts
###############################################

def get_variant_counts(genotype_element, list_of_sample_ids, allelecount):
  number_of_variants = []
  #print(extended_position)

  if len(genotype_element) != len(list_of_sample_ids):
    print('ERROR')

  for i in list(range(0, len(genotype_element), 1)):
    genotype = genotype_element[i].split(':')[0]
    
    if genotype == './.':
      variant_count = 'NA'
      #print(variant_count)
      number_of_variants.append(variant_count)
      continue

    else:
      total_depth = genotype_element[i].split(':')[3]
      allele_depth = genotype_element[i].split(':')[1]
    
      allele_depth_Ref = allele_depth.split(',')[0]
      allele_depth_Alt = allele_depth.split(',')[1]

      if total_depth == '.':
        variant_count = 'NA' 
        number_of_variants.append(variant_count)
        continue

      if (genotype == '0/1') and (float(total_depth) >= 10) and (int(allele_depth_Alt) != 0) and (int(allele_depth_Alt) >= (int(total_depth)*0.2)):          
        variant_count = 1         
      elif (genotype == '1/1') and (float(total_depth) >= 10) and (int(allele_depth_Alt) >= 8):
        variant_count = 2    
      else:
        variant_count = 0
        
      number_of_variants.append(variant_count)
  

      #print(number_of_variants)
  
  total_variants = 0
  for i in list(range(0, len(number_of_variants), 1)):
    if number_of_variants[i] == 'NA':
      continue
    else:
      total_variants = total_variants + int(number_of_variants[i])
  print(total_variants)



  if len(list_of_sample_ids) != len(number_of_variants):
    print('ERROR at line 176')
      #annotation_dictionary = {k:[v] for k,v in zip(list_of_sample_ids, firstAnnotation)}
  variant_count_dictionary = {k:v for k,v in zip(list_of_sample_ids, number_of_variants)}
  #print(variant_count_dictionary)
  #print('Total variants: ', sum(variant_count_dictionary.values()))
  #print('Allele count: ', allelecount)
  #if allelecount == sum(variant_count_dictionary.values()):
  #  print('Yes, equal')
  return(variant_count_dictionary)

def multiple_alts_get_variant_counts(genotype_element, list_of_sample_ids, allelecount, k, extended_position):
  number_of_variants = []
  k = k+1
  #print(genotype_element)
  #print('length of genotype element: ', len(genotype_element))
  #print('length of sample ids: ', len(list_of_sample_ids))
  if len(genotype_element) != len(list_of_sample_ids):
    print('ERROR')

  for i in list(range(0, len(genotype_element), 1)):
    #print(genotype_element)
    genotype = genotype_element[i].split(':')[0]
    
    if genotype == './.':
      variant_count = 'NA'
      #print(variant_count)
      number_of_variants.append(variant_count)
      continue
    if '|' in genotype:
      variant_count = 'NA'
      number_of_variants.append(variant_count)
      continue

    total_depth = genotype_element[i].split(':')[3]
    #print(total_depth)
    #print('\n')
    all_allele_depths = genotype_element[i].split(':')[1]

    # split into separate alt allele depths
    alt_allele_depth = all_allele_depths.split(',')[k]

    split_genotype = genotype.split('/')

    first_genotype = int(split_genotype[0])
    second_genotype = int(split_genotype[1])

    if first_genotype == k and (int(total_depth) >= 10) and ( (int(alt_allele_depth) >= (int(total_depth)*0.2)) ):
      first_variant_count = 1
    else:
      first_variant_count = 0
    # then check value of second genotype
    if second_genotype == k and (float(total_depth) >= 10) and ( (int(alt_allele_depth) >= (int(total_depth)*0.2)) ):
      second_variant_count = 1
    else:
      second_variant_count = 0
    #print('Second variant count: ', second_variant_count)
    
    total_variant_count = first_variant_count + second_variant_count
    #if total_variant_count != 0:
    #  print('Total variant score for this subject: ', total_variant_count)
    number_of_variants.append(total_variant_count)
  
  #print('sum of variants for this position: ', sum(number_of_variants))

  #print('length of number of variants: ', len(number_of_variants))
  #print('length of sample ids list: ', len(list_of_sample_ids))
  variant_count_dictionary = {key:v for key,v in zip(list_of_sample_ids, number_of_variants)}
  return(variant_count_dictionary)

###########################################################################
## function to obtain case / control variant counts at given position
## insert path to subject key
###########################################################################

def get_responder_nonresponder_counts(big_variant_dict):
  subjectDict = {}
  with open(' ', "r") as subjectKey:
    for line in subjectKey:
      split_subjectKey = line.strip().split('\t')
      #print(split_subjectKey)
      if split_subjectKey[0].startswith('record_id'):
        continue
      else:
        subjID = split_subjectKey[3]
        status = split_subjectKey[37]
        if status == "0":
          status = "Responder"
        elif status == "1":
          status = "Nonresponder"
        subjectDict[subjID] = status

    return (subjectDict)


###################
## body
###################

# open txt file with genes of interest
# parse each line for chromosome and position (and gene name)
# use chromosome number to gzip.open that chromosome's file
# for line in input file

with open(' ', 'r') as coordinateFile:
  for txtLine in coordinateFile:
    big_variant_dict = {}

    didNotGetUsed = []
    gotUsed = []
    countLines = 0

    split_ea_line = txtLine.strip().split(', ')
    nameOfGene = split_ea_line[0]
    chromosomeOfGene = split_ea_line[1]
    geneStartsAt = int(split_ea_line[2]) - 50000
    geneEndsAt = int(split_ea_line[3]) + 50000

    main_results_output = 'filtered_Variant_Counts_for_Gene_' + nameOfGene + '_CADD_' +  str(CADD_to_filter) + '_MAFgreaterthan_' + str(lower_MAF_to_filter) + '_butLessThan_' + str(upper_MAF_to_filter)
    print(main_results_output)

    log_output = 'log_output_Gene_' + nameOfGene + '_CADD_' + str(CADD_to_filter) + '_MAFgreaterthan_' + str(lower_MAF_to_filter) + '_butLessThan_' + str(upper_MAF_to_filter) + '.log'
    print(log_output)


    with gzip.open(' ' , "rt") as inputFile:
      with open (log_output, 'a') as logOutput:
        head = list(islice(inputFile, 3599))

        for line in inputFile:
          line_split2 = line.strip().split()
          if line_split2[0].startswith('#'):
            samp_ids = [line_split2[i] for i in list(range(9,len(line_split2), 1))]      
          else:
            descriptors = []
            chromosome = line_split2[0]
            position = int(line_split2[1])
            ref = line_split2[3]
            alts = line_split2[4]
            separate_alts = alts.split(',')
            how_many_alts = len(separate_alts)
            quality = line_split2[5]
            vqsrt = line_split2[6]
            if not vqsrt == "PASS":
              continue

            info_column = line_split2[7]
            split_genotypes = line_split2[9:len(line_split2)]
            
            # check to see if there are multiple alternate alleles
            ###################
            ## one alt allele
            ###################
            if how_many_alts == 1 and ( position >= geneStartsAt and position <= geneEndsAt ):        
              countLines = countLines + 1    
              extended_position = chromosome + '_' + str(position) + '_' + ref + '_' + alts
              
              sallelecount = re.search('AC=(.*?);', info_column)
              allelecount = sallelecount.group(1) 

              # filter on batch allele count, allow up to x counts in the batch
              if int(allelecount) < 200:
                #print('batch AC is: ', allelecount)
                pass
              else: 
                didNotGetUsed.append(extended_position)
                continue

              # pull out CADD score and filter on it
              sCADD = re.search(';CADD_scaled=(.*?);', info_column)
              CADD = sCADD.group(1)
              if CADD == 'NA' or CADD == '':
                didNotGetUsed.append(extended_position)
                continue
              #print(CADD)
              if float(CADD) > CADD_to_filter:
                #gotUsed.append(extended_position)
                pass
              else:
                didNotGetUsed.append(extended_position)
                continue
              sVEP = re.search('CSQ=(.*?\Z)', info_column)
              VEP = sVEP.group(1)
              # multiple vep annotations are split with ','
              vepSplitOnComma = VEP.strip().split(',')
              if len(vepSplitOnComma) == 1:  # there is one vep annotation for the single alternate allele
                vepSplit = VEP.strip().split('|')
              
              else:  # there is more than one vep annotation for the single alternate allele
                vepSplit = vepSplitOnComma
              annotation_dictionary = single_alt_pull_annotations(head, descriptors, vepSplit)
              # filter on MAF
              if not check_MAF(annotation_dictionary, vepSplitOnComma, extended_position) == 0:
                didNotGetUsed.append(extended_position)
                continue
              else:
                gotUsed.append(extended_position)
                big_variant_dict[extended_position] = annotation_dictionary

                # if variant made it through the filters, then add AC to its dictionary entry
                big_variant_dict[extended_position]["batch_AC"] = allelecount

                # if variant made it through the filters, then add CADD to its dictionary entry
                big_variant_dict[extended_position]["CADD"] = CADD

                variant_count_dictionary = get_variant_counts(split_genotypes, samp_ids, allelecount)
                big_variant_dict[extended_position]["subject_ids"] = variant_count_dictionary

            
            ###################
            ## >1 alt alleles
            ###################
            if how_many_alts > 1 and ( position >= geneStartsAt and position <= geneEndsAt ):
              sVEP = re.search('CSQ=(.*?\Z)', info_column)
              VEP = sVEP.group(1)

              ###############################################################################################
              ## * alleles
              ## in the loop below, n should be equal to the position within the vep annotation list
              ## w should be equal to the correct position in any list outside the vep annotation list 
              ###############################################################################################

              if '*' in separate_alts:
                # iterate through these multiple alternate alleles
                w = 0
                vepSplitOnComma = VEP.strip().split(',')

                while separate_alts[w] != '*':
                  countLines = countLines + 1
                  extended_position = chromosome + '_' + str(position) + '_' + ref + '_' + separate_alts[w]

                  sallelecount = re.search('AC=(.*?);', info_column)
                  allelecount = sallelecount.group(1)
                  allelecountSplit = allelecount.split(',')

                  # filter on batch allele count
                  if int(allelecountSplit[w]) < 200:
                    #print('batch AC is: ', allelecountSplit[w])
                    pass
                  else: 
                    didNotGetUsed.append(extended_position)
                    w+=1
                    continue

                  sCADD = re.search(';CADD_scaled=(.*?);', info_column)
                  CADD = sCADD.group(1)
                  CADDsplit = CADD.split(',')
                  
                  # filter on CADD
                  if CADDsplit[w] == 'NA' or CADDsplit[w] == '' or float(CADDsplit[w]) < CADD_to_filter:
                    didNotGetUsed.append(extended_position)
                    w+=1
                    continue
                  else:
                    vepSplit = vepSplitOnComma[w].strip().split(',')
                    annotation_dictionary = single_alt_pull_annotations(head, descriptors, vepSplit)
                    # filter on MAF
                    if not check_MAF(annotation_dictionary, vepSplitOnComma, extended_position) == 0:
                      didNotGetUsed.append(extended_position)
                      w+=1
                      continue
                    else:
                      gotUsed.append(extended_position)
                      big_variant_dict[extended_position] = annotation_dictionary
                  
                      big_variant_dict[extended_position]["batch_AC"] = allelecountSplit[w]
                      big_variant_dict[extended_position]["CADD"] = CADDsplit[w]

                      #variant_count_dictionary = get_variant_counts(split_genotypes, samp_ids, allelecount)
                      variant_count_dictionary = multiple_alts_get_variant_counts(split_genotypes, samp_ids, allelecountSplit[w], w, extended_position)
                      big_variant_dict[extended_position]["subject_ids"] = variant_count_dictionary
                      w+=1         
                
                # the following will be executed if separate_alts[w] == '*'
                if separate_alts[w] == '*':
                  w+=1
                  if w == how_many_alts:
                    continue
                  else:
                    pass

                else:
                  print('ERROR at line 543')

                # if there are more alleles after the '*' allele, they will enter this while loop
                while w != how_many_alts:
                  countLines = countLines + 1
                  extended_position = chromosome + '_' + str(position) + '_' + ref + '_' + separate_alts[w]
                  sallelecount = re.search('AC=(.*?);', info_column)
                  allelecount = sallelecount.group(1)
                  allelecountSplit = allelecount.split(',')

                  # filter on batch allele count, allow up to 3x counts in the batch
                  if int(allelecountSplit[w]) < 200:
                    #print('batch AC is: ', allelecountSplit[w])
                    pass
                  else: 
                    didNotGetUsed.append(extended_position)
                    w+=1
                    if w == how_many_alts:
                      break
                    else:
                      continue
                
                  sCADD = re.search(';CADD_scaled=(.*?);', info_column)
                  CADD = sCADD.group(1)
                  CADDsplit = CADD.split(',')
                ######################
                # below might be a problem in the future if I'm not filtering on CADD
                # will need to refigure code to line up vep annotations (or just leave in
                # the if CADDsplit[w] == 'NA' part)
                ######################
                  if CADDsplit[w] == 'NA' or CADDsplit[w] == '' or float(CADDsplit[w]) < CADD_to_filter:
                    didNotGetUsed.append(extended_position)
                    w+=1
                    if w == how_many_alts:
                      break
                    else:
                        continue

                  else:
                    n = w-1  
                    vepSplit = vepSplitOnComma[n].strip().split(',')
                    annotation_dictionary = single_alt_pull_annotations(head, descriptors, vepSplit)
                      # filter on MAF
                    if not check_MAF(annotation_dictionary, vepSplitOnComma, extended_position) == 0:
                      didNotGetUsed.append(extended_position)
                      w+=1
                      continue
                    else:
                      if extended_position == '22_18737130_T_G':
                        print('about to get added')
                      gotUsed.append(extended_position)
                      big_variant_dict[extended_position] = annotation_dictionary
                  
                      big_variant_dict[extended_position]["batch_AC"] = allelecountSplit[w]
                      big_variant_dict[extended_position]["CADD"] = CADDsplit[w]

                      #variant_count_dictionary = get_variant_counts(split_genotypes, samp_ids, allelecount)
                      variant_count_dictionary = multiple_alts_get_variant_counts(split_genotypes, samp_ids, allelecountSplit[w], w, extended_position)
                      big_variant_dict[extended_position]["subject_ids"] = variant_count_dictionary
                      w+=1

                ###################
                ## END * alleles
                ###################
                
              else: # there are no * alleles
                # loop through these multiple alternate alleles
                for k in list(range(0, how_many_alts, 1)):
                  extended_position = chromosome + '_' + str(position) + '_' + ref + '_' + separate_alts[k]
                  countLines = countLines + 1

                  vepSplitOnComma = VEP.strip().split(',')

                # check whether there is one vep annotation for each alt allele
                  if len(vepSplitOnComma) == how_many_alts:   # there is one vep annotation for each alt allele
                    
                    # pull out AC, can add some code for filtering
                    sallelecount = re.search('AC=(.*?);', info_column)
                    allelecount = sallelecount.group(1)
                    allelecountSplit = allelecount.strip().split(',')

                    # filter on batch allele count
                    if int(allelecountSplit[k]) < 200:
                      #print('batch AC is: ', allelecountSplit[k])
                      pass
                    else: 
                      didNotGetUsed.append(extended_position)
                      continue

                    # pull out CADD score and filter on it
                    sCADD = re.search(';CADD_scaled=(.*?);', info_column)
                    CADD = sCADD.group(1)
                    CADDsplit = CADD.strip().split(',')
                    if CADDsplit[k] == 'NA' or CADDsplit[k] == '':
                      didNotGetUsed.append(extended_position)
                      continue
                    if float(CADDsplit[k]) > CADD_to_filter:
                      #gotUsed.append(extended_position)
                      pass
                    else:
                      didNotGetUsed.append(extended_position)
                      continue

                    vepSplit = VEP.strip().split(',')[k].split(',')
                    annotation_dictionary = single_alt_pull_annotations(head, descriptors, vepSplit)
                    # filter on MAF
                    if not check_MAF(annotation_dictionary, vepSplitOnComma, extended_position) == 0:
                      didNotGetUsed.append(extended_position)
                      continue
                    else:
                      gotUsed.append(extended_position)
                      big_variant_dict[extended_position] = annotation_dictionary
            
                      # variant passed filtering, so add AC to the dictionary
                      big_variant_dict[extended_position]["batch_AC"] = allelecountSplit[k] 

                      # variant passed filtering, so add CADD to the dictionary
                      big_variant_dict[extended_position]["CADD"] = CADDsplit[k]

                      variant_count_dictionary = multiple_alts_get_variant_counts(split_genotypes, samp_ids, allelecount, k, extended_position)

                      big_variant_dict[extended_position]["subject_ids"] = variant_count_dictionary
                  
                  else: # there is more than one vep annotation for each alt allele              
                    # pull out AC, can add some code to filter on it
                    sallelecount = re.search('AC=(.*?);', info_column)
                    allelecount = sallelecount.group(1)
                    allelecountSplit = allelecount.strip().split(',')

                    # filter on batch allele count
                    if int(allelecountSplit[k]) < 200:
                      #print('batch AC is: ', allelecountSplit[k])
                      pass
                    else: 
                      didNotGetUsed.append(extended_position)
                      continue

                    # pull out CADD score and filter on it
                    sCADD = re.search(';CADD_scaled=(.*?);', info_column)
                    CADD = sCADD.group(1)
                    CADDsplit = CADD.strip().split(',')
                    if CADDsplit[k] == 'NA' or CADDsplit[k] == '':
                      didNotGetUsed.append(extended_position)
                      continue
                    if float(CADDsplit[k]) > CADD_to_filter:
                      pass
                    else:
                      didNotGetUsed.append(extended_position)
                      continue

                    annotation_dictionary = multiples_pull_annotations(head, descriptors, VEP, k, how_many_alts)
                    # filter on MAF
                    if not check_MAF(annotation_dictionary, vepSplitOnComma, extended_position) == 0:
                      didNotGetUsed.append(extended_position)
                      continue
                    else:
                      gotUsed.append(extended_position)
                      big_variant_dict[extended_position] = annotation_dictionary
            
                      # variant passed filtering, so add AC to the dictionary
                      big_variant_dict[extended_position]["batch_AC"] = allelecountSplit[k] 

                      # variant passed filtering, so add CADD to the dictionary
                      big_variant_dict[extended_position]["CADD"] = CADDsplit[k]

                      variant_count_dictionary = multiple_alts_get_variant_counts(split_genotypes, samp_ids, allelecount, k, extended_position)
                      big_variant_dict[extended_position]["subject_ids"] = variant_count_dictionary

    ###########################################################################
    ## big_variant_dict now generated, save it for future use
    ###########################################################################
    print('about to pickle dump')
    toSaveDict = 'big_variant_dict_for_' + nameOfGene + '_CADD_' + str(CADD_to_filter) + '_MAFgreaterthan_' + str(lower_MAF_to_filter) + '_butLessThan_' + str(upper_MAF_to_filter) + '.p'
    pickle.dump(big_variant_dict, open(toSaveDict, "wb"))

    

    ###########################################################################
    ## big_variant_dict now generated, next count variants in subjects
    ###########################################################################

    subjectDict = get_responder_nonresponder_counts(big_variant_dict)

    badIDs = []
    allVariantsDict = {}

    for position in big_variant_dict:
      responder_counts = []
      nonresponder_counts = []
      positionCount = {}
      for subject in big_variant_dict[position]['subject_ids']:
        #print(big_variant_dict[position]['subject_ids'])
        if subject in subjectDict:
          if subjectDict[subject] == 'Nonresponder':
            #print('Control: ', subject)
            if big_variant_dict[position]['subject_ids'][subject] == 'NA':
              continue
            else:
              nonresponder_counts.append(big_variant_dict[position]['subject_ids'][subject])
            #print(big_variant_dict[position]['subject_ids'][subject])
          elif subjectDict[subject] == 'Responder':
            if big_variant_dict[position]['subject_ids'][subject] == 'NA':
              continue
            else:
              responder_counts.append(big_variant_dict[position]['subject_ids'][subject])
        else:
          if subject not in badIDs:
            badIDs.append(subject)
          else:
            continue
      positionCount['Nonresponder'] = sum(nonresponder_counts)
      positionCount['Responder'] = sum(responder_counts)

      allVariantsDict[position] = positionCount

    ###########################################################################
    ## Quality Control for filtered variant dictionary
    ###########################################################################

    with open (log_output, 'a') as logOutput:
      for item in didNotGetUsed:
        if item in big_variant_dict:
          logOutput.write(item + 'yes, I am here: ' + big_variant_dict[item] + '\n')

      for item in gotUsed:
        if not item in big_variant_dict:
          logOutput.write('problem' + '\n')

      logOutput.write('Did not get used: ' + str(len(didNotGetUsed)) + '\n')
      logOutput.write('Got used: ' + str(len(gotUsed)) + '\n')
      logOutput.write ('line count: ' + str(countLines) + '\n')
      totalUsedAndNotUsed = len(didNotGetUsed) + len(gotUsed)
      if totalUsedAndNotUsed == countLines:
        logOutput.write('Used and not used equal line count' + '\n')
      else:
        logOutput.write('Used and not used ARE NOT EQUAL TO line count' + '\n')
      logOutput.write('These are the bad IDs: ' + '\n')
      for ids in badIDs:
        logOutput.write('%s\n' % ids)


    ###########################################################################
    ## End QC
    ###########################################################################


    ###########################################################################
    ## Output
    ###########################################################################

    variantCounter = 0
    with open (main_results_output, 'x') as mainOutput:
      mainOutput.write('POSITION' + '\t' + 'NONRESPONDER' + '\t' + 'RESPONDER' + '\n')
      for position in allVariantsDict:
        mainOutput.write(position + '\t' + str(allVariantsDict[position]['Nonresponder']) + '\t' + str(allVariantsDict[position]['Responder']) + '\n')
        variantCounter = variantCounter + 1

      #mainOutput.write('#After filtering, there are ' + str(variantCounter) + ' variants.')
      #print('After filtering, there are ', variantCounter, ' variants.' + '\n')

    with open (log_output, 'a') as logOutput:  
      logOutput.write('#After filtering, there are ' + str(variantCounter) + ' variants.\n')
      logOutput.write('\n\n\n\n\n\n\n')






