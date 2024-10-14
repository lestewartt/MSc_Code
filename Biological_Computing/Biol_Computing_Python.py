'''
Comment on University of Glasgow's AI Guidance and Policy:

Artificial intelligence should be used appropriately and responsibly, and it should not be used as a replacement for independent thinking.
All use of ChatGPT was as an information tool and not for directly generating code.
ChatGPT was used to assist in debugging, highlighted by **

Statement of Academic Integrity: This submission is a product of this student's work and approach to solving the assigned task. 

'''

#importing the necessary libraries 
import argparse
import logging 
import vcf
import gffutils
from Bio.Seq import Seq 
import os 
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import math 

#argparse for command line input 
parse = argparse.ArgumentParser(description='Please provide file names for analysis')
parse.add_argument('--vcf', required=True, help="Specify the VCF input file containing variant information")
parse.add_argument('--gff', required=True, help="Specify GFF input file for genome annotation information")
parse.add_argument('--fasta', required=True, help="Specify genome FASTA input file for sequence information")

args = parse.parse_args()
#defining file names to be called upon in later code
vcf_file = args.vcf 
gff_file = args.gff
fasta_file = args.fasta 

#setting up a logger to be used in error handling/providing helpful info for user guidance 
logger = logging.getLogger()
#logging level
logger.setLevel(logging.INFO)
#create stream handler and defining the format
sh = logging.StreamHandler()
#type of error, the date and time, and the associated message 
sh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(sh)

#logging info for the name of the files provided at the command line
logging.info(f'VCF File: {vcf_file}')
logging.info(f'GFF File: {gff_file}')
logging.info(f'FASTA File: {fasta_file}')


'''Marking the start of the analysis'''


#opening the vcf file with added error catching 
try:    
    read_vcf = vcf.Reader(filename = vcf_file)
except FileNotFoundError:
    logger.error(f"Error reading file. {vcf_file} is not found.")
    raise SystemExit(1)

#making a count for the variants with a quality of 20 or less 
quality_count_under20 = 0 
#creating a list that contains dicts with the corresponding information for variants with a quality over 20
quality_over20_info = []
#for loop to parse the vcf file 
for line in read_vcf:
    #pulling the relevant columns from the vcf file 
    CHROM = line.CHROM
    POS = line.POS
    REF = line.REF
    ALT = line.ALT
    QUAL = line.QUAL
    #using a conditional to delineate the variants into two "groups" based on quality scores 
    if QUAL <= 20: 
        #storing the incremental counts in the variable outside the for loop 
        quality_count_under20 += 1 
    elif QUAL > 20: 
        #storing the information for each variant that is relevant to later code and final output 
        quality_over20_info.append({'CHROM': CHROM, 'POS': POS, 'Ref': REF, 'Alt': ALT, 'QUAL': QUAL})
#creating a logging statement to deliver the results of the counts with a quality score below the threshold 
logging.info(f'Count of variants with quality score <= 20: {quality_count_under20}')



#reading the gff file -- first step is to create a database 
#naming the database file based on the gff input file (for better readability/ease of transfer)
db_gff = gff_file.replace('.gff', '.db')
#checking if the database already exists
if not os.path.exists(db_gff):
    #if the database is not already found, it will be created 
    logger.info(f'Creating Database {db_gff}')
    try: 
        #creating a new database: inputting the name of the gff file and keep_order=True will keep the order of attributes the same  
        db_gff = gffutils.create_db(gff_file, dbfn=db_gff, keep_order=True, sort_attribute_values= True)
    #general error catcher  
    except Exception as error: 
        logger.error(f'Unable to Create Database due to {error}')
        raise SystemExit(1) 
#if the database does already exist, this delivers a statement that a connection is being established to that database 
else: 
    logger.warning('Existing Database Found. Connecting.')
    #this creates a connection to the preexisting database, again maintaining the order of attributes 
    db_gff = gffutils.FeatureDB(db_gff, keep_order=True)


'''Several functions to improve modularization of code'''


#function for calculating the position of the variant in the CDS in which it is found
def calc_position_snp(SNP_position, child, feature):
    '''
    This conditional is taking advantage of the pre-established variables in the main code. 
    If the child (i.e. current CDS in the loop) is not the same as the feature (i.e. the CDS in which a variant is located),
    then find the total length of the CDS and add it to the count. 
    '''
    if child != feature: 
        cds_len = child.end - child.start + 1 
        SNP_position += cds_len
    return SNP_position

#function for connecting the CDS sequences 
def concatenate_cds(sequence, child, fasta_file):
    '''
    Takes in three parameters: sequence is the variable that stores the final sequence, child is the CDS in the loop,
    and the fasta file stores the sequence. 
    This statement below appends the current sequence to the end of the previous sequence.
    '''
    sequence = sequence + child.sequence(fasta_file, use_strand=True)
    return sequence 

#Function for creating the output file with relevant variant data 
def create_output_variant_table(variant_info_table, variants, base, variant_type, parent, protein_position, REF_aa, ALT_aa):
    
    '''
    Creating a dict within a list so that it can be easily referenced outside of the loop for the purposes of a summary output file.
    Chrom, Pos, and Ref, all originate from the dict created using the vcf file (Alt comes directly from the vcf file). 
    Type, Transcript, Protein Loc, Ref AA, and Alt AA are all derived from the following code. 
    '''
    variant_info_table.append({
        'Chrom': variants['CHROM'], 
        'Pos': variants['POS'],
        'Ref': variants['Ref'],
        'Alt': base,
        'Type': variant_type, 
        'Transcript': parent.id if variant_type == 'Synonymous' or variant_type == 'Non-synonymous' else 'NA',
        'Protein Location': protein_position if variant_type == 'Synonymous' or variant_type == 'Non-synonymous' else 'NA',
        'Ref AA': str(REF_aa) if variant_type == 'Synonymous' or variant_type == 'Non-synonymous' else 'NA',
        'Alt AA': str(ALT_aa) if variant_type == 'Synonymous' or variant_type == 'Non-synonymous' else 'NA'
    })


'''MAIN CODE'''


#creating an empty list to eventually store the dictionary information for each variant 
#each dict will be a row in the output tsv file 
variant_info_table = []
#creating a count for the types of mutations found
#this will be later used for the bar plot at the end of the script
synonymous_count = 0 
non_synonymous_count = 0 
non_coding_count = 0 
#now only filtering through variants over the quality threshold 
for variants in quality_over20_info: 
    #extracting the relevant information from the previously established list 
    CHROM, POS, REF, ALT = variants['CHROM'], variants['POS'], variants['Ref'], variants['Alt']
    #creating an empty list to append a tuple containing the variant, its calculated position, and the corresponding sequence
    #this will ultimately be important for determining the type of variant 
    variant_info = []
    #conditional below checks if the length of the list of features (of type=CDS) containing a variant is zero 
    #start and end are the same because the variant is a SNP (using the genomic position of the variant)
    if len(list(db_gff.region(seqid=CHROM, start=POS, end=POS, featuretype='CDS'))) == 0:
        #if there is no CDS, a count is added to the number of non coding features 
        non_coding_count += 1 
    else: 
        #if there are coding regions found, then the features are passed through a for loop for further processing 
        for feature in db_gff.region(seqid=CHROM, start=POS, end=POS, featuretype='CDS'):
            #print(f'{feature}\n')
            #print(f'Processing feature:{feature}\n')
            '''
            In order reconstruct the coding sequence, the other CDS's must be found within that transcript.
            GFF files use a hierarchial architecture, so we must first pull the scope of analysis from the individual CDS to its transcript.
            '''
            for parent in db_gff.parents(feature, featuretype='mRNA', order_by='start'):
                #iterating through  each CDS and adding them up along the way -- requires defining outside the loop 
                #print(f'Processing parent:{parent}\n')
                SNP_position = 0
                sequence = ''
                #positional analysis varies depending on the strand 
                if parent.strand == '+':
                #order_by start is important to ensure the features are in the correct order (i.e. not lexical order )
                    for child in db_gff.children(parent, featuretype='CDS', order_by='start'): 
                        '''
                        Obtaining the sequence of the CDS in the loop and passing it up to the function where it will be added to a growing 'overall' sequence.
                        Important note: this is not the full coding sequence. This sequence ends with the CDS in which the variant was located. 
                        This is important for determining the position of the variant within the coding sequence.
                        '''
                        sequence = concatenate_cds(sequence, child, fasta_file)

                        #Checks if the current CDS in the loop is the same CDS in which the variant is located 
                        if child == feature: 
                            #Find the distance between the start of the CDS and the position of the variant within the CDS
                            #Cuts off this CDS at the position of that variant, such that the final count is the location of the variant (relative to the actual coding sequence) 
                            #the +1 accounts for the indexing differences between python and the files
                            SNP_position += (int(POS) - child.start + 1)
                            #This loop is now considered complete, and the relevant variant information is appended as a tuple to a list outside the loop 
                            variant_info.append((variants, SNP_position, sequence))
                            #stops the current iteration 
                            break
                        #If the current CDS is not the CDS containing the variant then...
                        else: 
                            #calculating the total length of the CDS (see function)
                            SNP_position = calc_position_snp(SNP_position, child, feature) 
                
                #similar process for the negative parent strand 
                elif parent.strand == '-':
                    #same syntax as the positive strand, except it is reversed 
                    for child in db_gff.children(parent, featuretype='CDS', order_by='start', reverse=True):
                        
                        #this will provide the sequence of all the CDS's up to the CDS including the variant 
                        #(i.e. not the full coding sequence)
                        sequence = concatenate_cds(sequence, child, fasta_file) 
                        '''
                        Find the distance between the CDS and the position of the variant within the CDS
                        Cuts off this CDS at the position of that variant, such that the final count 
                        is the location of the variant (relative to the actual coding sequence) 
                        '''
                        if child == feature: 
                            #NB: because the negative strans is reversed, the calculation must also be reversed
                            #therefore, the position of the variant is treated as the 'start'
                            SNP_position += (child.end - int(POS) + 1)
                            #same convention as the prior loop
                            #once the CDS of the variant has been found, the list is updated and the next iteration begins
                            variant_info.append((variants, SNP_position, sequence))
                            break
                        else: 
                            #calculating the length of each CDS
                            #dont need to flip orientation of the maths because it has already been reversed in the above code  
                            SNP_position = calc_position_snp(SNP_position, child, feature)

        #pulling down the positions recorded above to be used to swap the bases at those designated positions in the sequence. 
        #Using a loop to evaluate one at a time.
            for info in variant_info:
                #indexing the tuple to isolate the variables 
                variant = info[0]
                position = info[1]
                seq = info[2]
                #using a 'sandwich' method, or slicing, to swap the reference for the alternative base
                #treating positions like an index in the sequence  
                for base in ALT:    
                    #first index includes all of the sequence leading up to the variant position
                    #str(base) is the specific base we want to insert
                    #the last index completes the rest of the sequence after the position of interest
                    variant_seq = seq[:position - 1] + str(base) + seq[position:]
                    #this command rounds up to give the coordinate of the amino acid in which the variant can be found 
                    protein_position = math.ceil(position / 3)

                    #defining the original and alt sequence 
                    REF_sequence = Seq(seq)
                    ALT_sequence = Seq(variant_seq)
    
                    try:
                        #translating the original sequence 
                        REF_protein_sequence = REF_sequence.translate()
                    except Exception as error:
                        logger.error(f'Translation error: {error}')
                        logger.error(f'Sequence related to error: {REF_protein_sequence}')
                        continue 
                    try:
                        #translating the mutated sequence 
                        ALT_protein_sequence = ALT_sequence.translate()
                    except Exception as error:
                        logger.error(f'Translation error: {error}')
                        continue 

                    #pulling the original amino acid that contains the variant position **
                    #genomic positions are 1 indexed and python is 0 indexed -> subtract 1 to account
                    try: 
                        REF_aa = REF_protein_sequence[protein_position - 1]
                    except IndexError as error:
                        print(f'Error in determining the reference amino acid:{error}')
                        #print(len(seq) #print(len(variant_seq))
                        #logger.error(f'Sequence related to error: {REF_protein_sequence}')
                        #print(variant) #print(f'Feature: {feature}')
                        continue
                    #pulling the amino acid that contains the variant 
                    try:
                        ALT_aa = ALT_protein_sequence[protein_position - 1]
                    except IndexError as error:
                        print(f'Error in Determining the Alt amino acid:{error}')
                        continue 

                    #comparing the amino acid in the original and alternate sequences to determine the type of mutation 
                    if REF_aa == ALT_aa:
                        #synonymous mutations do not change the amino acid 
                        #adds to the count of synonymous mutations
                        synonymous_count += 1 
                        #stores the type of mutation in a variable 
                        variant_type = 'Synonymous'
                    else: 
                        #non-synonymous mutations result in an amino acid change
                        #adds to the count of non-synonymous mutations 
                        non_synonymous_count += 1 
                        #stores the type of mutation in a variable 
                        variant_type = 'Non-synonymous'

                    #compiling all the information derived from this loop and feeding it into the function that will help create the output tsv file 
                    create_output_variant_table(variant_info_table, variants, base, variant_type, parent, protein_position, REF_aa, ALT_aa)


#write a plot for the proportion of the variants (QUAL>20) that are non-coding, synonymous, or non-synonymous 
#find the total count for all categories in order to find the proportion 
categories_total = non_coding_count + synonymous_count + non_synonymous_count
#use the above value to find the proportion for each of the criteria 
non_coding_prop = non_coding_count / categories_total
synonymous_prop = synonymous_count / categories_total
non_syn_prop = non_synonymous_count / categories_total

#organizing the x and y axes into separate lists 
categories = ['Non-Coding', 'Synonymous', 'Non-Synonymous']
proportions = [non_coding_prop, synonymous_prop, non_syn_prop]
#taking the above lists and putting them into a dataframe --> preparation for plotting 
plot_data = {"Categories": categories, 'Proportions': proportions}
data_frame = pd.DataFrame(plot_data)
#print(data_frame)
#command to make the boxplot using the dataframe created above 
sns.barplot(x='Categories', y='Proportions', data=data_frame).set(title='Proportion of Variant Type with Quality Score > 20')
#saves the plot as a PNG file for viewing 
plt.savefig('variant_type_barplot_2915700.png')

#converting the list of dictionaries i made in the 'create_output_variant_table' function into a table using pandas
variant_info_df = pd.DataFrame(variant_info_table)
variant_info_df.to_csv('2915700_variant_info.tsv', sep='\t', index=False)

#setting up a filepath that will be used in the logging statement at the end
#this just provides the path to the current working directory 
output_loc = os.getcwd()
#creating path components -- gives the path to each file 
barplot_loc = os.path.join(output_loc, 'variant_type_barplot_2915700.png')
tsv_loc = os.path.join(output_loc, '2915700_variant_info.tsv')
#this essentially combines the above information -- provides paths to both output files, as well as the current directory
#better than hardcoding the names  
logging.info(f'Script completed. Output files: {barplot_loc} and {tsv_loc}". Located in: {output_loc}')

