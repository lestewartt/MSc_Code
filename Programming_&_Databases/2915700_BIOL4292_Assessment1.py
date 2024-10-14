'''
Comment on University of Glasgow's AI Guidance and Policy:

Artificial intelligence should be used appropriately and responsibly, and it should not be used as a replacement for independent thinking.
All use of ChatGPT was as an information tool and not for directly generating code.
ChatGPT was used to assist in debugging two errors in logic, highlighted by **

Statement of Academic Integrity: This submission is a product of this student's work and approach to solving the assigned task. 

'''
#import libraries 
import re 
import sys
import logging 


#using sys.argv for command line configuration --> set to accept two input files
SAMfile = sys.argv[1]
TXTfile = sys.argv[2]


#create logger to help with error handling 
logger = logging.getLogger()
#selecting logging level 
logger.setLevel(logging.ERROR)
#create stream handler, define the format, add this to logger 
sh = logging.StreamHandler()
#provides information on the type of error, the date and time it occured, and message associated with the error 
sh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(sh)


"""""
Logic for the SAM file as follows: 
filter based on alignment (i.e. aligns once) 
select cigar reads with matches, deletions, and skipped regions 
add these values to the position where the alignment starts (POS) to create a boundary where the intron exists 
create a list that contains unique junction dictionaries and use this list to compile counts of supporting reads
this list will then be used to filter information from the text file 
"""

#writing a function for the block of code that parses the CIGAR string 
def cigar_info(CIGAR):
    '''
    The return function below contains the code for parsing the CIGAR string. 
    'int(match.group(1))' isolates the numeric count of the regular expression as an integer.
    'match.group(2)' isolates the character of the regular expression, which for our purposes will be selected for 'M', 'D', or 'N.'
    '''
    return [(int(match.group(1)), match.group(2)) for match in re.finditer(r'(\d+)([A-Z])', CIGAR)]

#function test: here, what would be 'CIGAR' is replaced with a test cigar string.  
#when the assert is run, the test string is passed through the function 'cigar_info' and will produce the result defined below
#if this was not the result, there would be an assertion error that interrupts the code (this is not the case here, code works)
test_cigar_info = cigar_info("2M2D2N2M")
assert test_cigar_info == [(2, 'M'), (2, 'D'), (2, 'N'), (2, 'M')]


#writing a function for the block of code that stores junction information and supporting reads 
def junction_info(junction_counts):
    found = False #initially set to false and will be used as a check to see if the current junction in the loop already exists or not 
    '''
    The subsequent loop checks all the junction information (existing_junction) that is contained within the list (junction_counts).
    If there is a match found for all variables, then a count is added to the existing junction and found is then changed to True.
    This signifies that a duplicate was identified.  
    After this loop, another conditional is introduced to check if 'found' is still False. 
    If it is, this means the junction is new and is then added to the dictionary (with a count of 1).
    '''
    for existing_junction in junction_counts: 
        if (existing_junction['RNAME'] == RNAME
            and existing_junction['start'] == start_position
            and existing_junction['end'] == end_position):
            existing_junction['count'] += 1 #adds a count to the existing dictionary for that junction 
            assert existing_junction['count'] > 1 #assertion to make sure the count was properly added (should always be greater than 1 here)
            found = True 
            continue 
    if found is False:
            junction_counts.append({'RNAME': RNAME, 'start': start_position, 'end': end_position, 'count': 1})
            
            
#junction list to store the information for each junction, including the count of supporting reads 
junction_counts = []

#added try/except catch here in the case that a non-existent file is input 
try:
    with open(SAMfile) as sam:
        for line in sam:
            #there is a variable of headers in the SAM file --> skipping them using a conditional 
            if line.startswith('@'): 
                pass
            else:
                #formatting the SAM file as a TAB delimited file 
                line = line.rstrip()
                line = line.split('\t')

                #pull the necessary columns from SAM file 
                RNAME = line[2] #chromosome 
                POS = int(line[3]) #position on chromosome where the alignment begins 
                CIGAR = line[5] #describes the alignment 
                ALIGN = line[-1] #how many times the read aligns 
                
                #created to keep track of the 'N' count within each cigar string
                n_count = CIGAR.count('N')

                #selecting for reads that align only once 
                if ALIGN != 'NH:i:1': 
                    continue 
                '''
                The below code is a continuation of parsing the CIGAR string, building from the function defined 'cigar_info.' 
                'current_position' is initially defined as POS, as this will be the baseline for counting the position of M, D, and N. 
                    For M and D in the cigar string, 'current position' is updated to be the position of the alignment plus the length of the match/deletion.
                    For N, the start positions and end positions are then defined and updated based on the length of the skipped region.
                    
                This loop also introduces the function 'junction_info,' which is described above. 
                This effectively stores the current junction information in the loop by either adding it to the count of a previously existing junction dict,
                or by adding a new junction if it does not already exist. 
                '''
                #setting a baseline for the length
                current_position = POS 
                #introducing the cigar function using a for loop, and defining the two match groups as length for match.group 1 and region for match.group 2 
                for length, region in cigar_info(CIGAR):

                    #set to move on to the next cigar string once the count of skipped regions in the current string is zero
                    #the count of matches and deletions after the skipped regions are not relevant for the purposes of this calculation 
                    if n_count == 0: 
                        continue
                    elif region == 'M' or region == 'D':
                        #updating the 'current position' as the POS and the M/D region added together
                        #useful for setting the starting boundary of the junction
                        current_position = current_position + length 
                    elif region == 'N': 
                        #defining the start position of the intron as the end of the M or D region before it 
                        #necessary to define as its own variable here outside of 'current position' so it is not overwritten
                        start_position = current_position 
                        #defining the end position of the intron using the established start position and the length of the skipped region 
                        end_position = start_position + length 
                        #the current position is then updated to the end of the intron. very important for strings with multiple skipped regions. 
                        current_position = current_position + length # **
                        #set to subtract one from 'n.count' for every skipped region found 
                        #when all the skipped regions in the cigar are counted, the count will then be zero.
                        n_count -= 1  

                        #putting the above information into the junction function to be added to the junction list, either as a count or a new junction
                        junction_info(junction_counts)

                    else: #This accounts for insertions (I) and soft-clipped regions (S) that were not considered in this calculation. 
                        #Assertion to make sure no skipped regions are missed 
                        assert region != 'N' 
                        pass 
except FileNotFoundError:
    logger.error(f"File {SAMfile} not found. Please check and try again. \n")
    raise SystemExit(1)

#final count of unique junctions is 227

""""
Logic for Tab File: 
Open text file with gene locations and store the relevant variables.
Using the junction list created in the sam file loop, match the junction locations within their gene boundaries.
Save the junction locations (and their supporting reads) for each gene in a new output file. 

"""

#creating a new text file (with the requested file name format)  
#added the same try/except condition as above, in case a file is not input correctly/does not exist
try:
    with open('2915700.txt', 'w') as junction_table:
        with open(TXTfile) as tab:
            #skipping the header line 
            header = next(tab) 
            #defining a variable that can be used inside the for loop to compare the current gene in the loop to the previous one 
            prev_gene_ID = None #Setting to 'None' outside of the loop, because in the first iteration there will be no previous gene to compare 

            for row in tab: 
                #formatting the text file 
                #replace function is used to remove the commas in the column for the genomic location 
                row = row.rstrip().replace(',','')
                row = row.split('\t')
                #defining relevant variables 
                gene_ID = row[0]
                #splitting the genomic location column into its relevant parts-- chromosome and positional information are separated by a colon 
                chromosome, gene_location = row[2].split(':')
                #using a regular expression to select the start and end positions of the genes 
                match = re.search(r'(\d+)\.\.(\d+)', gene_location)
                if match:
                    start = int(match.group(1))
                    end = int(match.group(2))
                else:
                    #Error catch if no match is found 
                    logger.error(f"Error identifying start and end positions for the gene: {gene_ID}")
                    continue

                '''
                The below conditional checks if it is the first gene in the loop  
                and if the current gene in the loop is different from the previous
                '''

                if prev_gene_ID is not None and gene_ID != prev_gene_ID: 
                    #if prev_gene_ID is not 'None' = there is a gene already iterated that the current gene should be compared to
                    #If the above conditions are met, then the current gene is different from the previous and a space between the different genes will be added  
                    junction_table.write('\n') 
        
                #For junctions associated with the gene, because there are no junctions iterated in the loop yet
                gene_written = False #used alongside 'previous_gene_ID' to keep track of the previous and current gene
                
                #pulling items from the previously created 'junction_counts' list 
                for junction in junction_counts: 
                    junction_chromosome = junction['RNAME']
                    junction_start = junction['start']
                    junction_end = junction['end']
                    count = junction['count']

                    #using a conditional to match not only the chromosomes but also ensuring the junction falls within the gene boundary 
                    if chromosome == junction_chromosome and junction_start > start and junction_end < end:
                        #if it satisfies the above conditions, it is written to the output file 
                        junction_table.write(f'{gene_ID}\t{junction_start}\t{junction_end}\t{count}\n')
                        gene_written = True #Set to 'True' once junctions are added 
                #using if gene was written to the output file (i.e. if it had junctions) because it is a unique, distinguishing factor between the genes
                if gene_written is True: 
                    #Gene_written is now true, which then updates the previous gene to be the gene that was just used in the loop. Sets up for the next iteration. 
                    prev_gene_ID = gene_ID
                else: 
                    prev_gene_ID = None
                    '''
                    **
                    The above code sets 'prev_gene_ID' back to 'None' if no junctions were added for the current gene in the loop (i.e. gene_written is still False). 
                    This accounts for instances where there are several rows of different genes but no junctions are found in those gene boundaries and therefore not added.
                    Otherwise, previous_gene_ID will continue to be stored as these genes with no junctions, and so the only things added to the output from the loop are the empty rows. 
                    Thus, this ensures proper formatting so that extra spaces are not added between rows. 
                    '''
except FileNotFoundError:
    logger.error(f"File {TXTfile} not found. Please check and try again. \n")
    raise SystemExit(1)

#output in 'junction_table' has 193 rows 

#user indicator that program has sucessfully run and output file was written
logger.info('Program complete. Junction table has been generated, designated by GUID 2915700.')        

            
       

            
            
        
        
           
