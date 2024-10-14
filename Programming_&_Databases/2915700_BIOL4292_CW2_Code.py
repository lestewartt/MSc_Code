'''
Message regarding University of Glasgow's AI Guidance and Policy:

Artificial intelligence should be used appropriately and responsibly, and it should not be used as a replacement for independent thinking.
All use of ChatGPT was as an information tool and not for directly generating code.
ChatGPT was used to assist in debugging, highlighted by **

Statement of Academic Integrity: This submission is a product of this student's work and approach to solving the assigned task. 

'''
#importing the relevant libraries 
import argparse 
import re
import logging
import sqlite3 
import seaborn as sns 
import pandas as pd 
import matplotlib.pyplot as plt

#argparse for command line input 
parse = argparse.ArgumentParser(description='2915700 CW2 Database Infrastructure')
#the action 'store_true' allows for these commands to be treated as optional, so they do not have to be called all at once
parse.add_argument('--createdb', action='store_true', help="Initialize Database Creation")
parse.add_argument('--loaddb', action='store_true', help="Specify relevant Files for Input into Database")
#nargs performs a similar purpose to store_true 
parse.add_argument('--querydb', type=int, nargs='?', help="Select a Query (n= num 1 through 9)")
parse.add_argument('database_file', nargs='?',help='SQLite Database File Name') 

args = parse.parse_args()

#logger for error handling 
logger = logging.getLogger()
#setting the logging level
logger.setLevel(logging.ERROR)
#create stream handler
sh = logging.StreamHandler()
#this tells the system to provide the type of error, the date and time, and the associated message when the command is initiated
sh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(sh)
    

'''
Created a class for Database Connection/Creation
This is a general class that can be recycled and adapted into other code for the purpose of 
locating the database, connecting to it, loading the data into tables, and executing database queries. 
''' 
class EstablishDatabase: 

    #constructor for the class 
    def __init__(self, database):

        #The constructor takes the database parameters and assigns it to the attribute below
        self.database = database 
        try:
            #establish a connection with the database 
            self.connection = sqlite3.connect(self.database)
            #creating a cursor to interact with the database
            self.cur = self.connection.cursor()
        except sqlite3.OperationalError as error: 
            logger.error(f'Could not establish connection to database: {error}')

    def create_tables(self, table_creation):
        #loop to iterate through each 'create table' statement in the create tables function 
        try:
            for creation in table_creation:
                self.cur.execute(creation)
        except sqlite3.DatabaseError as error: 
            logger.error(f'Error during table creation: {error}')

    def insert_db(self, sql_query, parameters): 
        try:
            #for inserting data into the tables 
            self.cur.executemany(sql_query, parameters)
        except sqlite3.IntegrityError as error:
            print(f"SQLite Integrity Error: {error}")
        self.connection.commit() 

    #The query_db function carries out the query commands 
    def query_db(self, sql_query, parameters): 
        try: 
            #using the provided parameters, the sql query is executed
            columns = self.cur.execute(sql_query, parameters).fetchall()
            logger.info("Query Successfully Completed")
        except sqlite3.ProgrammingError as error:
            logger.error(f"Query Error: {error}")
        return columns 


#simple function to store the tables being created 
'''
The organization of the database has been adapted from the theorized, original ERD. Here, visitID is stored
in the Sample table. 
One subject can have many samples, there can be many peak data associated with the samples, and there is a
many to many relationship between peak and metabolite, which necessitates a linkage table (peakannotation). 
'''
def get_tables(): 
    return [
        "CREATE TABLE IF NOT EXISTS Subject(SubjectID PRIMARY KEY, Sex INT, Age INT, BMI DECIMAL(2,2), InsulinStatus VARCHAR(2))",
        "CREATE TABLE IF NOT EXISTS Sample(SampleID, VisitID, TranscriptAbundance DECIMAL(2,8), SubjectID, PRIMARY KEY(SampleID), FOREIGN KEY(SubjectID) REFERENCES Subject(SubjectID))",
        "CREATE TABLE IF NOT EXISTS Peak(PeakNumber, SampleID, PRIMARY KEY(PeakNumber), FOREIGN KEY(SampleID) REFERENCES Sample(SampleID))",
        "CREATE TABLE IF NOT EXISTS Metabolite(MetaboliteID, PeakID, KEGG, Pathways, PRIMARY KEY(MetaboliteID))",
        "CREATE TABLE IF NOT EXISTS PeakAnnotation(PeakNumber, MetaboliteID, PRIMARY KEY(PeakNumber, MetaboliteID), FOREIGN KEY(PeakNumber) REFERENCES Peak(PeakNumber), FOREIGN KEY(MetaboliteID) REFERENCES Metabolite(MetaboliteID))"
    ]

#this function was created to process and print the results of the query in the designated format 
def format_results(result):
    #this loop iterates over each row of the answer obtained from the database 
    for answer in result:
        #this filters through each value within the row. Here, a tab separation is added between the ouput.
        #improves readability 
        for values in answer:  
            print(f'{values}', end='\t')
        #after iterating through that row, a new line will be printed for the next row 
        print()

'''
The next three functions all deal with filtering, cleaning, and loading the appropriate data into the tables 
that were established in the get_tables function. 
When argparse is called to 'loaddb', each file that has data necessary for the queries is opened, parsed, 
(cleaned/manipulated if necessary), and then inserted into the relevant table. 
'''
#opening the file that will populate the subject table 
def subject_data(subject_file): 
    try:
        with open(subject_file) as subject:
            #skip header
            next(subject)
            #setting up an empty list for the tuples containing each row of data to be appended to 
            parameter_list = [] 
            #formatting the file - only pulling columns that provide the necessary info to answer queries 
            for row in subject:
                row = row.rstrip().rsplit(',')
                SubjectID = row[0]
                Sex = row[2]
                Age = row[3]
                BMI = row[4]
                InsulinStatus = row[6]
                #setting up the structure of the parameters that will be stored in the table 
                parameters = (SubjectID, Sex, Age, BMI, InsulinStatus)
                parameter_list.append(parameters)
            #insert command that will take the parameterized values and populate the subject table 
            sql_command = "INSERT OR IGNORE INTO Subject(SubjectID, Sex, Age, BMI, InsulinStatus) VALUES(?, ?, ?, ?, ?)"
            #connects to the method inside the EstablishDatabase class that carries out the insert commands 
            db.insert_db(sql_command, parameter_list)
    except FileNotFoundError:
        logger.error(f'The file {subject_file} cannot be opened. Please try again.')

'''
Within this function are the tables that will be used to represent the many to many relationship between peak 
and metabolite. The function begins by opening the metabolite annotation file, cleaning it, and writing it 
to a new file. This new file, with the clean data, is then opened, and in a similar convention to the subject 
data, is iterated through to create a list of tuples with the desired values. The peak data is opened and 
prepared for insertion into the peak table; however, there is the addition of the unique metabolite key 
in this loop, and the two unique keys for peak and metabolite are added to their own table. This successfully 
creates a link between the two tables. 

'''

def annotation_data(metabolite_annotation, peak_file):
    try:
        #first the metabolite data must be cleaned 
        #writing the new file that the cleaned data will be sent to 
        with open('metabolite_annotation_cleaned.txt', 'w') as clean:
            with open(metabolite_annotation) as annotation:
                #establishing headers for the new file 
                clean.write('PeakID\tMetabolite\tKEGG\tPathway\tMetaboliteID\n')
                #skip header 
                next(annotation)
                counts = 0
                #formatting the annotation data to clean it 
                for line in annotation:
                    line = line.rstrip().split(',')
                    peakID = line[0]
                    metabolite = line[1]
                    KEGG = line[2] 
                    pathway = line[5]
                    #counting if there exists more than one KEGG value 
                    #the count is based on '|' because it only ecists if there is more than one KEGG 
                    count = int(KEGG.count('|'))
                    #the format for duplicate metabolite names is (#), so this finds and replaces them with blank spaces
                    #this way, when queries are called, these names are treated as the same 
                    match = re.sub(r'\(\d\)', ' ', metabolite)
                    #a unique key is created for the metabolite annotation table 
                    counts += 1 
                    unique_key = f'{peakID}-{counts}'
                    #if the count of KEGG's is greater than 0 
                    if count > 0: 
                        #split the KEGG and metabolite values based on the divider 
                        #did not have to also count metabolites because if there is more than one KEGG, there will be 
                        #more than one metabolite that needs to also be separated 
                        KEGG_values = KEGG.split('|')
                        metabolites = metabolite.split('|')
                        #zip pairs up the two corresponding elements in each list 
                        for KEGG_hit, metabolite_value in zip(KEGG_values, metabolites):
                            #makes sure that neither is empty **
                            if KEGG_hit and metabolite_value:
                                #writes the newly cleaned data to the output file 
                                clean.write(f'{peakID}\t{metabolite_value}\t{KEGG_hit}\t{pathway}\t{unique_key}\n')
                    #the duplicate names already exist on their own row so the cleaned names just need to be added to the output file 
                    if match: 
                        #this file path inputs the cleaned metabolite names that matched with the regular expression above
                        clean.write(f'{peakID}\t{match}\t{KEGG}\t{pathway}\t{unique_key}\n')
                    #if the data passes all the checks, then no cleaning necessary and its added to the output file 
                    else: 
                        clean.write(f'{peakID}\t{metabolite}\t{KEGG}\t{pathway}\t{unique_key}\n')
    except FileNotFoundError:
        logger.error(f'The file {metabolite_annotation} cannot be opened. Please try again.')
                    
    try:
        #input the cleaned metabolite data 
        with open('metabolite_annotation_cleaned.txt') as metabolite_cleaned: 
            #skip the headers 
            next(metabolite_cleaned)
            #establishing an empty list for table data 
            parameter_list_1 = []
            #filtering through data 
            for info in metabolite_cleaned: 
                info = info.rstrip().split('\t')
                peakID = info[0]
                KEGG = info[2]
                Pathways = info[3]
                MetaboliteID = info[4]
                #establishing the structure of the parameter list 
                parameters = (MetaboliteID, peakID, KEGG, Pathways)
                #storing the metaboliteID in a dictionary so it can be used to create the junction table with peak data 
                parameter_list_1.append(parameters)
                #print(f"MetaboliteID: {MetaboliteID}, PeakID: {peakID}, KEGG: {KEGG}, Pathways: {Pathways}")

            #command that inserts the relevant data (stored in the parameter list) into the metabolite table 
            sql_command = "INSERT OR IGNORE INTO Metabolite(MetaboliteID, PeakID, KEGG, Pathways) VALUES(?, ?, ?, ?)"
            #calls the method inside the class to carry out the above command with the corresponding data 
            db.insert_db(sql_command, parameter_list_1)
    except FileNotFoundError:
        logger.error(f'The file {metabolite_cleaned} cannot be opened. Please try again.')

    try:
        #opening the peak data. this will end with the input of both the peak and peakannotation tables  
        with open(peak_file) as peak:
            #skip header
            next(peak)
            #initializing empty lists/values outside the loop 
            count = 0
            parameter_list_1 = [] 
            parameter_list_2 = []
            #formatting the file 
            for lines in peak:
                lines = lines.rstrip().split('\t')
                SampleID_peaks = lines[0]
                #splitting the sampleID to extract the VisitID information 
                subject, VisitID = SampleID_peaks.split('-')
                #creating a unique key for the peak table 
                count += 1 
                peaknumber = f'{subject}-{count}'
                #creating the two sets of parameters that will be used for the peak and peakannotation tables 
                parameters = (peaknumber, SampleID_peaks)
                parameter_list_1.append(parameters)
                params = (peaknumber, MetaboliteID)
                parameter_list_2.append(params)

            #commands for inserting the data into the peak and peakannotation table 
            #followed by connection to the method inside the class that will facilitate inputting the values into these tables 
            sql_link_table = "INSERT OR IGNORE INTO PeakAnnotation(PeakNumber, MetaboliteID) VALUES(?,?)"
            db.insert_db(sql_link_table, parameter_list_2)
        
            sql_command = "INSERT OR IGNORE INTO Peak(PeakNumber, SampleID) VALUES(?, ?)"
            db.insert_db(sql_command, parameter_list_1)
    except FileNotFoundError:
        logger.error(f'The file {peak_file} cannot be opened. Please try again.')

#final function for filtering and inputting data 
#open the transcript data file 
def transcript_data(transcriptome_data): 
    try:
        with open(transcriptome_data) as transcript: 
            #skip header
            next(transcript)
            #initialize empty list 
            parameter_list = []
            #format the file 
            for rows in transcript:
                rows = rows.rstrip().split('\t')
                SampleID_transcript = rows[0]
                #splitting the sampleID to extract the VisitID and subjectID
                subject, VisitID = SampleID_transcript.split('-')
                AIBG_abundance = rows[1]
                #establishing the organization of the tuples inside the parameter list 
                parameters = (SampleID_transcript, VisitID, AIBG_abundance, subject)
                #appending these values to the empty list outside the loop 
                parameter_list.append(parameters)
            
            #command to insert the data collected in the parameter list into the Sample table 
            sql_command = "INSERT OR IGNORE INTO Sample(SampleID, VisitID, TranscriptAbundance, SubjectID) VALUES(?, ?, ?, ?)"
            #connects to the method inside the class to carry out the insertion of data 
            db.insert_db(sql_command, parameter_list)
    except FileNotFoundError:
        logger.error(f'The file {transcriptome_data} cannot be opened. Please try again.')


#reading in the files and establishing variables 
db = EstablishDatabase(args.database_file)
subject_file = "Subject.csv"
peak_file = "HMP_metabolome_abundance.tsv"
metabolite_annotation = "HMP_metabolome_annotation.csv"
transcriptome_data = "HMP_transcriptome_abundance.tsv"

try:
    #depending on the command that is input at the command line, argparse will initiate the corresponding behavior in the script
    #this command creates the tables within the database 
    if args.createdb:
        #calls the get tables function 
        table_creation = get_tables()
        #pushes the data within the get tables function through to the method in the EstablishDatabase class that will facilitate making the tables
        db.create_tables(table_creation)

    #this command loads the data into the tables 
    if args.loaddb:
            #the file names are then passed through to their corresponding functions 
            subject_data(subject_file)
            transcript_data(transcriptome_data)
            annotation_data(metabolite_annotation, peak_file)

    '''
    The block of code below serves to establish each query based on the number called in the command line. 
    First the select statement is made, followed by the corresponding parameters. Then the result is obtained 
    by passing the query statement and its parameters into the query method from the EstablishDatabase class.
    The results are then pushed through the 'format_results' function, which prints the answer according to 
    the outlined specifications. 
    '''
    if args.querydb == 1:
        sql_query1 = "SELECT SubjectID, Age FROM Subject WHERE Age > ? AND Age NOT LIKE ?"
        params1 = (70, 'NA')
        result1 = db.query_db(sql_query1, params1)
        format_results(result1)
    elif args.querydb == 2:
        sql_query2 = "SELECT SubjectID FROM Subject WHERE Sex = ? AND BMI BETWEEN ? AND ? ORDER BY BMI DESC"
        params2 = ('F', 18.5, 24.9)
        result2 = db.query_db(sql_query2, params2)
        format_results(result2)
    elif args.querydb == 3:
        sql_query3 ="SELECT VisitID FROM Sample WHERE SubjectID = ?"
        params3 = ('ZNQOVZV',)
        result3 = db.query_db(sql_query3, params3)
        format_results(result3)
    elif args.querydb == 4:
        sql_query4 ="SELECT DISTINCT Subject.SubjectID FROM Subject, Sample, Peak WHERE Subject.SubjectID = Sample.SubjectID and Sample.SampleID = Peak.SampleID AND Subject.InsulinStatus = ?"
        params4 = ('IR',)
        result4 = db.query_db(sql_query4, params4)
        format_results(result4)
    elif args.querydb == 5:
        sql_query5 ="SELECT KEGG FROM Metabolite WHERE PeakID IN (? , ? , ?, ?)"
        params5 = ('nHILIC_121.0505_3.5', 'nHILIC_130.0872_6.3', 'nHILIC_133.0506_2.3', 'nHILIC_133.0506_4.4')
        result5 = db.query_db(sql_query5, params5)
        format_results(result5)
    elif args.querydb == 6:
        sql_query6 ="SELECT MIN(Age), MAX(Age), AVG(Age) FROM Subject WHERE Age NOT LIKE ?"
        params6 = ('NA',)
        result6 = db.query_db(sql_query6, params6)
        format_results(result6)
    elif args.querydb == 7:
        sql_query7 ="SELECT Pathways, COUNT(*) AS Annotations FROM Metabolite GROUP BY Pathways HAVING COUNT(*) >= ? ORDER BY Annotations DESC "
        params7 = (10,)
        result7 = db.query_db(sql_query7, params7)
        format_results(result7)
    elif args.querydb == 8:
        sql_query8 ="SELECT MAX(TranscriptAbundance) From Sample WHERE Sample.SubjectID = ?"
        params8 = ('ZOZOW1T',)
        result8 = db.query_db(sql_query8, params8)
        format_results(result8)
    elif args.querydb == 9:
        sql_query9 ="SELECT Age, BMI FROM Subject WHERE Age IS NOT NULL AND Age NOT LIKE ? AND BMI IS NOT NULL AND BMI NOT LIKE ?"
        params9 = ('NA','NA')
        result9 = db.query_db(sql_query9, params9)
        format_results(result9)
        #query 9 also has a scatterplot that is made when initialized from the command line 
        #turms the data from the query results into a dataframe
        query9_data = pd.DataFrame(result9)
        #making the columns for the plot 
        query9_data.columns = ['Age','BMI']
        #combining the above variables to create the scatterplot 
        sns.scatterplot(data=query9_data, x='Age',y='BMI')
        #saves the scatterplot under this file name 
        plt.savefig('age_bmi_scatterplot.png')
    else: 
        pass 
finally:
    #ensures the connection to the database is closed
    db.connection.close()
