#### CHEM5042 ASSIGNMENT 2 ####
#StudentID: 2915700 

#importing libraries 
import sqlite3
import tkinter 
import tkinter as tk 
from tkinter import ttk 
from tkinter.ttk import *
from tkinter import *
from tkinter import PhotoImage
import rdkit 
from rdkit import Chem 
from rdkit.Chem import Descriptors, AllChem, Lipinski
from rdkit.Chem import Draw
import logging
import os 
import PIL
from PIL import Image, ImageTk
import io 
import csv


#create logger to help with error handling 
logger = logging.getLogger()
#selecting logging level 
logger.setLevel(logging.ERROR)
#create stream handler, define the format, add this to logger 
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(sh)

#### PARSE FILE 
#function for finding fused aromatic rings 
def finding_fused_aromatics(molecule):
    # Calculate aromaticity
    fused_aromatics_count = 0
    #only considering compounds that have 2 or more rings 
    if rings >= 2: 
        #iterates over each ring in the molecule 
        for ring in range(0, rings):
            #conditional statement checks if the bonds in the ring are aromatic 
            #bond index identifies coordinates for the bonds and then determines if the bond is aromatic 
            if all(molecule.GetBondWithIdx(bond_index).GetIsAromatic() for bond_index in molecule.GetRingInfo().BondRings()[ring]):
                #if the bond is aromatic it is then checked if its fused 
                if molecule.GetRingInfo().IsRingFused(ring) is True:
                    #if both conditions are met, then its added to the count of fused aromatic rings for that molecule 
                    fused_aromatics_count += 1
    return fused_aromatics_count 

#function for calculating if a molecule adheres to the lipinski rules 
def lipinski_rule_5(molwt,logP, hbd, hba):
    #basic conditional that includes the rules 
    if molwt <= 500 and logP <= 5 and hbd <= 5 and hba <= 10: 
        #boolean statement returned (which is converted to 0 or 1 in SQL)
        return True 
    else:
        return False 
    
#function for determining lead likeness for each molecule 
def lead_likeness(molwt, logD, rings, rotatable_bonds, hbd, hba):
    if molwt <= 450 and rings <= 4 and rotatable_bonds <= 10 and hbd <= 5 and hba <= 8 and -4 <= logD <= 4:
        return True
    else:
        return False 
    
#function that determines if the molecule in the current iteration is bioactive 
def bioavailability(molwt, logP, hbd, hba, rotatable_bonds, psa, aromatic_fused_count):
    bioavailability_count = 0 
    if molwt <= 500: 
        bioavailability_count += 1 
    if logP <= 5:
        bioavailability_count += 1 
    if hbd <= 5:
        bioavailability_count += 1 
    if hba <= 10: 
        bioavailability_count += 1 
    if rotatable_bonds <= 10: 
        bioavailability_count += 1 
    if psa <= 200:
        bioavailability_count += 1 
    if aromatic_fused_count <= 5:
        bioavailability_count += 1 
    #if 6 or more of the above conditions are met, then it is counted as bioactive 
    if bioavailability_count >= 6:
        return True 
    else:
        return False 
    
'''function that converts the smiles for each molecule into the actual image, which is then saved to a folder as a PNG file, 
but also stored as binary data for inserting into the SQL database'''
def molecule_structure(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    image = Draw.MolToImage(mol, mols_per_row=2, img_size=(200, 200))

    #convert the image to binary data (this is what is going to be stored in the database)
    image_binary = io.BytesIO() 
    image.save(image_binary, format='PNG')
    binary_data = image_binary.getvalue()

    #saving the structures to a file for easy referencing and error catching 
    output_path = os.path.join(output_folder, filename)
    image.save(output_path)

    #return the binary data for that molecule to be stored in the database 
    return(binary_data)

#creates the folder that will house the images (this is not actually used, just for QC purposes)
try:
    output_folder = "molecular_structures"
    os.mkdir(output_folder)
    logger.info(f'Folder {output_folder} sucessfully created. Located in: {os.getcwd()}')
except FileExistsError:
    logger.error('File already exists. Rename folder or continue with existing file.')


#initiating an empty list to be populated inside the loop 
parameter_list = [] 
#using rdkit to parse the sdf file 
molecules_file = Chem.ForwardSDMolSupplier('Molecules21.sdf')
for  molecule in molecules_file: 
    try:
        #error catch 
        if molecule is None: 
            logger.error('Loading molecule failed. Skipping.')
            ''''efonidipine hydrochloride is the only problematic molecule in my list.'''
            continue 
        else:
            #extracting cdid 
            cdid = molecule.GetProp("CdId")
            #extract the name 
            name = molecule.GetProp("Name")
            #formula 
            formula = molecule.GetProp("Formula")
            #molID
            molID = molecule.GetProp("Mol_ID")
            #extracting smiles 
            smiles = Chem.MolToSmiles(molecule)
            #extracting molecular weight 
            molwt = Descriptors.ExactMolWt(molecule)
            molwt = float("{:.2f}".format(molwt))
            #logp
            logP = Chem.Crippen.MolLogP(molecule)
            logP = float("{:.2f}".format(logP))
            #logD
            logD = float(molecule.GetProp("LogD"))
            #PSA
            psa = Chem.rdMolDescriptors.CalcTPSA(molecule)
            psa = float("{:.2f}".format(psa))
            #h-bond acceptors  
            hba = int(Chem.rdMolDescriptors.CalcNumLipinskiHBA(molecule))
            #h-bond donors 
            hbd = int(Chem.rdMolDescriptors.CalcNumLipinskiHBD(molecule)) 
            #calculate rings
            rings = int(molecule.GetRingInfo().NumRings())
            #rotatable bonds 
            rotatable_bonds = int(Chem.rdMolDescriptors.CalcNumRotatableBonds(molecule)) 
            #now calling the above functions to fill in the rest of the parameters needed
            #fused aromatic rings 
            aromatic_fused_count = finding_fused_aromatics(molecule)
            #lipinski's rule of 5 
            lipinski_status = lipinski_rule_5(molwt,logP, hbd, hba)
            #lead likeness 
            lead_likeness_status = lead_likeness(molwt, logD, rings, rotatable_bonds, hbd, hba)
            #bioavailability 
            bioavailability_status = bioavailability(molwt, logP, hbd, hba, rotatable_bonds, psa, aromatic_fused_count)
            #2D structure and BLOB
            #BLOB doesnt show up in the database -- this is computationally expensive so for the purposes of efficiency, say if this was a larger database, 
            #its probably best to just store the data and then call it back for the GUI to actually be viewed
            img_filename = f"{molecule.GetProp('CdId')}.png"
            structure = molecule_structure(smiles, img_filename)
            #taking the structures and getting binary data 
            #creating a list of desired info per molecule 
            parameters = (cdid, name, formula, molID, smiles, molwt, logP, logD, psa, hba, hbd, lipinski_status, lead_likeness_status, bioavailability_status, rings, rotatable_bonds, aromatic_fused_count, structure)
            #appending to the list outside the loop, creating a tuple 
            parameter_list.append(parameters)
    except AttributeError:
        logger.error(f'Error in {molecule.GetProp("Name")}. Skipping.')
        continue 


#### PUT INTO DATABASE 
#create and connect to database 
db = 'CHEM5042_Assignment.db'
try: 
    connection = sqlite3.connect(db)
    cur = connection.cursor() 
except sqlite3.OperationalError as error: 
    logger.error(f'Unable to establish connection to {db}: {error}')
#creating the table for data insertion 
try:
    entity = "CREATE TABLE IF NOT EXISTS Molecules(CdId INT PRIMARY KEY, Name TEXT, Formula VARCHAR(50), MoleculeID VARCHAR(50), Smiles VARCHAR(50), Mol_Weight DECIMAL(5,2), LogP DECIMAL(2,2), LogD DECIMAL(2,2), PSA INT, H_Bond_Acceptors INT, H_Bond_Donors INT, Lipinski_5 BOOLEAN, Lead_Likeness BOOLEAN, Bioavailability BOOLEAN, Rings INT, Rotatable_Bonds INT, Fused_Aromatics INT, Structure BLOB)"
    cur.execute(entity)
except sqlite3.DatabaseError as error: 
    logger.error(f'Error during table creation: {error}')
#inserting the data into the table 
try: 
    sql_command = "INSERT OR IGNORE INTO Molecules VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    cur.executemany(sql_command, parameter_list)
    connection.commit() 
except sqlite3.IntegrityError as error:
    logger.error(f"SQLite Integrity Error: {error}")

#### GUI
#establishing GUI 
#defining basic color scheme that will be applied to all elements of the GUI 
main_background_color = "#163020"
secondary_background_color = "#B6C4B6"
accent_text_color = "#163020"
#establishing the main window of the GUI 
window= tkinter.Tk()
window.title('Chemical Structure Database')
#setting parameters for the width and length of the GUI and its placement on the user's screen 
window.configure(background=main_background_color)
window.minsize(1300, 800)
window.maxsize(1300, 800)
window.geometry("1300x800+50+50")
#this is necessary for fitting the database display to the screen 
window.columnconfigure(1,weight=1)
window.rowconfigure(0, weight=1)
#basic aesthetic formatting
theme = ttk.Style()
theme.theme_use('default')
#used to ensure the images in the database can fit into the rows theyve been assigned to 
theme.configure('Treeview', rowheight=150, columnwidth= 50)

################## LEFT FRAME 

#establishing frames inside the GUI 
left_frame = Frame(window, width=100, height=600, background=secondary_background_color)
#placing frame inside the main window 
left_frame.grid(row=0, column=0, padx=0, pady=5)
#adding a window for filtering on this frame
selection_window = Frame(left_frame, width=100, height= 600, background=secondary_background_color)
selection_window.grid(row = 2, column=0, padx=5, pady=5)
#labels inside the filter window 
left_label = Label(left_frame, text="Filter Parameters", fg=accent_text_color, font=("Helvetica", 18), background=secondary_background_color).grid(row=0, column=0, padx=5, pady=5)
Label(selection_window, text="Molecular Weight",fg=accent_text_color,background=secondary_background_color).grid(row=2, column=0, padx=5, pady=5)
Label(selection_window, text="LogP",fg=accent_text_color,background=secondary_background_color).grid(row=3, column=0, padx=5, pady=5)
Label(selection_window, text="LogD",fg=accent_text_color,background=secondary_background_color).grid(row=4, column=0, padx=5, pady=5)
Label(selection_window, text="PSA",fg=accent_text_color,background=secondary_background_color).grid(row=5, column=0, padx=5, pady=5)
Label(selection_window, text="Acceptor Count",fg=accent_text_color,background=secondary_background_color).grid(row=6, column=0, padx=5, pady=5)
Label(selection_window, text="Donor Count",fg=accent_text_color,background=secondary_background_color).grid(row=7, column=0, padx=5, pady=5)
Label(selection_window, text="Ring Count",fg=accent_text_color,background=secondary_background_color).grid(row=8, column=0, padx=5, pady=5)
Label(selection_window, text="Rotatable Bond Count",fg=accent_text_color,background=secondary_background_color).grid(row=9, column=0, padx=5, pady=5)
Label(selection_window, text="Fused Aromatic Count",fg=accent_text_color,background=secondary_background_color).grid(row=10, column=0, padx=5, pady=5)
Label(selection_window, text="Select by Column", fg=accent_text_color,background=secondary_background_color).grid(row=0, column=0, padx=5, pady=5)
Label(selection_window, text="Enter Range:",fg=accent_text_color,background=secondary_background_color).grid(row=1, column=0, padx=5, pady=5)
Label(selection_window, text="Minimum Value",fg=accent_text_color,background=secondary_background_color).grid(row=1, column=1, padx=5, pady=5)
Label(selection_window, text="Maximum Value", fg=accent_text_color,background=secondary_background_color).grid(row=1, column=2, padx=5, pady=5)
#making a search bar (last minute)
#Label(selection_window, text="Search by Name:", fg=accent_text_color,background=secondary_background_color).grid(row=14, column=0, padx=5, pady=5)


#establishing empty variables that will be used to flag if a check button has been selected or not 
var1 = tk.IntVar()
var2 = tk.IntVar()
var3 = tk.IntVar()
#making buttons for the booleans 
#lipinski
lipinksi_button = Checkbutton(selection_window, text="Lipinski Rule of 5", variable=var1, onvalue= 1, offvalue=0, fg=accent_text_color,background=secondary_background_color)
lipinksi_button.grid(row=11, column=1, padx=5, pady=5)
#lead likeness
lead_button = Checkbutton(selection_window, text="Lead Likeness",variable=var2, onvalue= 1, offvalue=0, fg=accent_text_color,background=secondary_background_color)
lead_button.grid(row=11, column=2, padx=5, pady=5)
#bioavailability 
bioavailability_button = Checkbutton(selection_window, text="Bioavailability",variable=var3, onvalue= 1, offvalue=0,fg=accent_text_color,background=secondary_background_color)
bioavailability_button.grid(row=12, column=1, padx=5, pady=5)
#general button that will be used to apply all the filter selections at once 
apply_filter_button = Button(selection_window, text="Apply",fg=accent_text_color,background=secondary_background_color)
apply_filter_button.grid(row=13, column=2, padx=5, pady=5) 
#this will clear the filters and reset it to the default 
apply_reset_button = Button(selection_window, text="Clear",fg=accent_text_color,background=secondary_background_color)
apply_reset_button.grid(row=13, column=1, padx=5, pady=5) 

#dropdown menu to sort the database by the selected parameter
#string var stores the selection made by the user 
selected_option1 = StringVar()
options1 = ["None", "Mol_Weight", "LogP", "LogD", "PSA", "H_Bond_Acceptors", "H_Bond_Donors", "Rings", "Rotatable_Bonds", "Fused_Aromatics"]
dropdown_menu1 = OptionMenu(selection_window, selected_option1, *options1)
dropdown_menu1.grid(row=0, column=1, padx=5, pady=5) 
#order by ascending/descending for the column selected from option1
selected_option2 = StringVar()
options2 = ["None", "ASC", "DESC"]
dropdown_menu2 = OptionMenu(selection_window, selected_option2, *options2)
dropdown_menu2.grid(row=0, column=2, padx=5, pady=5) 

#creating entry options for filtering 
min_molwt = Entry(selection_window)
min_molwt.grid(row =2, column=1)
max_molwt = Entry(selection_window)
max_molwt.grid(row =2, column=2)

min_logp = Entry(selection_window)
min_logp.grid(row =3, column=1)
max_logp = Entry(selection_window)
max_logp.grid(row =3, column=2)

min_logd = Entry(selection_window)
min_logd.grid(row =4, column=1)
max_logd = Entry(selection_window)
max_logd.grid(row =4, column=2)

min_psa = Entry(selection_window)
min_psa.grid(row =5, column=1)
max_psa = Entry(selection_window)
max_psa.grid(row =5, column=2)

min_acceptor = Entry(selection_window)
min_acceptor.grid(row =6, column=1)
max_acceptor = Entry(selection_window)
max_acceptor.grid(row =6, column=2)

min_donor = Entry(selection_window)
min_donor.grid(row =7, column=1)
max_donor = Entry(selection_window)
max_donor.grid(row =7, column=2)

min_rings = Entry(selection_window)
min_rings.grid(row =8, column=1)
max_rings = Entry(selection_window)
max_rings.grid(row =8, column=2)

min_rotate= Entry(selection_window)
min_rotate.grid(row =9, column=1)
max_rotate = Entry(selection_window)
max_rotate.grid(row =9, column=2)


min_fused = Entry(selection_window)
min_fused.grid(row =10, column=1)
max_fused = Entry(selection_window)
max_fused.grid(row =10, column=2)

'''
#search bar entry
search_bar = Entry(selection_window)
search_bar.grid(row=14, column=1)

def search_by_name(cur):
    entry = search_bar.get()

    if entry == '':
       pull_from_db(cur, tree)
    else: 
        cur.execute(f"SELECT {entry} from Molecules")
        for row in cur.fetchall():
'''

############# RIGHT FRAME 

#making a label for the frame 
right_frame = Frame(window, width=1200, height=800, bg='#31304D')
right_frame.grid(row=0, column=1, sticky='nsew', padx=10, pady=5)

#putting the database viewer into the right frame of the GUI 
tree = ttk.Treeview(right_frame, style='Treeview', columns= ("Name", "Formula", "Mol_Weight", "LogP", "LogD", "PSA", "H_Bond_Acceptors", "H_Bond_Donors", "Lipinski_5", "Lead_Likeness", "Bioavailability", "Rings", "Rotatable_Bonds", "Fused_Aromatics"))
#have to write this as a for loop because the headings command adds the column header one at a time 
for col in tree["columns"]:
    tree.heading(col, text=col)

#creating a 'garbage' collection that stores the image data, since the images are only projected in the GUI, not stored 
tree.image_cache = []
#function to pull the information from DB 
def pull_from_db(cur, tree):
    cur.execute("SELECT * FROM Molecules")
    #iterating over each row from the results of the SQL command and pulling desired values 
    for row in cur.fetchall():
        image_data = row[-1]
        cdid = row[0]
        name = row[1]
        formula = row[2]
        moleculeID = row[3]
        smile = row[4]
        mol_weight = row[5]
        logp = row[6]
        logd = row[7]
        psa = row[8]
        h_bond_acceptors = row[9]
        h_bond_donors = row[10]
        lipinski_5 = row[11]
        lead_likeness = row[12]
        bioavailability = row[13]
        rings = row[14]
        rotatable_bonds = row[15]
        fused_aromatics = row[16]
        #converting the binary image data back to an image 
        convert_to_image = Image.open(io.BytesIO(image_data))
        #resizing the image to fit in the rows of the GUI 
        convert_to_image.thumbnail((150,150))
        #convert to tkinter photoimage 
        tk_image = (ImageTk.PhotoImage(convert_to_image))
        #storing the values pulled from above 
        values_for_gui = (name, formula, mol_weight, logp, logd, psa, h_bond_acceptors, h_bond_donors, lipinski_5, lead_likeness, bioavailability, rings, rotatable_bonds, fused_aromatics)
        #inserting these values into the tree 
        tree.insert('', 'end', image=tk_image, values=(values_for_gui))
        #updating the list outside the loop 
        tree.image_cache.append(tk_image)
        #updating the GUI display 
        window.update()
        #calculating the number of data entries as 'hits'
        number_hits = 0 
        for child in tree.get_children():
            number_hits += 1
        Label(selection_window, text=f"Hits: {number_hits}", fg=accent_text_color,background=secondary_background_color).grid(row=13, column=0, padx=5, pady=5)

#function to save the filter results to a csv file
query_number = 0 
def filter_to_file(table_values, number_hits):
    #to try to make the file name unique and to also make identifying the filter results a little easier, the query count is added to the filename 
    #made a global variable because its defined outside of the function, but i dont want to input it as a parameter of the function 
    global query_number 
    query_number += 1 
    filter_results_filename = f"filter_results_query{query_number}_hits{number_hits}.csv"
    try:
        #opening the file for data insertion 
        with open(filter_results_filename, 'w') as results_file:
            #writing the file with headers 
            write = csv.writer(results_file)
            headers = ["Name", "Formula", "Mol_Weight", "LogP", "LogD", "PSA", "H_Bond_Acceptors", "H_Bond_Donors", "Lipinski_5", "Lead_Likeness", "Bioavailability", "Rings", "Rotatable_Bonds", "Fused_Aromatics"]
            write.writerow(headers)
            #iterating over each row of values that are included in the filter parameters to be included in the csv file 
            for row in table_values:
                write.writerow(row)
        logger.info(f"Filter results saved to {filter_results_filename}")
    except FileExistsError:
        logger.error('File already exists. Rename file or continue with existing file.')

    
#function that will concatenate all the filter options into one SQL query
#another garbage collector for the tkinter images to be stored  
tree.image_garbage = []
#empty list to store the values of each row that applies to the outlined filter requirements. must be outside the loop to store/update properly 
filter_results_table_values = []
def apply_filters_button(cur, tree, min_molwt, max_molwt, min_logp, max_logp, min_logd, max_logd, min_psa, max_psa, min_acceptor, max_acceptor, min_donor, max_donor, min_rings, max_rings, min_rotate, max_rotate, min_fused, max_fused, var1, var2, var3, selected_option1, selected_option2):
    #making dynamic SQL query by appending the selection statement and the desired parameters to these empty lists 
    filter_conditions = []
    user_parameters = [] 

    #establishing default values that will include all the values for this dataset (so if nothing is selected it will select all essentially)
    #mol_weight
    default_min_molwt = 0
    default_max_molwt = 1000
    #will either store the user input values or the pre-defined default values 
    min_molwt = float(min_molwt.get()) if min_molwt.get() != '' else default_min_molwt
    max_molwt = float(max_molwt.get()) if max_molwt.get() != '' else default_max_molwt
    #appending the empty lists 
    filter_conditions.append("Mol_Weight BETWEEN ? AND ?")
    #using extend instead of append because adding multiple elements 
    user_parameters.extend([min_molwt,max_molwt])
    #same convention as above follows for all the valued parameters 
    #logP
    default_min_logp = -1000
    default_max_logp = 1000
    min_logP = float(min_logp.get()) if min_logp.get() != '' else default_min_logp
    max_logP = float(max_logp.get()) if max_logp.get() != '' else default_max_logp

    filter_conditions.append("LogP BETWEEN ? AND ?")
    user_parameters.extend([min_logP, max_logP])
    #logD
    default_min_logd = -1000
    default_max_logd = 1000
    min_logD = float(min_logd.get()) if min_logd.get() != '' else default_min_logd
    max_logD = float(max_logd.get()) if max_logd.get() != '' else default_max_logd

    filter_conditions.append("LogD BETWEEN ? AND ?")
    user_parameters.extend([min_logD, max_logD])
    #PSA
    default_min_psa = -1000
    default_max_psa = 1000
    min_PSA = float(min_psa.get()) if min_psa.get() != '' else default_min_psa
    max_PSA = float(max_psa.get()) if max_psa.get() != '' else default_max_psa

    filter_conditions.append("PSA BETWEEN ? AND ?")
    user_parameters.extend([min_PSA, max_PSA])
    #H_Bond_Acceptor Count 
    default_min_acceptor = 0 
    default_max_acceptor = 1000
    min_acceptor_count = float(min_acceptor.get()) if min_acceptor.get() != '' else default_min_acceptor
    max_acceptor_count = float(max_acceptor.get()) if max_acceptor.get() != '' else default_max_acceptor

    filter_conditions.append("H_Bond_Acceptors BETWEEN ? AND ?")
    user_parameters.extend([min_acceptor_count, max_acceptor_count])
    #H_Bond_Donor_Count 
    default_min_donor = 0 
    default_max_donor = 1000
    min_donor_count = float(min_donor.get()) if min_donor.get() != '' else default_min_donor
    max_donor_count = float(max_donor.get()) if max_donor.get() != '' else default_max_donor
    
    filter_conditions.append("H_Bond_Donors BETWEEN ? AND ?")
    user_parameters.extend([min_donor_count, max_donor_count])
    #Ring Count 
    default_min_ring = 0 
    default_max_ring = 1000
    min_ring_count = float(min_rings.get()) if min_rings.get() != '' else default_min_ring
    max_ring_count = float(max_rings.get()) if max_rings.get() != '' else default_max_ring

    filter_conditions.append("Rings BETWEEN ? AND ?")
    user_parameters.extend([min_ring_count, max_ring_count])
    #Rotatable Bonds Count 
    default_min_rotate = 0 
    default_max_rotate = 1000
    min_rotate_count = float(min_rotate.get()) if min_rotate.get() != '' else default_min_rotate
    max_rotate_count = float(max_rotate.get()) if max_rotate.get() != '' else default_max_rotate

    filter_conditions.append("Rotatable_Bonds BETWEEN ? AND ?")
    user_parameters.extend([min_rotate_count,max_rotate_count])
    #Fused Aromatics Count 
    default_min_fused = 0
    default_max_fused = 1000
    min_fused_aromatics = float(min_fused.get()) if min_fused.get() != '' else default_min_fused
    max_fused_aromatics = float(max_fused.get()) if max_fused.get() != '' else default_max_fused

    filter_conditions.append("Fused_Aromatics BETWEEN ? AND ?")
    user_parameters.extend([min_fused_aromatics,max_fused_aromatics])

    #convention is adjusted for the boolean parameters
    #the status of the checked boxes is stored (indicates if the user interacted with them or not)
    lipinski_condition = var1.get()
    lead_condition = var2.get()
    bioavailability_condition = var3.get()
    #if the buttons have been checked they will be stored as a '1' and therefore indicates that the user has selected them for filtering 
    if lipinski_condition == 1:
        #when selected for filtering the corresponding select statement is added to the overall SQL query 
        filter_conditions.append("Lipinski_5 = 1")
    if lead_condition == 1:
        filter_conditions.append("Lead_Likeness = 1")
    if bioavailability_condition == 1:
        filter_conditions.append("Bioavailability = 1")

    #this is the SQL query that will have been created from all the above information 
    base_sql_query = "SELECT * FROM Molecules WHERE " + " AND ".join(filter_conditions)

    #if the user decides to also display the results based on ASC/DESC of a specific parameter, then those conditions are also applied 
    #conditional statement checks that the user did in fact select a parameter from the database
    if selected_option1.get() != "" and selected_option1.get() != "None":
        #storing the parameter that was selected and if the order is by ASC or DESC 
        sort_column = selected_option1.get()
        order_by = selected_option2.get()
        
        dynamic_sql_query = f"{base_sql_query} ORDER BY {sort_column} {order_by}"
    #if the drop down menus were not interacted with, then the function proceedes unchanged 
    else: 
        dynamic_sql_query = base_sql_query

    #the SQL query is then executed 
    cur.execute(dynamic_sql_query, user_parameters) 
    filter_results = cur.fetchall()

    #get_children gets the list of all the items currently in the GUI and deletes them 
    for data in tree.get_children():
        tree.delete(data)
        window.update()

    #now iterates over each row that was generated from the results of the filter parameters
    for row in filter_results:
        cdid = row[0]
        name = row[1]
        formula = row[2]
        moleculeID = row[3]
        smile = row[4]
        mol_weight = row[5]
        logp = row[6]
        logd = row[7]
        psa = row[8]
        h_bond_acceptors = row[9]
        h_bond_donors = row[10]
        lipinski_5 = row[11]
        lead_likeness = row[12]
        bioavailability = row[13]
        rings = row[14]
        rotatable_bonds = row[15]
        fused_aromatics = row[16]
        image_data = row[-1]

        convert_to_image = Image.open(io.BytesIO(image_data))
        #resize 
        convert_to_image.thumbnail((150,150))
        #convert to tkinter photoimage 
        tk_image = (ImageTk.PhotoImage(convert_to_image))
        tree.image_garbage.append(tk_image)
        
        table_values = (name, formula, mol_weight, logp, logd, psa, h_bond_acceptors, h_bond_donors, lipinski_5, lead_likeness, bioavailability, rings, rotatable_bonds, fused_aromatics)
        tree.insert('', 'end', image=tk_image, values=(table_values))

        number_hits = 0 
        for child in tree.get_children():
            number_hits += 1
        #updates the hits specific to the filter specifications 
        Label(selection_window, text=f"Hits: {number_hits}", fg=accent_text_color,background=secondary_background_color).grid(row=13, column=0, padx=5, pady=5)
        
        #need to write the results outside the loop in order for the list to be written to the csv correctly
        filter_results_table_values.append(table_values)
    #calling the file output function to store the results of the hits for easy referencing 
    filter_to_file(filter_results_table_values, number_hits)

#when the apply button is clicked, the filter is called 
apply_filter_button.configure(command= lambda: apply_filters_button(cur, tree, min_molwt, max_molwt, min_logp, max_logp, min_logd, max_logd, min_psa, max_psa, min_acceptor, max_acceptor, min_donor, max_donor, min_rings, max_rings, min_rotate, max_rotate, min_fused, max_fused, var1, var2, var3, selected_option1, selected_option2))

#function to clear all the filters and reset the database 
def reset_filters_button(min_molwt, max_molwt, min_logp, max_logp, min_logd, max_logd, min_psa, max_psa, min_acceptor, max_acceptor, min_donor, max_donor, min_rings, max_rings, min_rotate, max_rotate, min_fused, max_fused, var1, var2, var3):
    min_molwt.delete(0, 'end')
    max_molwt.delete(0, 'end')
    min_logp.delete(0, 'end')
    max_logp.delete(0, 'end')
    min_logd.delete(0, 'end')
    max_logd.delete(0, 'end')
    min_psa.delete(0, 'end')
    max_psa.delete(0, 'end')
    min_acceptor.delete(0, 'end')
    max_acceptor.delete(0, 'end')
    min_donor.delete(0, 'end')
    max_donor.delete(0, 'end')
    min_rings.delete(0, 'end')
    max_rings.delete(0, 'end')
    min_rotate.delete(0, 'end')
    max_rotate.delete(0, 'end')
    min_fused.delete(0, 'end')
    max_fused.delete(0, 'end')

    var1.set(0)
    var2.set(0)
    var3.set(0)
    selected_option1.set("None")
    selected_option2.set("None")

    for data in tree.get_children():
        tree.delete(data)
        window.update()
    #calls the original function that makes the database (without any filters)
    return pull_from_db(cur, tree)
 
#calls the function above when the clear button is clicked 
apply_reset_button.configure(command=lambda: reset_filters_button(min_molwt, max_molwt, min_logp, max_logp, min_logd, max_logd, min_psa, max_psa, min_acceptor, max_acceptor, min_donor, max_donor, min_rings, max_rings, min_rotate, max_rotate, min_fused, max_fused, var1, var2, var3))


#inserting the database into the GUI using grid 
tree.grid(row=0, column=0, sticky='ns')

#creating a x scrollbar because i cant get the window and the database to work together 
scrollbar = Scrollbar(right_frame,orient= 'horizontal',command=tree.xview)
scrollbar.grid(row=1,column=0, sticky='ew')
tree.configure(xscrollcommand=scrollbar.set)
#creating a y scrollbar because i cant get the window and the database to work together 
scrollbar = Scrollbar(right_frame,orient= 'vertical',command=tree.yview)
scrollbar.grid(row=0,column=5, sticky='ns')
tree.configure(yscrollcommand=scrollbar.set)

#apparently it matters in what order you make the frame expandable...
right_frame.rowconfigure(0, weight=1)
right_frame.columnconfigure(0,weight=1)

#calling the DB function 
pull_from_db(cur, tree)
#applying everything to be generated into the GUI 
window.mainloop()
#closing connection to the DB 
connection.close() 
