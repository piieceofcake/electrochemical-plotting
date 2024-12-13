"""
A function library to work with electrochemical_plotting.py script from GitHub

In your main .py file, make sure that the following line is included:
    Import functions as fun

Run functions in another Python script by calling them from here,
e.g. fun.neware_version(YOUR_FILE)

--------------------------------------------------------
GitHub 13.12.2024

Contributors:
	Amalie Skurtveit (amalie.skurtveit@kjemi.uio.no)
---------------------------------------------------------
"""

############################################################
"""
SECTION1: Functions to run with electrochemical plotting
"""
############################################################

### DATATYPE CHECK AND BATTERY ID ###
def datatype(file, elchem=True):
    """
    Function to check datatype from a variety of files
    
    Supported for both electrochemical plotting and X-ray based plotting
    
    NB! Electrochemical data collected at Neware1 and Neware2 after January2023 are equivalent and will be defined as Neware2  

    Parameters
    --------------
    file : path or path-like file (electrochemistry) or list (X-ray)
    elchem : bool | default True | (electrochemistry)

    Returns 
    --------------
    datatype : str | options : "Neware2", "Biologic", Neware1", "XRD", "EXAFS", "XANES"
    battery_id : str | optional | variable

    """
    if elchem:
        if file.endswith(".csv"):
            battery_id = file.split(".csv")[0]
        
            with open(file, 'r') as f:
                lines = f.readlines()
                try:
                    float(lines[2].split(',')[0])
                    datatype = "Neware2"
                except ValueError:
                    datatype = "Neware1"
            return datatype, battery_id

        elif file.endswith(".mpt"):
            battery_id = file.split(".mpt")[0]
            datatype = "Biologic"
        
        return datatype, battery_id
        
    elif elchem == False:
        if file[0].endswith(".nor"):
            datatype = "XANES"
        elif file[0].endswith(".chir"):
            datatype = "EXAFS"
        elif file[0].endswith(".xye"):
            datatype = "XRD"
        return datatype


    else:
        raise Exception("Filetype not supported by this program")
    

    

### Cyclic voltammetry and galvanostatic cycling from Biologic ###
def read_biologic(data, cv=False, eis=False):
    """
    Function to read and apply a filter to Biologic electrochemical

    Parameters
    ----------
    data : path or path-like file, must have extention ".mpt"
    cv : bool | default False
    eis : bool | default False
    
    Returns
    ---------
    df : filtered pandas dataframe
    """
    import pandas as pd
    
    ### Passing a filter on the dataframe to only extract interested columns
    FILTER_CV = [
        'ox/red',
        'Ewe/V',
        '<I>/mA',
        'cycle number']
    
    FILTER_EIS = [
        'mode',
        'ox/red',
        'Ewe/V',
        'z cycle',
        'Re(Z)/Ohm',
        '-Im(Z)/Ohm',
        'cycle number',
        'Capacity/mA.h',
        'Efficiency/%'
    ]
    
    FILTER_GC = [
        'ox/red',
        'Ewe/V',
        'cycle number',
        'Capacity/mA.h',
        'Efficiency/%',
        'time/s'
    ]
    
    with open(data, 'rb') as f:
        lines = f.readlines()
    header_lines = int(lines[1].split()[-1]) - 1   

    df = pd.read_csv(data, sep='\t', skiprows=header_lines, encoding="cp1252")
    if cv == True:
        df = df.filter(FILTER_CV)
    elif eis == True:
        df = df.filter(FILTER_EIS)
    else:
        df = df.filter(FILTER_GC)

    return df

### NEWARE2 in pandas DataFrame ####
def file_to_df(data):
    """
    Function to import and apply a filter to Neware electrochemical files

    Parameters
    ----------
    data : path or path-like file

    Returns
    ---------
    df : filtered pandas dataframe 
    """
    import pandas as pd

    # Filtering the data based on columns that we're interested in
    FILTER_NEWARE2 = [
        'Cycle Index', 
        'Step Index',
        'Step Type',
        'Voltage(V)', 
        'Spec. Cap.(mAh/g)', 
        'dQm/dV(mAh/V.g)']
    

    df = pd.read_csv(data, sep=",", engine="python")
    df = df.filter(FILTER_NEWARE2)
   
    return df

def moles_ion(specific_capacity:list, 
              molar_mass:int, 
              F=26400):
    """
    Calculation of number of ions electrochemically transferred during (dis)charge
    
    Parameters
    ------------
    specific_capacity : list, unit mAhg-1
    molar_mass : int, unit gmol-1
    F : constant, unit mAhmol-1
    
    Returns
    -------------
    List of calculated number of ion transferred 
    """
    return (specific_capacity*molar_mass)/F

def neware_version(file):
    """
    Function to check which Neware your data is from
    param: file: .csv file with Neware data. 
            New Neware1 files reads as Neware2.
    """
    with open(file, 'r') as f:
        lines = f.readlines()
        try:
            float(lines[2].split(',')[0])
            datatype = "Neware2"
        except ValueError:
            datatype = "Neware1"
    

    return datatype



### OLD FUNCTIONS THANT MIGHT BE IN USE IN SOME SCRIPTS, INCLUDED IN CASE OF TROUBLES ###
def process_neware1_data(data):
    """
    Function to process Neware data
    """
    Ewe  = []        #Voltage
    C    = []        #Specific capacity running
    Cy   = []        #Cycle
    dqdv = []        #dQdV
    sind = []        #step index: Chg or DChg
    t    = []        #time
    Chg  = []        #Charge specific cap
    DChg = []        #Discharge specific cap
    mA   = []        #Current
    CE   = []        #COULOMBic efficency


    with open(data, 'r') as f:
        file = f.readlines()
        
            
    for line in file[3:]:
        words = line.split(',\t')
        if len(words) == 23:                        #Header 3
            Ewe.append(float(words[3]))
            C.append(float(words[7]))
            dqdv.append(float(words[20]))
            t.append((words[2]))
            mA.append(float(words[4]))
        if len(words) == 25:                        #Header 2
            sind.append((words[2]))
        if len(words) == 30:                        #Header 1
            Cy.append(float(words[0]))
            Chg.append(float(words[3]))
            DChg.append(float(words[4]))
            CE.append(float(words[5]))
    
 
    return Ewe, C, Cy, dqdv, sind, t, Chg, DChg, mA, CE
   
def process_neware2_data(data):
    """
    Function to process Neware data
    """
    import datetime

    Ewe  = []        #Voltage
    C    = []        #Specific capacity running
    Cy   = []        #Cycle
    dqdv = []        #dQdV
    sind = []        #step index: Chg or DChg
    t    = []        #time
    SI   = []        #Step index

    with open(data, 'r') as f:
        file = f.readlines()
  
    for line in file[2:]:
        column = line.split(',')
        Cy.append((column[1]))
        SI.append(column[2])
        sind.append((column[3]))
        t.append(column[5])
        Ewe.append(float(column[7]))
        C.append(float(column[9]))

        # bear in mind that somehting has happened to Neware1 data 
        # recently, so that dqdv data might be in either column 22 or 23
        dqdv.append(float(column[22]))
        
        # possbility of implementing this::::
        # ## 
        # date_time_str = column[20]  
        # date_time_obj = datetime.datetime.strptime(date_time_str, '%d.%m.%Y %H:%M')

        # # Check for a specific date and time
        # specific_date_time = datetime.datetime(2024, 4, 16, 12, 8)  # Replace this with the date and time you want to check
        # is_specific_date_time = date_time_obj == specific_date_time

        # # Use the is_specific_date_time variable in your code to proceed accordingly
        # if is_specific_date_time:
            
            
        # else:
        #     dqdv.append(float(column[23]))

    

        

    return Ewe, C, Cy, dqdv, sind, t, SI

def cycle_index(step_index, voltage, cathode=False):
    """
    Function to do some thingys with Neware2 data, because Neware2 data is screwed
    """
    cy = []
    cyc1 = []
    cy1 = 0


    for i in range(len(step_index)):
        if i !=0:
            cy.append(cy1)  
        if cathode==False:  
            if step_index[i-1] != step_index[i]:   
                if step_index == "Rest":
                    cy1 = 0
                if step_index[i] == "CC DChg":
                    cyc1.append(i)
                    cy1 +=1
                if step_index[i] == "CCCV DChg":
                    cyc1.append(i)
                    cy1 += 1   
        else:
            if step_index[i-1] != step_index[i]:
                if step_index == "Rest":
                    cy1 = 0
                if step_index[i] == "CC Chg":
                    cyc1.append(i)
                    cy1 +=1
                if step_index[i] == "CCCV Chg":
                    cyc1.append(i)
                    cy1 += 1  
            
    cyc1.append(len(voltage))

    return cyc1

def cycle_neware1(voltage, current, cathode=False):
    """
    Functon to return cycles for N1 data
    """
    cyc1 = []
    
    for i in range(len(voltage)):
        if cathode:
            if current[i] >= 0:
                if current[i-1] < 0:
                    cyc1.append(i)
        else:
            if current[i] <= 0:
                if current[i-1] > 0:
                    cyc1.append(i)
        
           
    cyc1.append(len(voltage))
    #print(cyc1)
    return cyc1

def cycle_neware2(capacity, sind):
    """
    Function to split cycles from Neware2 data correctly
    param Ewe = voltage
    param C = specific capacity
    param sind = Chg or DChg step index
    param SI = step index
    Returns DC (DCHg) and CH (Chg)
    """
   
    DChg_capacity = []
    Chg_capacity = []
    
                
    for i in range(len(sind)):          #Loop to split the cycles correct
        
        if sind[i-1] in ('CC DChg', 'CCCV DChg') and sind[i] in ('CC Chg', 'CCCV Chg'):
            DChg_capacity.append(float(capacity[i-1]))
     
        if sind[i-1] in ('CC Chg', 'CCCV Chg') and sind[i] in ('CC DChg', 'CCCV DChg'):
            Chg_capacity.append(float(capacity[i-1]))
            if i == 0:
                Chg_capacity.pop()
    
    return DChg_capacity, Chg_capacity

def number_of_cycles(discharge_N2, charge_N2):
    """
    Function to make a list of the total number of cycles in N2 data
    """
    discharge = len(discharge_N2)
    charge = len(charge_N2)

    if charge >= discharge:
        total_number_of_cycles = discharge
    elif charge <= discharge:
        total_number_of_cycles = charge
    elif charge == discharge:
        total_number_of_cycles = discharge
    else:
        total_number_of_cycles = discharge
            
    ### making a list of cycle number information based on the length of the discharge list (i.e. total number of cycles)
    cycle_number = [i for i in range(1, total_number_of_cycles+1)] 

    return cycle_number

def process_biologic_data(data, mass):
    """
    Function to process biologic data.
    """
    Ewe     = []            #Voltage
    Cap     = []           #Capacity
    S       = []            #cycle
    T       = []            #Time
    oxred   = []            #Oxidation or reduction
    Cou     = []            #COULOMBic efficiency
    

    with open(data, 'r') as f: 
        file = f.readlines()
        for line in file:
            if 'Nb header lines' in line:
                headerlines = eval(line.split(':')[-1])
                break
        # for line in file[headerlines-1::-1]:    #[0:headerlines]: Looks from bottom and up, due to Biologic adding changes in the bottom
        #     if 'Characteristic mass' in line:
        #         Am = float(line.split(': ')[-1].split(' ')[0].replace(',', '.'))
        #         print('Active mass:',Am)
        #         break
        # if Am != mass:
        #     Am = float(input("Enter the Characteristic mass in mg: "))
        #     #Am = mass
        lines = file[headerlines:]
        for line in lines:
            line = line.replace(',','.')
            words = line.split('\t')
            oxred.append(float(words[1]))
            T.append(float(words[7])) 
            Ewe.append(float(words[11]))
            Cap.append(float(words[23])) #changed from 25
            Cou.append(float(words[20]))
            S.append(float(words[-1]))
    
    return Ewe, Cap, S, T, oxred, Cou

def cycle_biologic(cycle, voltage, oxred):
    """
    Function to split and detect cycles for Biologic data
    """
    # Ts = []
    # tel = 0
    ### splitting cycles correctly based on an index
    Cy = [0]
    x = 0

    for i in range(len(cycle)):
        if cycle[i] != x:
            Cy.append(i)
            x = cycle[i]
        # tel += 30
        # Ts.append(tel)

    Cy.append(len(voltage))


    # split on reduction, i.e. counting how many cycles there are in total#
    S1 = []
    y = 0
    Cy1 = [0]
    for i in range(len(voltage)):
        S1.append(y)
        if oxred[i-1] == 1:
            if oxred[i] == 0:
                y += 1
                Cy1.append(i)
                if voltage[i] > voltage[i+1]:
                    voltage[i] = 0
    
    Cy1.append(len(voltage))

    return Cy, Cy1

  
def format_GC(gc_file, batsmall=False):
    """
    Function to extract potential and time from the GC file and saving only these extracted values to a file
    """
    import numpy as np

    if batsmall:
        start = 3
    else:
        start = 0

    time = []
    voltage = []
    with open(gc_file, "r") as f:
        data = f.readlines()[start:] # Data is from defined start and out
        for line in data:
            words = line.split('\t')
            time.append(float(words[0])*3600)
            voltage.append(float(words[1]))

    #time = np.array([t/3600 for t in time])
    #time = np.array(time)
    voltage = np.array(voltage)

    return np.array([t/3600 for t in time]), voltage 


##############################################################################
"""
SECTION2: Functions to run with plotting of X-ray based files
"""
##############################################################################

### X-RAY BASED FILES TO DATAFRAME ### UPDATE2024: CAN HANDLE MORE DATATYPES
def X_files_to_df(files, datatype=None, two_d=False):
    """
    Function to import X-ray based files into a pandas dataframe
    
    Parameters
    ---------
    files : filepath_or_buffer
    datatype : str | default None | options: "XRD", "EXAFS", "XANES"

    Returns
    -------
    df : pandas dataframe; index is either 2theta, FT magnitude, or energy
    """
    import pandas as pd
    from IPython.display import display
    import numpy as np

    df = pd.DataFrame()
    column_names = []

    if datatype == "XRD":
        delimiter = " "
        skip = 0
        n = 1
    elif datatype == "EXAFS":
        delimiter = "\s+"
        skip = 37
        n = 3
    elif datatype == "XANES":
        delimiter = "\s+"
        skip = 37
        n = 1
    else: 
        datatype = "XRD"

    if two_d:
        df = pd.read_csv(files, sep=delimiter, skiprows=skip)

        filter_xanes = ["#","e"]
        filter_exafs = ["#", "r"]
        filter_xrd = []

        if datatype =="XANES": df.filter(filter_xanes)
        elif datatype == "EXAFS": df.filter(filter_exafs)
    else:
        for i, file in enumerate(files):
            column_names.append(i)
            data = pd.read_csv(file, sep=delimiter, skiprows=skip, dtype="float64")
            display(data)
            data.index = data.iloc[:, 0]
            data_of_interest = data.iloc[:, n]
            df = pd.concat([df, data_of_interest], axis=1)
            df.columns = column_names
        

        #data.index = data.iloc[:, 0]
        #df.astype("float64", copy=False)
    return df



def convert_xye(files, folder, id):
    """
    Convert dat files into xye files
    """
    import os
    count = 0
    for filename in files:
        old_name = filename
        new_name = folder + f"\{id}_{str(count).zfill(5)}.xye"
        os.rename(old_name, new_name)
        count += 1
    
def differential(df, col):
    """
    Function for differential XRD or PDF 
    param df: pandas dataframe of scans
    param col: int, which column should be subtracted from each scan
    Returns a differential pandas dataframe, df_differential
    """
    import pandas as pd

    df_differential = pd.DataFrame()
    num_columns = df.shape[1]

    for columns in range(0, num_columns):
        diff_col = df.iloc[:, columns] - df.iloc[:, col]
        #df_differential = pd.concat([df_differential, diff_col], axis=1)
        df_differential[f'{columns}'] = diff_col
        
    return df_differential
