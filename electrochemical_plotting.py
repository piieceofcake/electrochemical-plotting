import os
import glob
import numpy as np
import functions as fun
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import scipy.fftpack as fft
import pandas as pd
from IPython.display import display


"""
Electrochemical visualization tool 2024

This code is for plotting .csv and .mpt files extracted from Neware and Biologic battery testers, respectively.
Change the parameters underneath to fit your data

GitHub 1st commit 13.12.2024
    update 1: 22.09.2025

Contributors:
	Amalie Skurtveit (amalie.skurtveit@kjemi.uio.no)
    Casper Skautvedt (casper.skautvedt@smn.uio.no)
    Erlend Tiberg North (e.t.north@smn.uio.no)
    Magnus Lid Andresen
"""

### Standard format of figures ##
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "figure.figsize": "12.8, 4.8",  # figure size in inches. This is twice as large as normal (6.4, 4.8)
    "figure.constrained_layout.use": True,
    "lines.markersize": 9, # marker size, in points. 6 standard. 9 for QCM
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.bottom": True,
    "xtick.top": False,
    "xtick.minor.visible": True,
    "ytick.left": True,
    "ytick.right": False,
    "ytick.minor.visible": True,
    "xtick.major.size": 6, # major tick size in points (3.5 default)
    "xtick.minor.size": 4, # minor tick size in points (2 default)
    "ytick.major.size": 8, # major tick size in points (3.5 default)
    "ytick.minor.size": 6, # minor tick size in points (2 default)
    "xtick.major.width": 1.5,
    "xtick.minor.width": 1.5,
    "ytick.major.width": 1.5,
    "ytick.minor.width": 1.5
})

### USER SPECIFIC PARAMETERS ###
ELCHEM_FOLDER   = r"PATH TO YOUR FOLDER  HERE"
ELCHEM_FILETYPE = "csv"
MIN_VOLTAGE     = 0            # User defined, lower limit of voltage, e.g. 0 for anode
MAX_VOLTAGE     = 2            # User defined, upper limit of voltage, e.g., 2 for anode
ION             = "Na"

PLOT_ONLY_NEW   = False        # True if you don't want to plot already existing data in a folder
CATHODE         = False        # Full cell or cathode = True, Anode = False
CV              = False        # Are you plotting CV data?
COUNT_ION       = False        # Would you like to calculate how many IONS are transferred during (dis)charge?
REPORT_CPC      = False        # Save a .csv file for capacity per cycle plotting

CYCLES          = [1,2,3,5,10,20,50,100]
CPC_MAX         = 500        # Right limit of x-axis in capacity per cycle plot (CPC)
GC_PLOT_RIGHT   = 800        # Right limit of x-axis in GC plot
LINEWIDTH       = 2          # Thickness of line in GC plot
DPI             = 300        # dpi of the figure 
CO_BOTTOM       = 70         # Coulombinc efficiency (CO) plot minimum value on y-axis (should be as high as possible)
CO_TOP          = 101        # Coulombinc efficiency (CO) plot maximum value on y-axis (101, 105 is typically OK values)

ntick           = 4          # How many tick marks do you want. 4 Gives 5 tick marks. Ticks are the number marks on the axis.



DEST_FOLDER = ELCHEM_FOLDER.replace(f'{ELCHEM_FOLDER}', f'{ELCHEM_FOLDER}\\Figures')

if not os.path.exists(DEST_FOLDER):
    os.makedirs(DEST_FOLDER)


# Making the color palette for plotting GC curves
# You can change colors to whatever you would like as long as you have a hex code '#000000'
color_palette = ["#991045", "#E95135", "#F99959", "#57B1AB", "#387CB7", "#524295", "#272046"]
# Cycling through the colors automatically
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=color_palette) 
cycle = [x - 1 for x in CYCLES]
    

colors = color_palette
cpc_color = "k" # "k" = black, other colors: https://matplotlib.org/stable/gallery/color/named_colors.html


def read_data(file):
    """
    Function to read datafile
    """
    if datatype == "Neware1":
        voltage, specific_capacity, cycle_index, dqdv, step_index, time, Chg, DChg, current, coulombic_efficiency = fun.process_neware1_data(file)
        cyc1 = fun.cycle_neware1(voltage, current, CATHODE)

        GC_cpc_plot_N1(voltage, specific_capacity, cyc1, cycle_index, DChg, Chg, coulombic_efficiency)
                    
                
    elif datatype == "Neware2":
        df = fun.file_to_df(file)
        if CATHODE:
            df = fun.replace_cycle_index_N2(df)
        GC_cpc_plot_N2(df)

       
    elif datatype == "Biologic":
        df = fun.read_biologic(file, cv=CV, eis=False)
        if CV:
            CV_plot(df)
        else:
            ACTIVE_MASS = float(input(f"Active mass for {BAT_ID} in g:"))
            GC_cpc_plot_Biologic(df, ACTIVE_MASS)
            
           
    ## Comment in/out for plotting of dQdV data from N1 and N2
    # if datatype == "Neware1" or "Neware2":
    #     valid_cycles = [cycle for cycle in CYCLES if cycle <= len(cyc1) and cycle > 0]
    #     plot_dqdv(voltage, dqdv, valid_cycles, cyc1)
    
    if COUNT_ION:
        plot_moles_ion(df)

def CV_plot(df):
    """
    Function to plot CV data from a dataframe
    """
    plt.rcParams.update({
    "figure.figsize": "6.2, 4.8"
    })
    fig, cv = plt.subplots(constrained_layout=True)

    for index in CYCLES:
        df_cycle = df.loc[df["cycle number"] == index]
        #df_cycle = df_cycle.where(df_cycle["<I>/mA"] > 0)

        cv.plot(df_cycle.loc[:, "Ewe/V"],
                df_cycle.loc[:, "<I>/mA"], label=index)
    

    cv.set(xlabel=f'Voltage (V)', ylabel= 'Current (mA)')
    cv.set_xlim(left=MIN_VOLTAGE, right=MAX_VOLTAGE)
    cv.set_ylim(bottom=-2, top=2)
    
    cv.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    cv.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    
   
    plt.legend(bbox_to_anchor=(0.99, 1.03), loc=1, title='Cycle', frameon= False)
   
    
    cv.locator_params(nbins=ntick)

    fig.savefig(f'{DEST_FOLDER}/{BAT_ID}_CV_plot.png', dpi=DPI, bbox_inches='tight')


def GC_cpc_plot_N2(df):   
    """
    Galvanostatic cycling plot, i.e., voltage vs. specific capacity
    """
    plt.rc('lines', lw= LINEWIDTH)
    

    fig, (gc, cpc) = plt.subplots(1, 2, constrained_layout=True)
    step_dictionary = fun.neware_step_counter(df) #not necessary
    
    #######  GC  #######
    # Splitting the dataframe into cycles and plotting only these cycles
    for index in CYCLES:
        if CATHODE: #update 22.09.2025, fixing the cathode bug
            df_cycle = df.loc[df["Updated Cycle Index"] == index]
            df_cycle = df_cycle.where(df_cycle["Spec. Cap.(mAh/g)"] > 0)
            gc.plot(df_cycle.loc[:, "Spec. Cap.(mAh/g)"], 
                    df_cycle.loc[:, "Voltage(V)"], linewidth=LINEWIDTH, label=index)
        else:
            df_cycle = df.loc[df["Cycle Index"] == index]

            df_cycle = df_cycle.where(df_cycle["Spec. Cap.(mAh/g)"] > 0)
            gc.plot(df_cycle.loc[:, "Spec. Cap.(mAh/g)"], 
                    df_cycle.loc[:, "Voltage(V)"], linewidth=LINEWIDTH, label=index)
    
    gc.set(xlabel='Specific capacity (mAh $\\mathdefault{g^{-1}}$)', 
           ylabel= 'Voltage (V)')
           
    
    gc.locator_params(nbins=ntick)
    gc.set_xlim(left=0, right=GC_PLOT_RIGHT)
    gc.set_ylim([MIN_VOLTAGE, MAX_VOLTAGE])
    gc.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    gc.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    
    if CATHODE:
        gc.legend(bbox_to_anchor=(0.99, 1.03), loc=1, title='Cycle', frameon= False)
    else:
        gc.legend(bbox_to_anchor=(0.99, 1.03), loc=1, title='Cycle', frameon= False)

    gc.spines['bottom'].set_linewidth(2); gc.spines['left'].set_linewidth(2); gc.spines['top'].set_linewidth(2); gc.spines['right'].set_linewidth(2)


    #######  CPC #######
    n_cycles = df["Cycle Index"].iloc[-1]
    
    DChg = [] # this is specific capacity just discharge
    Chg = [] # this is specific capacity just charge

    ### for better calculation of Coulombic efficiency
    #cap_DChg = [] # capacity for discharge
    #cap_Chg = []  # capacity of charge

    resting = True if df["Step Type"].iloc[0] == "Rest" else False

    for index in range(1, n_cycles + 1):        #index from 1 to n
        df_cycle = df.loc[df["Cycle Index"] == index]
        #display(df_cycle)

        ### Step index is not based on if its a reduction/oxidation process ongoing, but rather if it's the 1st, 2nd, 3rd, ..., nth step in the program. Hence, have to check if the battery is resting
        if resting:
            try:
                df_step2 = df_cycle.loc[df_cycle["Step Index"] == 2] 
                df_step3 = df_cycle.loc[df_cycle["Step Index"] == 3]
                DChg.append(df_step2["Spec. Cap.(mAh/g)"].iloc[-1])
                Chg.append(df_step3["Spec. Cap.(mAh/g)"].iloc[-1])

                #cap_DChg.append(df_step2["Capacity(mAh)"].iloc[-1])
                #cap_Chg.append(df_step3["Capacity(mAh)"].iloc[-1])

            except IndexError:
                try:
                    df_step5 = df_cycle.loc[df_cycle["Step Index"] == 5] 
                    df_step6 = df_cycle.loc[df_cycle["Step Index"] == 6] 
                    DChg.append(df_step5["Spec. Cap.(mAh/g)"].iloc[-1])
                    Chg.append(df_step6["Spec. Cap.(mAh/g)"].iloc[-1])
                except IndexError:
                    continue
        else:
            try:
                df_step2 = df_cycle.loc[df_cycle["Step Index"] == 1] 
                df_step3 = df_cycle.loc[df_cycle["Step Index"] == 2]
                DChg.append(df_step2["Spec. Cap.(mAh/g)"].iloc[-1])
                Chg.append(df_step3["Spec. Cap.(mAh/g)"].iloc[-1])

                #cap_DChg.append(df_step2["Capacity(mAh)"].iloc[-1])
                #cap_Chg.append(df_step3["Capacity(mAh)"].iloc[-1])
            except IndexError:
                continue


    ### WARNING: Unpretty part of the code below ### 
    # Checking if DChg and Chg are the same length
    # We would get troubles of plotting the data because we often export data before it is complete
            
    len_DChg = len(DChg)
    len_Chg = len(Chg)

    if len_Chg >= len_DChg:
        total_number_of_cycles = len_DChg
    elif len_Chg <= len_DChg:
        total_number_of_cycles = len_Chg
        DChg.pop(-1)
    elif len_Chg == len_DChg:
        total_number_of_cycles = len_DChg
    else:
        total_number_of_cycles = len_DChg


    print("Total number of cycles", total_number_of_cycles)    
    cycle_number = [x for x in range(1, total_number_of_cycles + 1)]

    if CATHODE:
        cpc.scatter(cycle_number, Chg, c="k", marker="o")
        cpc.scatter(cycle_number, DChg, edgecolors="k", facecolors="None", marker="o")
        coulombic_efficiency = [(charge/discharge)*100 for charge, discharge in zip(Chg, DChg)]
    else:
        cpc.scatter(cycle_number, DChg, c="k", marker="o")
        cpc.scatter(cycle_number, Chg, edgecolors="k", facecolors="None", marker="o")
        coulombic_efficiency = [(charge/discharge)*100 for charge, discharge in zip(Chg, DChg)]

    if REPORT_CPC:
        df_report = pd.DataFrame()
        df_report["Cycle No."] = (pd.Series(cycle_number)).values
        if CATHODE:
            df_report["Spec.cap Chg (mAhg-1)"] = (pd.Series(DChg)).values # Chg capacity is encoded in the list called DChg.
            df_report["Spec.cap DChg (mAhg-1)"] = (pd.Series(Chg)).values # DChg capacity is encoded in the list called Chg.
            df_report.to_csv(f"{DEST_FOLDER}/{BAT_ID}_CPC_output.csv") 
        else:
            df_report["Spec.cap DChg (mAhg-1)"] = (pd.Series(DChg)).values
            df_report["Spec.cap Chg (mAhg-1)"] = (pd.Series(Chg)).values
            df_report.to_csv(f"{DEST_FOLDER}/{BAT_ID}_CPC_output.csv") 


    co = cpc.twinx()
    co.scatter(cycle_number, coulombic_efficiency, c="g", marker="^")

    co.set_ylim([CO_BOTTOM, CO_TOP]); co.set_ylabel("Coulombic efficiency (%)", c="g")
    co.locator_params(nbins=ntick)
    co.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    co.tick_params(axis='y', colors='green')
    cpc.locator_params(nbins=5)

    cpc.set_xlim(left=0, right=CPC_MAX)
    cpc.set_xlabel("Cycle number")   
    cpc.set_ylim(bottom=0, top=GC_PLOT_RIGHT)
    cpc.set_ylabel('Specific capacity (mAh $\\mathdefault{g^{-1}}$)')
    cpc.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    cpc.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    cpc.spines['bottom'].set_linewidth(2); cpc.spines['left'].set_linewidth(2); cpc.spines['top'].set_linewidth(2); cpc.spines['right'].set_linewidth(2)

    fig.savefig(f'{DEST_FOLDER}/{BAT_ID}_GC_cpc.png', dpi=DPI, bbox_inches='tight')
    print("GC and cpc plotted")
    plt.show()


def GC_cpc_plot_N1(voltage, specific_capacity, cyc1, cycle_number, discharge, charge, coulombic_efficiency):
    """
    Galvanostatic cycling plot, i.e., voltage vs. specific capacity
    """
    plt.rc('lines', lw= LINEWIDTH)

    
    voltage= np.array(voltage); voltage = np.ma.masked_outside(voltage, (MIN_VOLTAGE-0.001), (MAX_VOLTAGE+0.001)) #masking data outside of the voltage window to suppress lines 
    specific_capacity = np.array(specific_capacity); specific_capacity = np.ma.masked_where(specific_capacity < 0.05, specific_capacity)
   
    fig, (ax1, ax2) = plt.subplots(1, 2, constrained_layout=True)

    for i in range(len(cyc1)-1):
        if CYCLES != "all":
            if i in cycle:
                ax1.plot(specific_capacity[cyc1[i]:cyc1[i+1]], voltage[cyc1[i]:cyc1[i+1]], label=f"{i+1}")
        else:
            ax1.plot(specific_capacity[cyc1[i]:cyc1[i+1]], voltage[cyc1[i]:cyc1[i+1]])


    ax1.set(xlabel='Specific capacity (mAh $\\mathdefault{g^{-1}}$)', ylabel= f'Potential vs {ION}/'f'{ION}''$\\mathdefault{^+}$(V)')
    ax1.set_xlim(left=0, right=GC_PLOT_RIGHT)
    ax1.set_ylim([MIN_VOLTAGE, MAX_VOLTAGE])
    
    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    
    if CYCLES != 'all':
        ax1.legend(bbox_to_anchor=(0.99, 1.03), loc=1, title='Cycle', frameon= False)
    else:
        ax1.legend(bbox_to_anchor=(0.99, 1.03), loc=1, frameon=False)
    
    ax1.locator_params(nbins=ntick)

    ### cpc specific plotting ###

    """
    Function to plot the specific capacity per cycle number plot for Neware1 data, along with coulombic efficiency
    """
  
    co = ax2.twinx()

    ax2.scatter(-100, 100, c="k", marker="o")

    for i in range(len(cycle_number)):
        ax2.scatter(cycle_number[i], charge[i], edgecolors=cpc_color, facecolors="None", marker="o")
        ax2.scatter(cycle_number[i], discharge[i], c=cpc_color, marker="o")
      
        co.plot(cycle_number[i], coulombic_efficiency[i], c="g", marker="^")
    

    co.set_ylim([80, 105]); co.set_ylabel("Coulombic efficiency (%)", c="g")
    co.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    co.tick_params(axis='y', colors='green')
    co.locator_params('y', nbins=ntick)


    ax2.set_xlim(left=0, right=CPC_MAX); ax2.set_xlabel("Cycle number (#)")    
    ax2.set_ylim(bottom=0, top=GC_PLOT_RIGHT); ax2.set_ylabel('Specific capacity (mAh $\\mathdefault{g^{-1}}$)')
    ax2.locator_params('x', nbins=ntick) 
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    fig.savefig(f'{DEST_FOLDER}/{BAT_ID}_GC_cpc.png', dpi=DPI, bbox_inches='tight')
    print("GC and cpc plotted")
    #plt.show()


def GC_cpc_plot_Biologic(df, active_mass):
    """
    Function to plot GC and CPC plot for Biologic data
    """

    ### NB! only does capacity per now -- will need to find a function to do math on an entire column and save to a new column (should be quite straightforward!)


    fig, (gc, cpc) = plt.subplots(1, 2, constrained_layout=True)


    for index in CYCLES:
        df_cycle = df.loc[df["cycle number"] == index]
        capacity = np.array(df_cycle["Capacity/mA.h"].tolist())
        specific_capacity = capacity/active_mass
        #df_cycle = df_cycle.where(df_cycle["<I>/mA"] > 0)

        gc.plot(specific_capacity,
                df_cycle.loc[:, "Ewe/V"], label=index)
    

    gc.set(xlabel=f'Potential vs {ION}/'f'{ION}''$\\mathdefault{^+}$(V)', ylabel= 'Current (mA)')
    gc.set_xlim(left=0, right=GC_PLOT_RIGHT)
    gc.set_ylim(bottom=MIN_VOLTAGE, top=MAX_VOLTAGE)
    
    gc.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    gc.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    
   
    plt.legend(bbox_to_anchor=(0.99, 1.03), loc=1, title='Cycle', frameon= False)
   
    
    gc.locator_params(nbins=ntick)

    n_cycles = df['cycle number'].iloc[-1]

    DChg = []
    Chg = []

    for index in range(1, n_cycles + 1):
        df_cycle = df.loc[df['cycle number'] == index]
        df_red = df_cycle.loc[df_cycle['ox/red'] == 0]
        df_ox = df_cycle.loc[df_cycle['ox/red'] == 1]

        DChg.append(df_red["Capacity/mA.h"].iloc[-1])
        Chg.append(df_ox["Capacity/mA.h"].iloc[-1])

    cycle_number = [x for x in range(1, 1+n_cycles)]

    cpc.scatter(cycle_number, DChg, c="k", marker= "o")
    cpc.scatter(cycle_number, Chg, edgecolors="k", facecolors="None", marker="o")

    co = cpc.twinx()

    coulombic_efficiency = [(charge/discharge)*100 for charge, discharge in zip(Chg, DChg)]  
    co.plot(coulombic_efficiency, c="g", marker="^")

    co.set_ylim([70, 105]); co.set_ylabel("Coulombic efficiency (%)", c="g")
    co.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    co.tick_params(axis='y', colors='green')

    cpc.set_xlim(left=0, right=CPC_MAX)
    cpc.set_xlabel("Cycle number (#)")   
    cpc.set_ylim(bottom=0, top=GC_PLOT_RIGHT)
    cpc.set_ylabel('Specific capacity (mAh $\\mathdefault{g^{-1}}$)')
    cpc.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))


    fig.savefig(f'{DEST_FOLDER}/{BAT_ID}_GC_cpc_plot.png', dpi=DPI, bbox_inches='tight')


def smooth_dqdv(voltage, dqdv):
    """
    Smoothing function for dqdv
    """
    # Calculate sample spacing
    sample_spacing = voltage[1] - voltage[0]
    #print(sample_spacing)


    low_cutoff_frequency  = 0
    high_cutoff_frequency = 0.7
    

    # Perform the Fast Fourier Transform (FFT)
    y_fft = fft.fft(dqdv)

    # Calculate the FFT frequencies
 
    fft_frequencies = fft.fftfreq(len(dqdv), sample_spacing)
    #print(fft_frequencies)
  
   
    # Create band-pass filter mask: True between low and high cutoff frequencies
    band_pass_mask = (fft_frequencies >= low_cutoff_frequency) & (fft_frequencies <= high_cutoff_frequency)
    
    # Apply the mask to the y_fft array to retain only the desired frequencies
    y_fft_bandpass = y_fft.copy()
    y_fft_bandpass[~band_pass_mask] = 0
    
    y_smoothed = fft.ifft(y_fft_bandpass).real
    

    # Perfrom the inverse FFT to the smoothed signal
    y_smoothed = fft.ifft(y_fft_bandpass).real
    #print(y_smoothed)
    
    return y_smoothed
  
def plot_dqdv(voltage, dqdv, valid_cycles, cyc1):
    """
    Function to plot (smoothed) dQdV data. The function calls on the smooth_dqdv described above.

    """    

    fig, ax = plt.subplots(figsize=(6.4, 4.8))

    for cycle_number in valid_cycles:
        # Convert cycle_number to 0-index
        i = cycle_number - 1
        
        start_index = cyc1[i]
        # Handle the last slice of data to the end of the list
        next_start_index = cyc1[i+1] if (i + 1) < len(cyc1) else len(dqdv)
        
        if start_index >= len(voltage) or start_index == next_start_index:
            print(f"dQdV: No data for cycle {cycle_number}. Skipping...")
            continue

        cycle_voltage = voltage[start_index:next_start_index]
        cycle_dqdv = dqdv[start_index:next_start_index]
        
        
        # Ensure that the slice is not empty
        if len(cycle_dqdv) == 0:
            print(f"Empty data slice for cycle {cycle_number}")
            continue
        
        #sample_spacing = np.mean(np.diff(voltage)) # Use the average spacing of the current cycle
        
        #print(sample_spacing)
        smoothed_dqdv = smooth_dqdv(cycle_voltage, cycle_dqdv)
        
        
        # Plot the smoothed data for the selected cycle
        ax.plot(cycle_voltage, smoothed_dqdv, label=f"{cycle_number}")
    
    ax.set(xlim=[MIN_VOLTAGE, MAX_VOLTAGE], xlabel= f'Potential vs {ION}/'f'{ION}''$\\mathdefault{^+}$(V)',
        ylabel= 'dQ $\\mathdefault{dV^{-1}}$ (mAh $\\mathdefault{Vg^{-1}}}$)')
    
    ax.locator_params(nbins=ntick)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2)); ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))  
    plt.legend(bbox_to_anchor=(0.99, 1.03), loc=1, title='Cycle', frameon= False)
    fig.savefig(f'{DEST_FOLDER}/{BAT_ID}_dQdV_plot.png', dpi=DPI, bbox_inches='tight')
    print("dQdV plotted")
    #plt.show()


def plot_moles_ion(df):
    """
    Function to plot voltage vs # ION
    """
    
    molar_mass = float(input(f"Molar mass for your material in g/mol:"))
    
    fig, ax = plt.subplots()

    for index in CYCLES:
        if datatype == "Biologic":
            active_mass = float(input(f"Active mass for {BAT_ID} in g:"))
            df_cycle = df.loc[df["cycle number"] == index]
            capacity = np.array(df_cycle["Capacity/mA.h"].tolist())
            specific_capacity = capacity/active_mass

    
            x_Na = fun.moles_ion(specific_capacity, molar_mass)
            x_Na = np.ma.masked_outside(x_Na, (x_Na.max()-0.1), (x_Na.min()- 10))

            ax.plot(x_Na, df_cycle.loc[:, "Ewe/V"], linewidth=2)

        elif datatype == "Neware2":
            df_cycle = df.loc[df["Cycle Index"] == index]
            df_cycle = df_cycle.where(df_cycle["Spec. Cap.(mAh/g)"] > 0)

            specific_capacity = np.array(df_cycle["Spec. Cap.(mAh/g)"].tolist())
            x_Na = fun.moles_ion(specific_capacity, molar_mass)
            x_Na = np.ma.masked_outside(x_Na, (x_Na.max()-0.1), (x_Na.min()- 10))

            ax.plot(x_Na, df_cycle.loc[:, "Voltage(V)"], linewidth=2)
    
    if datatype == "Neware1":
        """
        Neware1 data (i.e., old format) is weird, so have to to it like this...
        """
        voltage, specific_capacity, cycle_index, dqdv, step_index, time, Chg, DChg, current, coulombic_efficiency = fun.process_neware1_data(file)
        cyc1 = fun.cycle_neware1(voltage, current)
        specific_capacity = np.array(specific_capacity)

        for i in range(len(cyc1)-1):
            x_Na = (((specific_capacity))*molar_mass)/26400
    
   

        voltage = np.ma.masked_outside(voltage, (MAX_VOLTAGE-0.001), (MIN_VOLTAGE+0.001))
        specific_capacity = np.ma.masked_less(specific_capacity, 0.001)
        x_Na = np.ma.masked_outside(x_Na, (x_Na.max()-0.1), (x_Na.min()- 10))


        for i in range(len(cyc1)-1): 
            if i in CYCLES:
                ax.plot(x_Na[cyc1[i]:cyc1[i+1]], voltage[cyc1[i]:cyc1[i+1]], label=f"{i}", color="k", linewidth=2)
    
    print(f"Max number of {ION}", x_Na.max())

    ax.set(xlabel="Number of Na (#)", ylabel=f'Potential vs {ION}/'f'{ION}''$\\mathdefault{^+}$(V)')
    ax.set_xlim(0, 15)
    ax.set_ylim([MIN_VOLTAGE, MAX_VOLTAGE])
    ax.locator_params(nbins=4)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))


    fig.savefig(f"{DEST_FOLDER}/{BAT_ID}_x{ION}.png", dpi=600, bbox_inches="tight")
    #plt.show()
        




if __name__ == "__main__":
    os.chdir(rf"{ELCHEM_FOLDER}")

    BAT_FILES = sorted(glob.glob(f'*{ELCHEM_FILETYPE}'))
    number_of_files = len(BAT_FILES)
    print("Number of battery files in folder to be plotted is: ", number_of_files)
 

    for i, file in enumerate(BAT_FILES):
        
        print("Currently reading: ", file)
        datatype, BAT_ID = fun.datatype(file)
        print(f"{BAT_ID} is from {datatype}")
        print("----------")

        if PLOT_ONLY_NEW:
            if os.path.isfile(fr"{DEST_FOLDER}/{BAT_ID}_GC_cpc.png"): 
                print(f"{BAT_ID} already plotted, moving on to the next file")
            else:
                read_data(file)
        else:
            read_data(file)
        

    
        
       
        plt.close()
        

   

        
        
print("--------------")
print("---Job done---")
print("--------------")  


    



