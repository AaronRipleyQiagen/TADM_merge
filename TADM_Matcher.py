import os
import pandas as pd
import numpy as np
import warnings
import pytz
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfilenames, asksaveasfilename, askdirectory
import uuid
warnings.filterwarnings('ignore')

class UI(Frame):    
    def __init__(self, master=None):
        
        Frame.__init__(self, master)
        
        self.myParser = nmdx_file_parser()
        
        self.raw_data = pd.DataFrame()
        self.tadm_data = pd.DataFrame()

        self.bakFrame = tk.Frame(master, bg='white')
        self.bakFrame.place(relx=0, rely=0, relheight=1, relwidth=1, anchor='nw')

        self.GetDataButton = tk.Button(self.bakFrame, text="Select Raw Data Files", bg='white', command=self.load_raw_data)
        self.GetDataButton.place(relx=0.04, rely=0.05, anchor='nw', relwidth=0.2, relheight=0.25)

        self.GetTADMReferenceButton = tk.Button(self.bakFrame, text="Select TADM Reference Directory", bg='white', command=self.get_tadm_data)
        self.GetTADMReferenceButton.place(relx=0.28, rely=0.05, anchor='nw', relwidth=0.2, relheight=0.25)

        self.ProcessDataButton = tk.Button(self.bakFrame, text="Process Data", bg='white', command=self.process_data)
        self.ProcessDataButton.place(relx=0.52, rely=0.05, anchor='nw', relwidth=0.2, relheight=0.25)
        
        self.SaveDataButton = tk.Button(self.bakFrame, text="Save Data", bg='white', command=self.save_data)
        self.SaveDataButton.place(relx=0.76, rely=0.05, anchor='nw', relwidth=0.2, relheight=0.25)

        self.resetButton = tk.Button(self.bakFrame, text="Clear Data", bg='white', command=self.clearData)
        self.resetButton.place(relx=0.28, rely=0.35, anchor='nw', relwidth=0.48, relheight=0.20)

        ##Row 1 of Options
        self.includeXPCR = IntVar()
        self.includeXPCR.set(0)
        self.includeXPCRButton = tk.Checkbutton(self.bakFrame, text="Include XPCR Module Info", onvalue=1, offvalue=0, variable=self.includeXPCR, bg='white')
        self.includeXPCRButton.place(relx=0.025, rely=0.575, anchor='nw', relwidth=0.30, relheight=0.1)

        self.includeHM = IntVar()
        self.includeHM.set(0)
        self.includeHMButton = tk.Checkbutton(self.bakFrame, text="Include Heater Module Info", onvalue=1, offvalue=0, variable=self.includeHM, bg='white')
        self.includeHMButton.place(relx=0.35, rely=0.575, anchor='nw', relwidth=0.30, relheight=0.1)

        self.includeTS = IntVar()
        self.includeTS.set(0)
        self.includeTSButton = tk.Checkbutton(self.bakFrame, text="Include Test Strip Position Info", onvalue=1, offvalue=0, variable=self.includeTS, bg='white')
        self.includeTSButton.place(relx=0.675, rely=0.575, anchor='nw', relwidth=0.30, relheight=0.1)

        ##Row 2 of Options
        self.includeSP = IntVar()
        self.includeSP.set(0)
        self.includeSPButton = tk.Checkbutton(self.bakFrame, text="Include Sample Processing Info", onvalue=1, offvalue=0, variable=self.includeSP, bg='white')
        self.includeSPButton.place(relx=0.025, rely=0.70, anchor='nw', relwidth=0.30, relheight=0.1)

        self.includeCL = IntVar()
        self.includeCL.set(0)
        self.includeConsumableLotButton = tk.Checkbutton(self.bakFrame, text="Include Consumable Lot Info", onvalue=1, offvalue=0, variable=self.includeCL, bg='white')
        self.includeConsumableLotButton.place(relx=0.35, rely=0.70, anchor='nw', relwidth=0.30, relheight=0.1)

        self.includeCR = IntVar()
        self.includeCR.set(0)
        self.includeCRButton = tk.Checkbutton(self.bakFrame, text="Include Channel Result Info", onvalue=1, offvalue=0, variable=self.includeCR, bg='white')
        self.includeCRButton.place(relx=0.675, rely=0.70, anchor='nw', relwidth=0.30, relheight=0.1)

    def load_raw_data(self):
        self.ReadingLabels = tk.Label(self.bakFrame, text="Parsing Raw Data", bg='blue', fg='white')
        self.ReadingLabels.place(relx=0, rely=0.825, anchor='nw', relwidth=1, relheight=0.20)
        
        files = [('XLSX', '*.xlsx')] 

        files = askopenfilenames(filetypes = files, defaultextension = files)
        for file in files:
            self.raw_data = pd.concat([self.raw_data,self.myParser.scrapeFile(file=file, filename='test')])
        
        self.myTadmHelper = TadmHelper(self.raw_data)
        self.ReadingLabels.destroy()  
    def get_tadm_data(self):
        self.ReadingLabels = tk.Label(self.bakFrame, text="Parsing TADM Data", bg='blue', fg='white')
        self.ReadingLabels.place(relx=0, rely=0.825, anchor='nw', relwidth=1, relheight=0.20)
        self.bakFrame.update()
        try:
            file_dir = askdirectory()
            self.myTadmHelper.get_tadms(file_dir)
        except:
            print("Failed to collect TADM Data")
        self.ReadingLabels.destroy()
    def process_data(self):
        self.ReadingLabels = tk.Label(self.bakFrame, text="Matching TADM data with NMDX Data", bg='blue', fg='white')
        self.ReadingLabels.place(relx=0, rely=0.825, anchor='nw', relwidth=1, relheight=0.20)
        self.bakFrame.update()
        try:
            self.myTadmHelper.tadm_hunter()
        except:
            print("Failed to process data.")
        self.ReadingLabels.destroy()
    def save_data(self):
        print(self.includeXPCR, self.includeHM, self.includeTS, self.includeSP, self.includeCL, self.includeCR)
        self.ReadingLabels = tk.Label(self.bakFrame, text="Matching TADM data with NMDX Data", bg='blue', fg='white')
        self.ReadingLabels.place(relx=0, rely=0.825, anchor='nw', relwidth=1, relheight=0.15)
        self.bakFrame.update()
        try:
            tadm_output = self.myTadmHelper.tadm_merger(include_XPCR_info=self.includeXPCR.get(),include_HM_info=self.includeHM.get(),include_TS_info=self.includeTS.get(),include_SP_info=self.includeSP.get(),include_ConsLot_info=self.includeCL.get(), include_ChannelResult_info=self.includeCR.get())
            output_dir = asksaveasfilename(title="Choose where to save TADM Data", defaultextension=".xlsx", initialfile="TADM_output", filetypes=[("CSV", "*.csv")])
            tadm_output.to_csv(output_dir)
        except:
            print("Failed to save Data.")
        self.ReadingLabels.destroy()
    def clearData(self):
        self.ReadingLabels = tk.Label(self.bakFrame, text="Clearing Data", bg='blue', fg='white')
        self.ReadingLabels.place(relx=0, rely=0.825, anchor='nw', relwidth=1, relheight=0.20)
        self.bakFrame.update()
        try:
            self.raw_data = pd.DataFrame()
            del self.myTadmHelper
        except:
            print("Failed to clear Data from App.")
        self.ReadingLabels.destroy()

class nmdx_file_parser:
    """
    A class used to read raw data file(s) and convert to flat format.

    Methods
    -------
    scrapeFile(file=None, env=None)
        Scrapes data from one raw data file.
    """
    def __init__(self):
        self.file_data = {}

    def readChannelData(file, sheet, channel):

        channelData_all = pd.read_excel(io=file,sheet_name=sheet)
        if len(channelData_all) > 0:
            ChannelRawStart = channelData_all[channelData_all['Sample ID']=='Raw'].index.values[0] + 1
            ChannelRawEnd = channelData_all[channelData_all['Sample ID']=='Normalized'].index.values[0] - 2
            ChannelRaw = channelData_all.loc[ChannelRawStart:ChannelRawEnd]
            ChannelRaw['Processing Step'] = 'Raw'

            ChannelNormStart = channelData_all[channelData_all['Sample ID']=='Normalized'].index.values[0] + 1
            ChannelNormEnd = channelData_all[channelData_all['Sample ID']=='SecondDerivative'].index.values[0] - 2
            ChannelNorm = channelData_all.loc[ChannelNormStart:ChannelNormEnd]
            ChannelNorm['Processing Step'] = 'Normalized'

            Channel2ndStart = channelData_all[channelData_all['Sample ID']=='SecondDerivative'].index.values[0] + 1

            if 'Modulated' in channelData_all['Sample ID'].unique():
                Channel2ndEnd = channelData_all[channelData_all['Sample ID']=='Modulated'].index.values[0] - 2
                ChannelModulatedStart = channelData_all[channelData_all['Sample ID']=='Modulated'].index.values[0] + 1
                ChannelModulated = channelData_all.loc[ChannelModulatedStart:ChannelModulatedStart+len(ChannelRaw)]
                ChannelModulated['Processing Step'] = 'Modulated'
                Channel2nd = channelData_all.loc[Channel2ndStart:Channel2ndEnd]
                Channel2nd['Processing Step'] = '2nd'

                if len(ChannelRaw) == len(ChannelNorm) and len(ChannelRaw) == len(Channel2nd) and len(ChannelRaw) == len(ChannelModulated):

                    ChannelFinal = pd.concat([ChannelRaw, ChannelNorm, Channel2nd, ChannelModulated],axis=0)
                    ChannelFinal['Channel'] = channel
                    ChannelFinal.set_index(['Test Guid', 'Replicate Number'],inplace=True)
                else:
                    print("Error in parsing Datablocks")
            else:
                Channel2nd = channelData_all.loc[Channel2ndStart:Channel2ndStart+len(ChannelRaw)]
                Channel2nd['Processing Step'] = '2nd'
                #if len(ChannelRaw) == len(ChannelNorm) and len(ChannelRaw) == len(Channel2nd):
                ChannelFinal = pd.concat([ChannelRaw, ChannelNorm, Channel2nd],axis=0)
                ChannelFinal['Channel'] = channel
                ChannelFinal.set_index(['Test Guid', 'Replicate Number'],inplace=True)


        else:
            ChannelFinal = pd.DataFrame()



        return ChannelFinal
    
    def readRawData(file):
        channelDict = {'Green_470_510':'Green',
                    'Yellow_530_555':'Yellow',
                    'Orange_585_610':'Orange',
                    'Red_625_660':'Red',
                    'Far_Red_680_715':'Far_Red'}

        Summary_Tab = pd.read_excel(io=file,sheet_name='Summary',header=2)
        COC_Tab = pd.read_excel(io=file,sheet_name='Chain of Custody')
        Summary_COC_Data = Summary_Tab.set_index(['Test Guid', 'Replicate Number']).join(COC_Tab.set_index(['Test Guid', 'Replicate Number']).loc[:, [x for x in COC_Tab.columns if x not in Summary_Tab.columns]])


        channelDataDict = {}
        for channel in channelDict:
            channelDataDict[channel] = nmdx_file_parser.readChannelData(file, channel, channelDict[channel])
        channelDataFinal = pd.concat([channelDataDict[df] for df in channelDataDict if len(channelDataDict[df])>0],axis=0)

        
        channelDataFinal.set_index(['Target Result Guid', 'Processing Step', 'Channel'],append=True,inplace=True)
        for i in range(1,256):
            if "Readings "+ str(i) not in channelDataFinal.columns:
                channelDataFinal["Readings "+str(i)] = np.nan
        channelDataFinal_readings = channelDataFinal.loc[:, ['Readings '+str(i) for i in range(1,256)]]
        channelDataFinal_summary = channelDataFinal.swaplevel(3,0).swaplevel(3,1).swaplevel(3,2)
        channelDataFinal_summary = channelDataFinal_summary.loc['Raw'].drop(['Readings '+str(i) for i in range(1,256)],axis=1)

        return Summary_COC_Data, channelDataFinal_summary, channelDataFinal_readings
    
    def retrieveConsumableLots(data, consumable_types=['Pcr Cartridge', 'Capture Plate', 'Test Strip NeuMoDx', 'Buffer', 'Release Reagent', 'Wash Reagent']):
        """
        Retrieves Lot information for NMDX Consumables from Barcode String
        :param consumable_types: list-like List of Consumables to get Data For.
        """
    
        for consumable_type in consumable_types:
            data[consumable_type+" Lot"] = data[consumable_type+" Barcode"].str[18:24]

        return data

    def retrieveConsumableSerials(data, consumable_types=['Pcr Cartridge', 'Capture Plate', 'Test Strip NeuMoDx', 'Buffer', 'Release Reagent', 'Wash Reagent']):
        """
        Retrieves Consumable Serial information for NMDX Consumables from Barcode String
        :param consumable_types: list-like List of Consumables to get Data For
        """
        
        for consumable_type in consumable_types:
            data[consumable_type+" Serial"] = data[consumable_type+" Barcode"].str[27:32]

        return data

    def retrieveConsumableExpiration(data, consumable_types=['Pcr Cartridge', 'Capture Plate', 'Test Strip NeuMoDx', 'Buffer', 'Release Reagent', 'Wash Reagent']):
        """
        Retrieves Expiration Date information for NMDX Consumables from Barcode String
        :param consumable_types: list-like List of Consumables to get Data For.
        """
    
        for consumable_type in consumable_types:
            data[consumable_type+" EXP Date"] = data[consumable_type+" Barcode"].str[-6:].apply(lambda x: pd.to_datetime(arg=x, format="%y%m%d"))

        return data

    def getRawMinusBlankCheckReads(self, data):
        """
        A Function used to calculate the Difference between the First three Raw Readings and Blank Check Values for each target result included in dataset provided
        Parameters
        ----------
        data (pandas.DataFrame) = DataFrame to be used for Calculation.
        """
        RawReadsMinusBlankCheckFrame = data.reset_index()[['Processing Step', 'Test Guid', 'Replicate Number', 'Target Result Guid']+['Readings 1', 'Readings 2', 'Readings 3', 'Blank Reading']].copy()
        RawReadsMinusBlankCheckFrame.set_index(['Processing Step', 'Test Guid', 'Replicate Number', 'Target Result Guid'],inplace=True)
        RawReadsMinusBlankCheckFrame_Raw = RawReadsMinusBlankCheckFrame.loc['Raw']
        RawReadsMinusBlankCheckFrame_Raw['Blank Check - 1st 3 Reads'] = RawReadsMinusBlankCheckFrame_Raw[['Readings 1', 'Readings 2', 'Readings 3']].mean(axis=1) - RawReadsMinusBlankCheckFrame_Raw['Blank Reading']
        RawReadsMinusBlankCheckFrame = RawReadsMinusBlankCheckFrame.join(RawReadsMinusBlankCheckFrame_Raw[['Blank Check - 1st 3 Reads']])
        data['Blank Check - 1st 3 Reads'] = RawReadsMinusBlankCheckFrame['Blank Check - 1st 3 Reads'].values
    
    def channelParametersFlattener(self, data, stats=['Target Name', 'Localized Result', 'Ct', 'End Point Fluorescence', 'EPR', 'Max Peak Height', 'Baseline Slope', 'Baseline Y Intercept', 'Blank Check - 1st 3 Reads']):
        """
        Retrieves Channel Specific stats and returns them all channels in one-dimmensional column.
        stats:  Which Stats to flatten.
        """
        channel_stats = data.reset_index().drop_duplicates(['Test Guid', 'Channel', 'Replicate Number']).set_index(['Test Guid', 'Replicate Number']).loc[:, stats+['Channel']]

        channel_stats = channel_stats.reset_index().pivot(columns='Channel',values=stats,index=['Test Guid', 'Replicate Number'])
        channel_stats.columns = [y+" "+x for (x,y) in channel_stats.columns]
        data = data.reset_index().set_index(['Test Guid', 'Replicate Number']).join(channel_stats)
        return data

    def scrapeFile(self, file, filename):
           
        #time = pd.Timestamp.now()

        summary_coc, channel_summary, channel_readings = nmdx_file_parser.readRawData(file)
        for col in channel_summary.columns:
            if 'Barcode' in col:
                channel_summary[col] = channel_summary[col].astype(str)
                channel_summary[col] = channel_summary[col].str.replace("_x001D_", " ")
        channel_summary = channel_summary.astype(object).where(pd.notna(channel_summary), None)


        for col in summary_coc.columns:
            if 'Barcode' in col:
                summary_coc[col] = summary_coc[col].astype(str)
                summary_coc[col] = summary_coc[col].str.replace("_x001D_", " ")
            if 'ADP Position' in col:
                summary_coc[col] = summary_coc[col].astype(str)
        summary_coc = summary_coc.astype(object).where(pd.notna(summary_coc), None)
        for col in summary_coc.loc[:, [col for col in summary_coc if 'Date' in col]].columns:
            summary_coc[col] = pd.to_datetime(summary_coc[col], utc=False).apply(lambda x: x.replace(tzinfo=pytz.utc))
        
        channel_readings = channel_readings.astype(object).where(pd.notna(channel_readings), None)

        channel_summary['File Source'] = filename
        channel_readings['File Source'] = filename
        summary_coc['File Source'] = filename
        summary_coc.rename({'Flags':'Summary Flags'},axis=1,inplace=True)
        channel_summary.rename({'Flags':'Channel Flags'},axis=1,inplace=True)
        summary_coc = nmdx_file_parser.retrieveConsumableLots(summary_coc)
        summary_coc = nmdx_file_parser.retrieveConsumableSerials(summary_coc)
        summary_coc = nmdx_file_parser.retrieveConsumableExpiration(summary_coc)

         
        
        

        flat_data = summary_coc.join(channel_summary.loc[:, [x for x in channel_summary.columns if x not in summary_coc.columns]]).join(channel_readings.loc[:, [x for x in channel_readings.columns if x not in channel_summary.columns]])
        self.getRawMinusBlankCheckReads(flat_data)
        flat_data = self.channelParametersFlattener(flat_data)
        ##Add Target Result / Localized Result columns if not in flat_data columns
        if 'Localized Result' not in flat_data.columns:
            flat_data['Localized Result'] = np.nan
        
        if 'Target Result' not in flat_data.columns:
            flat_data['Target Result'] = np.nan

        return flat_data.reset_index()

class TadmHelper:

    def __init__(self, raw_data):
        """
        Parameters
        ----------
        raw_data pd.DataFrame: A Raw Data DataFrame in Flat format.
        """
        self.tadm_data = pd.DataFrame()
        self.raw_data = raw_data.copy()
        self.channels = sorted(raw_data['Channel'].unique()) 
        self.raw_data.drop_duplicates(subset=['Test Guid', 'Replicate Number'], inplace=True)
        self.raw_data_liquid_handle_processes = self.raw_data[['Test Guid', 'Replicate Number', 'Sample ID', 'Start Date/Time', 'LHPA Start Date Time', 'LHPB Start Date Time', 'LHPC Start Date Time', 'PCR Start Date Time', 'LHPA ADP Position', 'LHPB ADP Position', 'LHPC ADP Position']]

        self.processGroups = {}
        for process in ['LHPA', 'LHPB', 'LHPC']:
            self.processGroups[process] = self.raw_data[[process+' Start Date Time']].drop_duplicates([process+' Start Date Time']).dropna().sort_values(process+' Start Date Time').reset_index(drop=True)

    def get_tadms(self, file_dir):
        """
        A function used to prepare a the tadm_data by merging together data found within a directory.

        Parameters
        ----------
        file_dir (str): Name of File Directory to Search for Files within.
        """
        
        attempt_id = uuid.uuid4()
        files = [file_dir+'/'+x for x in os.listdir(file_dir)]
        files
        curvefile = [x for x in files if 'Curves' in x][0]
        curves_df = pd.read_csv(curvefile).set_index(['CurveID', 'Sheet'])
        pressurevalues = {}
        ##Get unique pressurevalues
        for pressurevalue in curves_df.index.unique(1):
            df = pd.read_csv([x for x in files if pressurevalue in x][0])
            df = df.set_index('Time').transpose()
            df['Sheet'] = pressurevalue
            df.index.names = ['CurveID']
            df.reset_index(inplace=True)
            df['CurveID'] = df['CurveID'].astype(int)
            df.set_index(['CurveID','Sheet'], inplace=True)
            df = df.join(curves_df)
            pressurevalues[pressurevalue] = df

        attempt_data = pd.concat([pressurevalues[df] for df in pressurevalues],axis=0).set_index(['LiquidClassName',
                                                                                    'Volume',
                                                                                    'StepType',
                                                                                    'Channel',
                                                                                    'Time',
                                                                                    'StepNumber',
                                                                                    'TadmMode',
                                                                                    'TadmError'],append=True)
        attempt_data['ParserID'] = attempt_id
        attempt_data.set_index(['ParserID'],append=True)    
        self.tadm_data = pd.concat([self.tadm_data, attempt_data])

    def closest_match(self, sample, main_process, max_time_offset=15, min_time_offset=0):
        """
        A function used to apply fuzzy logic to find tadms associated with a NeuMoDx Sample

        Parameters:
        ----------
        Sample (pd.DataFrame):  a slice of one row of data from NeuMoDx Raw Data
        main_process (str): Main Liquid Handling Process (LHPA, LHPB, LHPC) to use as time reference. 
        max_time_offset (int): An offset in seconds to apply to the maximum time bound applied to TADM search range.
        min_time_offset (int): An offset in seconds to apply to the minimum time bound applied to TADM search range.
        """
        def find_tadms(channel, repeat_offset=0):
            ##Get Test Guid of Sample
            test_guid = sample['Test Guid'].values[0]
            rep_number = sample['Replicate Number'].values[0]
            
            ##Convert Associated Start Date Time to be utc agnostic
            sample[main_process+' Start Date Time'] = sample[main_process+' Start Date Time'].apply(lambda x: x.replace(tzinfo=pytz.utc))
            
            ##Get Time of associated Sample
            time = sample[main_process+' Start Date Time'].astype('datetime64[ns]').values[0]
            
            ##Determine the processing group sample is associated with
            processGroupTimes = self.processGroups[main_process]
            processGroupTimes[main_process+" Start Date Time"] = processGroupTimes[main_process+" Start Date Time"].astype('datetime64[ns]')
            processGroupTimes['Reference Time'] = time
            processGroupTimes['Delta Time'] = abs(processGroupTimes[main_process+" Start Date Time"]-processGroupTimes['Reference Time'])
            
            ##Determine Minimum and Maximum Bounds for time allowed to search within.
            minimum_time_bound_index = processGroupTimes.loc[processGroupTimes['Delta Time']==processGroupTimes['Delta Time'].min(), main_process+" Start Date Time"].index.values[0]
            
            ##Apply a -1 run group offset in the case of a first time repeated sample.
            minimum_time_bound_index = minimum_time_bound_index - repeat_offset

            ##Get Value of minimum_time_bound
            minimum_time_bound = processGroupTimes.loc[minimum_time_bound_index, main_process+" Start Date Time"] - np.timedelta64(min_time_offset, 's')
            time = processGroupTimes.loc[minimum_time_bound_index, main_process+" Start Date Time"]

            if minimum_time_bound_index+1 < len(processGroupTimes):
                maximum_time_bound = processGroupTimes.loc[minimum_time_bound_index+1, main_process+" Start Date Time"] + np.timedelta64(max_time_offset, 's')
            else:
                maximum_time_bound = minimum_time_bound + np.timedelta64(5, 'm')
            

            
            ##Filter TADM Reference to only be for the Channel and Time Range allowed to search within
            tadm_reference_channel = self.tadm_data.reset_index(['Channel', 'LiquidClassName', 'StepType','Time'])
            tadm_reference_channel['Time'] = tadm_reference_channel['Time'].astype('datetime64[ns]')
            tadm_reference_channel = tadm_reference_channel[((tadm_reference_channel['Channel']==channel)&
                                                            (tadm_reference_channel['Time']>minimum_time_bound)&
                                                            (tadm_reference_channel['Time']<maximum_time_bound))]
            
            
            ##Determine Delta Time from Time observed for sample process
            tadm_reference_channel['Reference Time'] = time
            tadm_reference_channel['Delta Time'] = (tadm_reference_channel['Time']  - tadm_reference_channel['Reference Time']).astype('timedelta64[s]')

            ##Add Test Guid to TADM reference
            tadm_reference_channel['Test Guid'] = test_guid
            tadm_reference_channel['Replicate Number'] = rep_number
            tadm_reference_channel = tadm_reference_channel.reset_index()[['ParserID', 'CurveID', 'Test Guid', 'Replicate Number', 'Delta Time', 'LiquidClassName', 'StepType']].sort_values('Delta Time')
            
            ##Filter to make sure that we are only grabbing TADMs that we would expect based on process.
            if main_process == 'LHPB':
                tadm_reference_channel = tadm_reference_channel[((tadm_reference_channel['LiquidClassName'].str.contains(main_process))|(tadm_reference_channel['LiquidClassName'].str.contains('High')))]
            else:
                tadm_reference_channel = tadm_reference_channel[((tadm_reference_channel['LiquidClassName'].str.contains(main_process)))]

            ##Drop any duplicates that may have been found, keep lowest time delta.
            tadm_reference_channel.drop_duplicates(['LiquidClassName','StepType'],keep='first',inplace=True)
            
            
            return tadm_reference_channel


        channel = sample.loc[:, main_process+' ADP Position'].values[0]
        

        ##Determine which Channel to work with and if a sample is a aborted, or repeated sample.
        if pd.isnull(channel) or "nan" in channel:
            print("channel not found error")
            return

        elif "," in channel:
            channel_1 = pd.to_numeric(channel[-1])
            set1 = find_tadms(channel_1)
            channel_2 = pd.to_numeric(channel[0])
            set2 = find_tadms(channel_2, repeat_offset=1)
            return pd.concat([set1, set2],axis=0).drop_duplicates(['ParserID','CurveID'],keep='first')
        else:
            channel = pd.to_numeric(channel)
            set1 = find_tadms(channel)
            return set1

    def tadm_hunter(self):
        self.conversion_frame = pd.DataFrame(columns=['ParserID', 'CurveID', 'Test Guid', 'Replicate Number'])
        for id in self.raw_data_liquid_handle_processes.index.values:
            for process in ['LHPA', 'LHPB', 'LHPC']:
                self.conversion_frame = pd.concat([self.conversion_frame, self.closest_match(self.raw_data_liquid_handle_processes.loc[[id]], process)],axis=0)
                
        self.conversion_frame = self.conversion_frame[['Test Guid', 'Replicate Number', 'ParserID', 'CurveID']].set_index(['Test Guid','Replicate Number', 'ParserID', 'CurveID'])

    def tadm_merger(self, include_XPCR_info=0, include_HM_info=0, include_TS_info=0, include_SP_info=0, include_ConsLot_info=0, include_ChannelResult_info=0):
        raw_data_file_columns = ['Test Guid','Sample ID', 'Replicate Number','Overall Result','N500 Serial Number']
        
        if include_XPCR_info == 1:
            raw_data_file_columns = raw_data_file_columns + ['XPCR Module Serial','XPCR Module Index','Pcr Cartridge Lane']

        if include_HM_info == 1:
            raw_data_file_columns = raw_data_file_columns +['Heater Module Serial','Heater Module Index','Capture Plate Well']

        if include_TS_info == 1:
            raw_data_file_columns = raw_data_file_columns +['Test Strip NeuMoDx Carrier', 
                                                            'Test Strip NeuMoDx Carrier Position', 
                                                            'Test Strip NeuMoDx Well',
                                                            'Test Strip LDT Primer Probe Well',
                                                            'Test Strip LDT Primer Carrier',
                                                            'Test Strip LDT Primer Carrier Position',
                                                            'Test Strip LDT Master Mix Carrier',
                                                            'Test Strip LDT Master Mix Carrier Position',
                                                            'Test Strip LDT Master Mix Well']

        if include_SP_info == 1: 
            raw_data_file_columns = raw_data_file_columns + ['Sample Type', 
                                                            'Sample Specimen Type', 
                                                            'Test Specimen Type', 
                                                            'Specimen Tube Type', 
                                                            'Assay Name', 
                                                            'Result Code', 
                                                            'Status']
        if include_ConsLot_info == 1:
            raw_data_file_columns = raw_data_file_columns + ['LDT Test Strip Primer Probe Lot',
                                                            'LDT Test Strip Master Mix Lot',
                                                            'Pcr Cartridge Lot',
                                                            'Capture Plate Lot',
                                                            'Test Strip NeuMoDx Lot',
                                                            'Buffer Lot',
                                                            'Release Reagent Lot',
                                                            'Wash Reagent Lot']
        if include_ChannelResult_info == 1:
            for channel in self.channels:
                raw_data_file_columns = raw_data_file_columns+[channel + " " + x for x in ['Localized Result', 'Ct', 'End Point Fluorescence', 'EPR', 'Max Peak Height', 'Baseline Slope', 'Baseline Y Intercept', 'Blank Check - 1st 3 Reads'] if channel + " " + x in self.raw_data.columns]

        raw_data_index = self.raw_data.loc[:, raw_data_file_columns].set_index(['Test Guid', 'Replicate Number'])
        raw_data_index = raw_data_index.join(self.conversion_frame)
        merged_data = raw_data_index.join(self.tadm_data.reset_index().set_index(['ParserID', 'CurveID']))
        return merged_data.reset_index().set_index([x for x in self.tadm_data.index.names] + raw_data_file_columns)

window_width = 1000
window_height = 400
windowsize = str(window_width)+"x"+str(window_height)
root = Tk()
root.title("TADM Matcher v0.1")
root.geometry(windowsize)
my_gui = UI(root)
root.mainloop()
