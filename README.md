# NeuroDAP
 NeuroDAP (Neuroscience Data Analysis Pipeline) for recording, synchronization, preprocessing, and analysis of neuropixel, photometry, camera, and behavior recordings

 ## General structure

 NeuroDAP are built with **four major stages** in mind with variable levels of customization. 
 1. Recording phase (customize based on rig configuration, such as defining sync pulse is stored in channel 1, licking in channel 2, GCaMP in channel 3 etc)
 2. Preprocessing phase (synchronize between acquisition systems through common sync pulse, assign a common timestamp to every sample of each recording system)
 3. Session analysis (extract trial table for each session; align signals and perform basic analysis of these aligned signals; plot session summary)
 4. Experiment analysis (pooled all sessions across all animals; perform data analysis)

 Below, key functions and places for customization will be described. Detail implementation please refer to specific code.

 ## Recording phase
 1. Key functions: depend on specific rig. For data acquisition using the labjack system, please see run_labjack.mat and below for details
 2. Key files: session folder that contains respective data format for each acquisition system
 
 ### Key points: 
 - a common sync pulse of random inter-pulse interval should be sent to all acquisition system. Otherwise synchronization during preprocessing will not be possible
 - If labjack-based recording is used, there're two important scripts/functions:
    - run_labjack.mat: this script defines the labjack settings (see below) and runs the acquisition of labjack during recording. It will automatically save inside the session folder (ie sessionName/Photometry)
    - concatLabjack_setupName.mat: this function should be customized/edited based on individual rig setup. This defines the content of each channels and fills in empty labjack fields for analysis later.
 - In run_labjack.mat
    - There should be a labjack struct that contains following information. These information will be used during preprocessing phase for analysis steps like demodulation/detrending/z-score.
        - labjack.name: name of each recorded channel (eg {'NAc','LHb','PMT'} or {'NAc-green','NAc-Isosbestic'})
        - labjack.record: whether the channel is recorded/use for analysis (eg [1,0,0] if I don't have LHb recording)
        - labjack.mod: whether the channel is amplitude modulated (eg [1,1,0])
        - labjack.modFreq: what is the frequency of amplitude modulation for each channel (eg [200,250,nan])
        - labjack.LEDpower1: average power output of LED1 (calibrate every day)
        - labjack.LEDpower2: average power output of LED2 (calibrate every day)
        - labjack.LEDpowerMin1: minimal power output of LED1 (~5uW for me)
        - labjack.LEDpowerMin2: minimal power output of LED2 (~5uW for me)
    - User can define at which frequency the signal is modulated at. The script will takes into account to lowest power that user provided (to avoid absolutely no signal) and autocalculate the maximal amplitude of modulation to prevent clipping. The average power will be the same as labjack.LEDpower
 - In concatLabjack_setupName.mat
    - As this function is heavily dependent on rig setup, please refer to the corresponding code

 ## Preprocessing phase
 1. Key functions: loadSessions(sessionpath,options)
 2. Key files: data_sessionName.mat, timeSeries_sessionName.mat, sync_sessionName.mat, behavior_sessionName.mat

 ### Key points: 
 - 

 ## Session analysis
 1. Key functions: analyzeSessions_XXX(sessionpath,options)
 2. Key files: data_sessionName.mat, timeSeries_sessionName.mat, sync_sessionName.mat, behavior_sessionName.mat
 
 ### Key points: 
 - 

 ## Experiment analysis
 1. Key functions: analyzeExperiments_XXX(sessionpath,options)
 2. Key files: analysis_sessionName.mat; data_sessionName.mat, timeSeries_sessionName.mat, sync_sessionName.mat, behavior_sessionName.mat
 
 ### Key points: 
 - 