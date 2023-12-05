# NeuroDAP
 NeuroDAP (Neuroscience Data Analysis Pipeline) for recording, synchronization, preprocessing, and analysis of neuropixel, photometry, camera, and behavior recordings

 ## General structure

 NeuroDAP are built with **four major stages** in mind with variable levels of customization. 
 1. **Recording phase** (customize based on rig configuration, such as defining sync pulse is stored in channel 1, licking in channel 2, GCaMP in channel 3 etc)
 2. **Preprocessing phase** (synchronize between acquisition systems through common sync pulse, assign a common timestamp to every sample of each recording system)
 3. **Session analysis** (extract trial table for each session; align signals and perform basic analysis of these aligned signals; plot session summary)
 4. **Experiment analysis** (pooled all sessions across all animals; perform data analysis)

 Please refer to the wiki section for more detial documentation.

 ## Intro of sample data set

The sample data set is recorded by Shun Li in 2023. It contains 4 animals, with 1 animals with off-target expression ('SL137'). dLight signals in NAc, pupil/Eye area, and lick are simultaneously recorded for all sessions.

There are 5 major phases:
1. Random: water, airpuff, EP stim, and tone (75dB) are delivered randomly
2. Reward1/2: where EP stim and tone are paired with water
3. Punish1/2: where EP stim and tone are paired with airpuff
4. Timeline: Random (2 sessions) -> Reward1 (3 sessions) -> Punish1 (3 sessions) -> Reward2 (3 sessions) -> 1 week rest -> Punish2 (3 sessions but 3 animals)


 ## Recording phase
 1. Key functions: depend on specific rig. For data acquisition using the labjack system, please see ```run_labjack.mat``` and below for details
 2. Key files: session folder that contains respective data format for each acquisition system
 
 ### Key points: 
 - a common sync pulse of random inter-pulse interval should be sent to all acquisition system. Otherwise synchronization during preprocessing will not be possible
 - If labjack-based recording is used, there're two important scripts/functions:
    - ```run_labjack.mat```: this script defines the labjack settings (see below) and runs the acquisition of labjack during recording. It will automatically save inside the session folder (ie sessionName/Photometry)
    - ```concatLabjack_setupName.mat```: this function should be customized/edited based on individual rig setup. This defines the content of each channels and fills in empty labjack fields for analysis later.


 ## Preprocessing phase
 1. Key functions: ```loadSessions(sessionpath,options)```
 2. Key files: ```data_sessionName.mat, timeSeries_sessionName.mat, sync_sessionName.mat, behavior_sessionName.mat```


 ### Key points: 
 - Use ```loadSessions()``` to perform synchronization and necessary signal preprocessing.
 - ```loadSessions()``` contains following steps
    1. **Detects how many acquisition system exists within the session folder by searching for following files.**
        - NIDAQ and Neuropixel: .ap.bin or .imec.bin
        - Labjack photometry: /Photometry folder
        - Camera: cam1_XXX.avi and times_cam1_XXX.csv
    2. **Run signal preprocessing for each individual acquisition system**
        - For NIDAQ:
            - digital: assigns each digital channel to corresponding events (packaged into function later)
            - analog: threshold and extracts rising edge if licking is recorded, or perform downsample (to ```options.downsampleFs```) if photometry/movement traces are recorded
        - For labjack:
            - run ```concatLabjack_setupName.mat```. Within it, assign digital/analog channels
            - For analog&recorded channels (stored in rows of ```labjack.raw```), perform demodulation (```demodulatePhotometry()```) if corresponding ```labjack.mod(i) == true```, or downsample (```downsamplePhotometry()```) otherwise. Save corresponding result in ```timeSeries```.
            - The ```timeSeries``` struct will contain following fields and will be compatible with sabatini lab datajoint pipeline:
                - ```timeSeries.name```: name of each signal, same as ```labjack.name```
                - ```timeSeries.finalFs```: final sample frequency (should be equal to ```options.downsampleFs```)
                - ```timeSeries.system```: acquisition system (can be 'NI','LJ', or 'Cam')
                - ```timeSeries.demux_freq```: 
                - ```timeSeries.behavior_offset```: 
                - ```timeSeries.system```: 
                - ```timeSeries.system```: 
        - For camera:
            - read times_cam1_XXX.csv file. Extract sync pulse and camera readings from Bonsai
            - If there is Bonsai-analyzed data like pupil/eye area, perform downsample (```downsamplePhotometry()```) and optional detrending based on user input
    3. **Synchronization**
        - Calculate inter-pulse interval of sync pulse recorded in each acquisition system
        - Perform cross-correlation to find the first common sync pulse
        - Assign each sample of each acquisition system with a common timestamp. This timestamp will be used to cross-reference between different acquisition system.
            - For example: water is delivered at NIDAQ sample 10003240, which corresponding to common time 10.1s. To find the photometry traces in Labjack aligned to this event, I can find the closest timestamp in labjack system (eg. 10.1003s) and plot the traces
        - Calculate time offset comparing to behavior (used in sabatini lab datajoint pipeline)
        - Save everything in ```params.sync``` fields and save ```params``` struct in ```sync_sessionName.mat```
- ```params``` struct in ```sync_sessionName.mat``` will be the struct that stores ALL the key params related to synchronization, acquisition, and behavior type of the session. Following session analysis phase will heavily utilize and refer to this structure.


 ## Session analysis
 1. Key functions: ```analyzeSessions_XXX(sessionpath,options)```
 2. Key files: ```analysis_sessionName.mat; data_sessionName.mat, timeSeries_sessionName.mat, sync_sessionName.mat, behavior_sessionName.mat```
 
 ### Key points: 
 - The general workflow of ```analyzeSessions_XXX(sessionpath,options)``` is as follows:
    - Load all session-related data ```data_sessionName.mat, timeSeries_sessionName.mat, sync_sessionName.mat, behavior_sessionName.mat```
    - Filled in key information of the session if needed. (See ```analyzeSessions_OptoPair(sessionpath,options)``` for examples)
    - (Optional) further signal preprocessing
    - Generate trial table (```getTrialTable_XXX(sessionpath,options)```)
        - There should be a seperate function for each person/behavior tasks
        - In essence, it is a place to store all relevant behavior events for future alignment and analysis.
        - Also, create ```eventTable``` and ```blockTable``` for integrating with sabatini lab datajoint pipeline
    - Extract all important events worthy of analysis (```analysisEvents```). Also find their corresponding trial number (```eventTrialNum```) and define their corresponding name (```analysisLabels```). For example:
        - ```analysisEvents = {waterIdx, airpuffIdx, stimIdx, toneIdx}``` (where ```waterIdx=[213,3242,534234,436534]```)
        - ```analysisLabels = {'Water', 'Airpuff', 'Stim', 'Tone'}```
        - ```eventTrialNum = {waterTrials, airpuffTrials, stimTrials, toneTrials}``` (where ```waterTrials=[1,4,33,43]```)
    - Run ```analyzeTraces(analysisEvents,analysisLabels,eventTrialNum=eventTrialNum)```
        - This function will automatically align every signal in ```timeSeries``` struct to all input events and stored them in ```analysis_sessionName.mat```
        - The reason of creating a seprate ```analysis_sessionName.mat``` is as follows:
            - In following experiment analysis phase, loading all the raw data and the ```params``` struct for referencing across acquisition system is extremely slow and inefficient, therefore, ```analysis_sessionName.mat``` will store all the aligned events and can be efficiently concatenated in the next phase.
        - Notable options of ```analyzeTraces()```
            - ```timeRange```: seconds before and after events, default is 15s (```[-15,15]```)
            - ```stageTime```: define subtrial stages for further analysis. For example, you can define 2s before the event as baseline, 0.5s after the event as Cue period, and everything after that is Response period (```[-2,0;0,0.5;0.5,15]```). ```analyzeTraces()``` will calculate the average, max, and min for each stages and stored in field ```stageAvg.data```,```stageMax.data```,```stageMin.data``` respectively. Moreover, statistical analysis will also conducted within each stage to calculate the best fit line for the current stage across trials (saved in ```stageAvg.fit```), and bootstrap analysis will be done to calculate the p-value of calculated slope and intercept.
    - Design session summary plots and loop through all signals
        - This is highly dependent on specific experiments and analysis needs. See ```analyzeSessions_OptoPair(sessionpath,options)``` for examples. 

 ## Experiment analysis
 1. Key functions: ```analyzeExperiments_XXX(sessionpath,options)```
 2. Key files: ```analysis_sessionName.mat; data_sessionName.mat, timeSeries_sessionName.mat, sync_sessionName.mat, behavior_sessionName.mat```
 
 ### Key points: 
 - To start, copy this matlab file and replace template with specific experiments. In theory, there will be a specific analyzeExperiments file for each individual experiments, as the specific needs for analysis varies between different experiments.
 - The analysis pipeline are as follows:
    1. Select whether to load a previously ```animals``` struct (described below) or select individual session to combine.
    2. After selecting ALL SESSIONS from an experiments, the pipeline will automatically concatenate ```analysis_sessionName.mat``` for each recording sessions. Rename properties as needed in order to facilitate further analysis.
    3. Run ```getAnimalStruct(summary)``` function to recreate ```animals``` struct from ```summary``` struct. This combines all sessions from the same animals together while cutoffs between individual sessions are also recorded.
    4. Save ```animals``` struct if needed. Note: saving ```summary``` struct will take extremely long (>5hrs) so while saving ```animals``` struct is much shorter  (~2min). ```animals``` struct should contain information that satisfies MOST plotting requirements so saving ```summary``` struct is not needed.
    5. Data analysis and plotting. This part is designed to vary across experiments. Thus, following codes are just for demonstration of essential functions.

