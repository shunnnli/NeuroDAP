# NeuroDAP
 NeuroDAP (Neuroscience Data Analysis Pipeline) for recording, synchronization, preprocessing, and analysis of neuropixel, photometry, camera, and behavior recordings

 Test

 ## General structure

 NeuroDAP are built with *four major stages* in mind with variable levels of customization. 
 1. Recording phase (customize based on rig configuration, such as defining sync pulse is stored in channel 1, licking in channel 2, GCaMP in channel 3 etc)
 2. Preprocessing phase (synchronize between all recorded system through common sync pulse, assign a common timestamp to every sample of each recording system)
 3. Session analysis (extract trial table for each session; align signals and perform basic analysis of these aligned signals; plot session summary)
 4. Experiment analysis (pooled all sessions across all animals; perform data analysis)

 Below, key functions and places for customization will be described. Detail implementation please refer to specific code.

 ## Recording phase

 ## Preprocessing phase

 ## Session analysis

 ## Experiment analysis