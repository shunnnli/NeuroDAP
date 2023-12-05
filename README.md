# NeuroDAP
 NeuroDAP (Neuroscience Data Analysis Pipeline) for recording, synchronization, preprocessing, and analysis of neuropixel, photometry, camera, and behavior recordings

 ## General structure

 NeuroDAP are built with **four major stages** in mind with variable levels of customization. 
 1. **Recording phase** (customize based on rig configuration, such as defining sync pulse is stored in channel 1, licking in channel 2, GCaMP in channel 3 etc)
 2. **Preprocessing phase** (synchronize between acquisition systems through common sync pulse, assign a common timestamp to every sample of each recording system)
 3. **Session analysis** (extract trial table for each session; align signals and perform basic analysis of these aligned signals; plot session summary)
 4. **Experiment analysis** (pooled all sessions across all animals; perform data analysis)

 Please refer to the wiki section for more detial documentation.

 ## Introduction of sample data set

The sample data set is recorded by Shun Li in 2023. It contains two parts: 1）a sample recording session for showing preprocessing/session analysis phase 2）a sample `animals` struct for showing experiment analysis phase. The `animals` struct contains 4 animals, with 1 animals with off-target expression ('SL137'). dLight signals in NAc, pupil/Eye area, and lick are simultaneously recorded for all sessions.

There are 5 major phases:
1. Random: water, airpuff, EP stim, and tone (75dB) are delivered randomly
2. Reward1/2: where EP stim and tone are paired with water
3. Punish1/2: where EP stim and tone are paired with airpuff
4. Timeline: Random (2 sessions) -> Reward1 (3 sessions) -> Punish1 (3 sessions) -> Reward2 (3 sessions) -> 1 week rest -> Punish2 (3 sessions but 3 animals)

 ***

