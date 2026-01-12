import os
import sys
from pathlib import Path

o2_folder = Path('/n/data1/hms/neurobio/sabatini/Shun/DataJoint/Inbox/')
# root_folder = Path('/Volumes/MICROSCOPE/Shun/Project valence/Recordings/') # change the last folder name if necessary
root_folder = Path('/n/files/Neurobio/MICROSCOPE/Shun/Project\ valence/Recordings/')

experiment = '202310-Photometry-paAIP2'
session = '20231108-SL140-R12_g0' #sys.argv[1]
session_folder = os.path.join(root_folder,experiment,session)

# Get animal name
animal = session.split('-')[1]
o2_session_folder = os.path.join(o2_folder, animal, session)
o2_photometry_folder = os.path.join(o2_session_folder,'Photometry')
o2_behavior_folder = os.path.join(o2_session_folder,'Behavior')
os.system("mkdir -p " + o2_photometry_folder)
os.system("mkdir -p " + o2_behavior_folder)

# Move session to O2
print("Ongoing: moving session " + session + " to O2...")
print("     Server path: " + session_folder)
print("     O2 path: " + o2_session_folder)

for mat in list(Path(session_folder).glob("data*.mat")):
    # temp = f'"{str(mat)}"'
    # file_path = f"'{temp}'"
    file_path = mat
    print(mat)
    # print("rsync -r --mkpath " + file_path + " " + o2_photometry_folder)
    # os.system("mkdir " + o2_photometry_folder)
    print("Finished: rsync data.mat file")
    os.system("cp " + file_path + " " + o2_photometry_folder)

for mat in list(Path(session_folder).glob("timeseries*.mat")):
    # temp = f'"{str(mat)}"'
    # file_path = f"'{temp}'"
    file_path = mat
    # print("rsync -r " + file_path + " " + o2_photometry_folder)
    print("Finished: rsync timeseries.mat file")
    os.system("cp " + file_path + " " + o2_photometry_folder)

for mat in list(Path(session_folder).glob("*.parquet")):
    # temp = f'"{str(mat)}"'
    # file_path = f"'{temp}'"
    file_path = mat
    print("Finished: rsync *.parquet file")
    # print("rsync -r --mkpath " + file_path + " " + o2_behavior_folder)
    os.system("cp " + file_path + " " + o2_behavior_folder)