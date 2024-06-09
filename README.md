# Running the code

## Setup
1. Download the SPDZ distribution (https://github.com/data61/MP-SPDZ/releases) and unpack it
2. In /Programs/Source clone our git repo (https://github.com/Amaimaya/MPC-EFC)
3. Navigate to the root folder of the repo (~/mp-spdz-0.3.8/Programs/Source/EFC-MP-SPDZ) 
4. Run the **./run.sh** script that does the following:
   1. Pull from git
   2. Navigate to place from which the programs are compiled, compile and run it 
   3. Navigate back to the EFC-MP-SPDZ folder

Now anytime you want to run the **main.mpc** code that is in the root folder of the repository, just run **./run.sh**

## Running experiments

To run the experiments, do the following:
1. Navigate to the folder containing files for a desired dataset
2. Run **./run.sh** file that is in that folder