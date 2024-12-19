# Overview
- A privacy-preserving solution for the Energy-Based Flow Classifier (EFC), a machine learning algorithm for detecting network intrusions (https://github.com/EnergyBasedFlowClassifier/EFC-package). 
- Privacy is ensured using secure multi-party computation (MPC) implemented with the MP-SPDZ framework, a modified version of SPDZ-2 (https://github.com/data61/MP-SPDZ).
- Optimization techniques have been used to balance privacy and performance, resulting in a significantly more efficient algorithm compared to its naive version.

# Running the code

## Setup
1. Download the SPDZ distribution (https://github.com/data61/MP-SPDZ/releases) and unpack it
2. In /Programs/Source clone our git repo (https://github.com/Amaimaya/MPC-EFC)
3. Navigate to the root folder of the repo (~/mp-spdz-0.3.8/Programs/Source/MPC-EFC) 
3. Navigate to the folder containing the desidered version of the EFC algorithm (plaintext, base, optimized or optimized-rabbit) or the Rabbit protocol benchmarks (rabbit)
4. Run **./run.sh** or **./run_online.sh** (for the online phase only, follow the fake offline generation instructions in the MP-SPDZ repository) file that is in that folder, which does the following:
   1. Pull from git
   2. Navigate to place from which the programs are compiled, compile and run it 
   3. Navigate back to the MPC-EFC folder

Now anytime you want to run the **x.mpc** code that is in the desidered folder of the repository, just run **./run.sh** or **./run_online.sh**
