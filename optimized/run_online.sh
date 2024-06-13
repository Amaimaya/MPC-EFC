git pull
/bin/python3 /home/mpspdz/mp-spdz-0.3.8/Programs/Source/MPC-EFC/optimized/test_selectors.py
cd ../../../..
./compile.py MPC-EFC/optimized/optimized.mpc -F 31
Scripts/mascot.sh -F optimized
cd Programs/Source/MPC-EFC/optimized/