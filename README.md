# Citrine
A Code for the simulate the dynamic characterics of cable


Author's Information:

github Personal homepage: **https://github.com/Lin912**

ORCid (0000-0003-3820-3199) **[link](https://orcid.org/)**


## NewRelease Version(5.2)

**content:
For linux Server we release the Citrine 5.2. This code used for calculate the dynamic behaivor of cable.
In Citrine 5.2, the Porocess Folder contains the PreProcess and PostProcess that include the Jacabian Matrix 
calculation procedure, Data Processing and Plotting part. Then, the folder ForceBoundary and VelocityBoundary:

(ForceBoundary)  <Top Vel & Bottom Forces>

(VelocityBoundary)  <Top Vel & Bottom Vel>





## 1.0  Theory





## 2.0   Usage process
The code is written by the C++, before download the code you need:
'''
1. Prepare a compilable c++ environment. With download the GNU complier, like gcc, g++ etl;
2. Download and install the Eigen library;
3. Install Cmake and compile the source files according to CmakeLists in the Citrine5 folder.
'''
Running the Citrine5.exe and view the calculations result in the csv folder. The output.csv is the result file.

Then viewing the graphical results with using the code in Matlab folder, the PostProcess folder includes the Force and Profile result.

Note: 
1. Change the file "setting.json" to specify the path to the compiler.
2. Specify the CMake source path and the CMake executable path.


## 3.0 Showcase:
We did a set of related pendant calculations based on the above theory and code, and the specific distribution state and mechanical properties of the flexible body are shown below:

<div style="display: inline-block; text-align: center;">
  <img src="https://github.com/Lin912/Citrine5/blob/main/ResultShow/Experiment1.jpg" alt="Experment procedure 1" width="500"/>
  <img src="https://github.com/Lin912/Citrine5/blob/main/ResultShow/Experiment2.jpg" alt="Experment procedure 2" width="500"/>
</div>

## Copyright && Contact Me
Email:  Z0802816@gmail.com
