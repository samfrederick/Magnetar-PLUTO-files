# Magnetar-PLUTO-files
This repository provides files for modeling magnetar structure in PLUTO, an astrophysical dynamics software package.   

## Dependencies 
In order to utilize these files, the following packages must be installed:
 * PLUTO v4.30
 >* Python
 >* C Compiler 
 >* GNU Make
* Visualiztion software utility: Numerous options are compatable with PLUTO including Gnuplot, IDL, pyPLUTO, Mathematica, Para View, or VisIt. 

The code is largely written in C, with the exception of a few python scripts written in 2.7 under the folder "Python Files". 

## Installation Procedure

### 1. PLUTO

The PLUTO Code must be downloaded from the following link: http://plutocode.ph.unito.it/download.html. 

Download and extract the file to your home directory. It is recommended to create the enivronment variable "PLUTO_DIR" so that the user can quickly access the PLUTO directory. In bash, this is accomplished via the command:
```console
user@computer:~$ export PLUTO_DIR=/home/user/PLUTO
```

The python script "setup.py" included with PLUTO can be used to quickly configure simulation parameters. Settings are modified by arrow keys and selected using "Enter". 

To verify proper installation of PLUTO, it is recommended that the user run one of the various test problems included with the PLUTO base by navigating to:
```console
user@computer:~$ $PLUTO_DIR/Test_Problems
```
More information including step-by-step guides to test problems is included in the PLUTO User Manual, available via:
http://plutocode.ph.unito.it/files/userguide.pdf (See pages 6-9 for test problems)

### 2. VisIt
Although a varity of visualization software packages are compatable with PLUTO, this work has implemented VisIt v2.13.3, for which we provide detailed installation instructions. 

Visit may be downloaded from the following link: https://wci.llnl.gov/simulation/computer-codes/visit/executables

Installation for Windows and Mac is fairly straightforward, as .exe and .dmg installer programs are provided for the respective operating system. 

#### Installing VisIt on Unix Systems 
* Run the VisIt install script. The script for v2.13.3 is available at the following link: http://portal.nersc.gov/project/visit/releases/2.13.3/visit-install2_13_3 

* Add a bin directory after the location of the installed directory, i.e ( /user/local/visit/bin ). This directory will contain the script responsible for launching VisIt. This may be accomplished by adding it to the users .cshrc folder:
```console
user@computer:~$ echo "set path = ($path /usr/local/visit/bin)" >> .cshrc
```
*More details regarding VisIt installation is available at the following link:* http://portal.nersc.gov/project/visit/releases/2.13.3/INSTALL_NOTES

*It is recommended for systems with NVidia graphics cards to update to the latest drivers. Failure to do so may result in frequent program crashes.*

The four separate files: definitions.h, makefile, pluto.ini, and init.c are required for compiling the PLUTO code. 
The code is compiled via the command: 
```console 
user@computer:~$ make 
```
Subsequently, the user must type the command 
```console
user@computer:~$ ./pluto
```

