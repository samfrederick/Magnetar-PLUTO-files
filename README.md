# Magnetar-PLUTO-files
This repository provides files for modeling magnetar structure in PLUTO, an astrophysical dynamics software package.   

This work is in active development during the 2018-2019 academic year for the completion of an honors thesis in physics from Davidson College. 

## Dependencies 
In order to utilize these files, the following packages must be installed:
 * PLUTO v4.30
 >* Python
 >* C Compiler 
 >* GNU Make
* Visualiztion software utility: Numerous options are compatable with PLUTO including Gnuplot, IDL, pyPLUTO, Mathematica, Para View, or VisIt. 

The code provided in this repository is largely written in C, with the exception of a few python scripts written in 2.7 under the folder "Python Files". 

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
user@computer:~$ cd $PLUTO_DIR/Test_Problems
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

## Running Magnetar Simulations 
The following files provide the core of each PLUTO simulation:
* *definitions.h*
* *init.c*
* *makefile*
* *pluto.ini*

#### These files are included in the folder *"Magnetar Main"*

*Additional in-depth information regarding the function and manipulation of each of these scripts is available via the PLUTO User Manual*

### *definitions.h* 
A header file specifying parameters for a given simulation. Most of these values can be specified via *setup.py*, however, the user may specify additional simulation constraints which are detailed further in the PLUTO User Manual. 
* The *definitions.h* file particular to these magnetar simulations uses PLUTO's magnetohydrodyanmic (MHD) module and specifies a 3-dimensional spherical computational domain.  

### *init.c*
The init.c file collects most of the user-supplied functions useful for problem configuration. 
* The *init.c* file particular to these simulations provide initialization for simulation variables, such as velocity, pressure, density, magnetic field components, etc. 
* Boundary Conditions are specified in the UserDefBoundary() function. 
* Gravitational Potential is assigned via the BodyForcePotential() function. 

### *Makefile*
In order to compile modifications to the code, the user must run the following command
```console
user@computer:~/PLUTO$ make
```
### *pluto.ini* 
An initialization file which sets simulation attributes. Upon execution of the PLUTO Code via the command
```console
user@computer:~/PLUTO$ ./pluto
```
the code reads this initialization file to determine the specified simulation parameters. This begins the process of computation.
* The *pluto.ini* file particular to these simulations specifies a spherical domain: radially (0<r<=2.0) in normalized units of stellar radius, (0 < theta < pi), and ( 0 < phi < 2pi) where we use the physics convention specifying phi as the azimuthal angle. 
