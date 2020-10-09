# Easy Osc3++ - v3.1

The core of this simplified version of osc3++ is given by the classes in the ```event``` folder.
The main code in ```app/fitter.cpp``` uses derived classes inheriting from ```Sample.h``` to load systematics and build distributions, whereas the ChiSquared class computes and and minimises the chi-squared.
Support to beam sample and atmospheric sample is available.
The oscillation space is scanned over and it is defined and managed by the ```physics/ParameterSpace``` class.
Other important classes are ```physics/Oscillation``` to deal with oscillation physics and ```physics/Atmosphere```
to compute the Earth density profile seen by an atmospheric neutrino.

The space is specified by combinations of oscillation parameters defined inside ```cards/oscillation.card```.
For each point of the oscillation parameter space, the observables are created and oscillated with those oscillation parameters.
The default set is build around the T2K *Asimov A* point.
The true point is defined by the *Point* variable; if not defined, the default point is used and it is computed as explained in the card.
The cards ```cards/beam_sample.card``` and ```cards/atmo_sample.card``` defines the input files to build beam and amtospheric samples.
The ```cards/fit_options.card``` regulates fitting parameters.

A set of script is used to facilitate the execution of the code and extraction of results.
In executing the following scripts, the structure of the folder is important.
We will refer to where the systematic files, outputs of the fit, and plots files will be as the ```root``` folder.

## Differences with previous version

* Compilation across te cluster is now available
* Since version v2.0 the energy scaling error is handled more exactly because it is implemented as an anlaytic function.  Hence, the mathematical exact Jacobian and Hessians are used during minimisation.
* The Oscillator and ChiSquared classes have been remodelled to improve optimization and maximise use of vector instructions.
* The chi-squared penalisation term is automatically added and it can be specified in the ```cards/oscillation.card```.
* The atmospheric sample is now supported.
* The CardDealer class, to import parameters dynamically from text files, has been rewritten as a template class to support as many data types as possible.

## Getting started

The information here are aimed at a quick installation of the code and how to lunch the fitter rapidly.
The user is advised to look at the documentation first.

## Requirements

The requirements for the code to be compiled and run are
* *make*;
* *gcc-c++* version 4.8 or higher (for C++11);
* *ROOT* version 5.34/38 or higher;
* *Eigen* version 3.3 or higher;
* *HTCondor* or *SLURM* as workload managers, for distributed computing.

[ROOT][https://root.cern.ch/] should be installed and properly linked. To test if true, simply run
```
root-config --cflags --glibs
```

You must download and extract the [Eigen 3 library](https://eigen.tuxfamily.org/dox/index.html) and set an environmant variable ```EIGEN``` to point to the Eigen top folder, or edit the Makefile and manually change the ```EIGENINC``` variable, or copy the folder ```Eigen/``` from the Eigen installation folder to the local ```include/``` folder.


## Build

Running
```
make
```
will create the executables in the ```bin/``` folder.
If only one binary is needed, then this can be compiled alone with
```
make APP = <main>
```
where main any main file under the local app folder. The .cpp extensions must be omitted. 

If using a distributed computing  cluster, its nodes could have all different architectures.
Specific compilation on each node can be achieved by using
```
./cross-compile_c . sh # if HTCondor is the manager
./cross-compile_s . sh # if Slurm is the manager
```
It is not a cross-compilation properly speaking, as both scripts first determine the nodes on the
cluster, then ssh into each node and only if a new architecture is found a new, optimized binary
is compiled for that particular architecture. This workaround works as long as access to the file
system is shared on the cluster and that each node can be accessed via ssh by the user.


## Folder structure

After building the executables, run
```
./setup.sh
```
to build the folder structure required by the framework.
The script also downloads the input files from http://hep.lancs.ac.uk/~tdealtry/oa/ and https://pprc.qmul.ac.uk/~tboschi/HK/atmo.

You can specify a specific path with the ```-p prefix``` option, as specified in the usage. The default value is the folder ```errorstudy/``` in the current working directory 

The subfolders ```errorstudy/reconstruction_beam``` and ```errorstudy/reconstruction_atmo``` are created to contain reconsturction files to build the data samples.



## Prepare the beam systematics

This step is not necessary if you have run ```setup.sh``` before.

Creates correlation matrix of systematic parameters by combining matrices found in files of matrixN.root.
The script expects to find the systematics folder under root with the spline files to be processed (renaming of files and histograms).
It uses ```app/purifysystematics.cpp``` and ```app/addmatrix.cpp```.
```
./prepare_systematics.sh -r errorstudy/root matrix1.root [matrix2.root ...]
```

## Running the fitter

The fitter is the most important tool of the framework, whereas the other executables are
mostly validation or printing tools. Before running it, the user has to decide which sample to
fit:
* the beam sample;
* the atmospheric sample;
* the beam and atmospheric samples together.
Next, the true and fitted neutrino mass hirerachy should be chosen and whether or not fit the
systematic errors. The systematic files should be provided in any case.
A new folder in the analysis directory should be created whenever the samples or the sys-
tematic models are changed. For example,
```
mkdir errorstudy/first_run
```
Then, a sub folder containing the systematic files should be created (or copied) in this new
directory. The nominal T2K 2018 systematic model is found under ```errorstudy/0```, therefore
```
cp -r errorstudy /0/systematics errorstudy/first_run
```
Finally, the fitter can be launched on the distributed computing system with the trisens_c.sh
(HTCondor) or the trisens_s.sh (Slurm) utilities. For example,
```
# with HTCondor
./trisens_c.sh -r errorstudy/first_run -d comb -1 NH -2 NH -N 500
# with Slurm
./trisens_s.sh -r errorstudy/first_run -d comb -1 NH -2 NH -N 500
```
will launch 500 jobs (```-N 500```) on the cluster, fitting both the beam and atmospheric samples (```-d comb```) with true and fitted normal mass hierarchies (```-1 NH -2 NH```).
Other important options are
* ```-s``` to do a statistics only fit, i.e. no systematics;
* ```-f <scan_type>``` to perform a scan fit, i.e. chaning the true point at each iteration as it is done for CPV sensitivity studies;
* ```-v <verbosity>``` to change the verbosity of the log files.
The full list of options can be shown with
```
./trisens_c.sh -h # for HTCondor
./trisens_s.sh -h # for Slurm
```
