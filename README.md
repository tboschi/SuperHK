# SuperHK - an easier approach to Osc3++

The ```osc3++``` framework is a bundle of computational and plotting utilities helpful to asses the sensitivity of a future experiment, like Hyper-Kamiokande, to neutrino oscillation parameters.

The workflow of the framework can be structured as follows
* a combination of oscillation parameters is chosen to be the "true" (observed) one;
* beam and atmospheric predictions at the detector are built twice, using the "true" combination of parameters (observed events) and using a combination under examination (expected events);
* a chi2 between observed and expected events is created as a function of the systematic parameters;
* the minimum value of chi2 is found with respect to the systematic parameters;
* the process is repeated at different combinations of oscillation parameters.

Please refer to the documentation for a full description of the software.


### Current version

The latest release can be found in [releases](https://github.com/tboschi/SuperHK/releases).

Differences with previous version:
* Sources and script have been moved to the `src` directory. This is in view of moving to an API+library style.
* Sparse matrices are used now to compure Jacobian and Hessian improving speed by a factor of 5.

Differences with v2.0:
* Atmospheric sample is supported.
* Cross-compilation across the cluster is now available; see the documentation on how to use it.
* The energy scaling error is handled more exactly as it is implemented as an anlaytic function.
* The ```Oscillator``` and ```ChiSquared``` classes have been remodelled to improve optimization and maximise use of vector instructions.
* The chi2 penalisation term is automatically added.
* The ```CardDealer``` class has been rewritten as a template class.


## Getting started

The information here is aimed at a quick installation of the code and how to lunch the fitter rapidly.
The user is advised to look at the documentation first.

## Requirements

The requirements for the code to be compiled and run are
* **make**;
* **gcc-c++** version 4.8 or higher (for C++11);
* **ROOT** version 5.34/38 or higher;
* **Eigen** version 3.3 or higher;
* for running on the cluster **HTCondor**, **SLURM**, **SGE**, or **PBS** as workload managers for distributed computing; if you use another manager you should change the script accordingly

[ROOT](https://root.cern.ch/) should be installed and properly linked. To test if true, simply run
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
make APP=<main>
```
where main any main file under the local app folder. The .cpp extensions must be omitted. 
The `make` command also copies some bash scripts from the source directory to the bin directory.
It is strongly advised to apply any edits to the source scripts first and then run `make` again.

By default, SMIX and AVX are turned on for Eigen optimization.
To deactivate them compile like this
```
make ARCH=
```

If using a distributed computing cluster, its nodes could have all different architectures.
Extremely optimized compilation on each node can be achieved by using
```
bin/cross-compile
```
It is not a cross-compilation properly speaking, as both scripts first determine the nodes on the
cluster, then ssh into each node and only if a new architecture is found a new, optimized binary
is compiled for that particular architecture. This workaround works as long as **access to the file
system is shared** on the cluster and that **each node can be accessed via ssh by the user**.


### Documentation

#### Requirements
* a quite complete **TeX live** distribution

The usual *ams maths* packages are required to build the documentation. Other packages are standard and should be shipped with any basic TeX live distribution.
The compilation uses ```pdflatex``` and ```bibtex```.

The documentation is built with
```
make doc
```
which creates the ```doc/doc.pdf``` file.


## Folder structure

After building the executables, run
```
bin/download
```
to build the folder structure required by the framework.
The script also downloads the input files from https://pprc.qmul.ac.uk/~tboschi/HK/atmo.

You can specify a specific path with the ```-p prefix``` option, as explained in the usage. The default value is the folder ```errorstudy/``` relative to the current working directory 

The subfolders ```errorstudy/reconstruction_beam``` and ```errorstudy/reconstruction_atmo``` are created to contain reconstruction files to build the data samples.
The subfolders ```errorstudy/systematics_beam``` and ```errorstudy/systematics_atmo``` contain various systematic models that can be used straight away in the analysis.


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
Finally, the fitter can be launched on the distributed computing system with the `trisens`
utilities: it automatically detects the HPC manager.
For example
```
bin/trisens -r errorstudy/first_run -d comb -1 NH -2 NH -N 500
```
will launch 500 jobs (```-N 500```) on the cluster, fitting both the beam and atmospheric samples (```-d comb```) with true and fitted normal mass hierarchies (```-1 NH -2 NH```).
Other important options are
* ```-s``` to do a statistics only fit, i.e. no systematics;
* ```-f <scan_type>``` to perform a scan fit, i.e. chaning the true point at each iteration as it is done for CPV sensitivity studies;
* ```-v <verbosity>``` to change the verbosity of the log files.
The full list of options can be shown with
```
bin/trisens -h
```

## Everything else

Please refer to the documentation.
Most of the scripts have the ```-h``` option to show usage.

## Problems?

If you spot anything off, please reach out to [boschi.tommaso@gmail.com](mailto:boschi.tommaso@gmail.com) or open an [issue](https://github.com/tboschi/SuperHK/issues/new).
