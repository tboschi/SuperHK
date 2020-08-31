# Easy Osc3++ - v3.0

The core of this simplified version of osc3++ is given by the classes in the ```event``` folder.
The main code in ```app/fitter.cpp``` uses derived classes inheriting from ```Sample.h``` to load systematics and build distributions, whereas the ChiSquared class computes and and minimises the chi-squared.
Support to beam sample and atmospheric sample is available.
The oscillation space is scanned over and it is defined/managed by the ```physics/ParameterSpace``` class.
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

* Since version v2.0 the energy scaling error is handeled more exactly because it is implemented as an anlaytic function.  Hence, the mathematical exact Jacobian and Hessians are used during minimisation.
* The Oscillator and ChiSquared classes have been remodelled to improve optimization and maximise use of vector instructions.
* The chi-squared penalisation term is automatically added and it can be specified in the ```cards/oscillation.card```.
* The atmospheric sample is now supported.
* The CardDealer class, to import parameters dynamically from text files, has been rewritten as a template class to support as many data types as possible.

## Requirements

ROOT should be installed and properly linked. To test if true, simply run
```
root-config --cflags --glibs
```

You must download and extract the [Eigen 3 library](https://eigen.tuxfamily.org/dox/index.html) and modify the ```EIGENINC``` variable in the ```Makefile``` to point to the Eigen folder.

## Build

Running
```
make
```
will create the executables in the ```bin/``` folder.


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
Creates correlation matrix of systematic parameters by combining matrices found in files of matrixN.root.
The script expects to find the systematics folder under root with the spline files to be processed (renaming of files and histograms).
It uses ```app/purifysystematics.cpp``` and ```app/addmatrix.cpp```.
```
./prepare_systematics.sh -r errorstudy/root matrix1.root [matrix2.root ...]
```

**N.B** this is done automatically by the ```setup.py``` script!


# Run the fit

The syntax is
```
./bin/fitter id all $output/this_sensitivity.card
```
The fitter is meant to be used with some parallel computing, like a batch system.
For this reason, the first two inputs are telling the process which set of files (or points) to study ```id``` and how many processes are running ```all``` in total; this helps the executable to know which points study.
The last parameter is a configuration file which specifies links to cards for oscillation, samples, fitter, ecc.

The easiest way to run the fit is using the following script, which handles parallel computing, folder structure, and configuration files.
If your manager is HTCondor
```
./trisens_c.sh -r errorstudy/root/asim -1 [NH | IH] -2 [NH | IH] [-f] [-s] [-v <verbosity>]
```
or if your manager is SLURM
```
./trisens_s.sh -r errorstudy/root/asim -1 [NH | IH] -2 [NH | IH] [-f] [-s] [-v <verbosity>]
```
where options 1 and 2 specify the mass hierarchies of respectively the observed and expected event samples.
The script uses the ```card/oscillation.card``` to determine what is the point to fit.
The option -s does a stats only fit (no systematics). The option -f performs the scan of multiple true points *WARNING: this feature has not been tested yet*

The study multiple sets at the same time, the scheduler ```launch_trisens_c.sh``` working with HTCondor or ```launch_trisens_s.sh``` was devised.
Manually modify the model array and the point variable.

Assuming for example, a fit on the asimov A point 49773 with normal hierarchy for both the observed and expected samples,
the fit result will be organised under the following folder
```
errorstudy/root/NH_NH/sensitivity/point_49773/
```
The folder will contain also log files of the fitter, the cards used, and the batch jobs scripts, besides the output root files.
There will be a root file per job, each containing a tree with the computed chi-squared.
The default output name will be *SpaghettiSens.T2HK.XXX.root*, where XXX is a serial file number.



## Extract results


To build the contour for point 49773 of a NH vs NH fit, run the script
```
./contours.sh -p -r errorstudy/root/NH_NH point_49773
```
the name ```point_49773``` is the name of a folder (automatically generated), it is not a string.
The -p option specifies to use penalised files.
The contours are saved under
```
errorstudy/asim/NH_NH/contours/point_49773/
```

To create the sensitivity plot
```
./excludeall.sh -p -r errorstudy/root/asim/NH_NH
```
The exclusion file is saved under
```
errorstudy/asim/NH_NH/exclusions/
```

These operations could be done more easily and on multiple sets with the scheduler ```launch_contours_c.sh```.
It uses HTCondor and it works similarly to the fitter scheduler.

## Producing plots

In order to produce plots with gnuplot, the following script should be called first
```
./dropchi2.sh -p point_49773 -n name_plots -r errorstudy/ root0/asim/NH_NH/ [root1/asim/NH_NH/  ... ]
```
This will create text file with all the data necessary to make chi-squared profiles, contours, and sensitivity plots.
The -p option specifies which fitted point to use, the -n option is an output name and a folder with the same name will be created to collect all the files.
The -r option specifies the root path, after which a list of root folders can be specified: the created text file will contain also relative data with respect to the first root folder. This is useful when comparing different systematic sets.

Run 
```
./create_plots.sh name_plots
```
to execute gnuplot scripts.
By default pdf plots are generated with the text exported to latex files.
