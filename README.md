# Easy Osc3++

The simplified version of osc3++ is described by the ```event/ChiSquared``` class.
The main code in ```app/fitter.cpp``` uses the class to load systematics, build distributions and minimise the chi-squared.
It fits only the beam sample.  The ```app/atmofitter.cpp``` processes the atmospheric sample, but it is under development.
The oscillation space is scanned over and it is defined/managed by the ```physics/ParameterSpace``` class.

The space is specified by combinations of oscillation parameters defined inside ```cards/fit.card```.
For each point of the oscillation parameter space, a set of histograms is created and oscillated with those oscilliation parameters.
The main set available is *Asimov A*, labelled as `asim`.
The true point is defiend by the *Point* variable: in the default space, point 49773 is Asimov A.

A set of script is used to facilitate the execution of the code and extraction of results.
In exectuing the following scripts, the structure of the folder is important.
We will refer to paths where input files to be fitted are created as ```global``` (data only for atmospheric sample, but configuration files for beam sample too), and as ```root``` to where the outputs of the fitted file will be.

## Folder structure

Let's say the main path for the studies is the ```errorstudy/``` folder.
Create a folder for the global inputs and parameter set (e.g. asim)
```
mkdir -p errorstudy/global/asim
mkdir -p errorstudy/global/reconstruction
```
When dealing with atmospheric samples, ```asim``` will contain subfolders with oscillated files to be fitted, but most importantly ```.info``` files with a list of true points to fit: **point.info** to generate a simple chi-squared profile against one single true point, **scan.info** with all the points for a sensitivity sweep.
The folder ```reconstruction``` keeps configuration files and root files used to build the beam observables. From VALOR.


Create also folder for each set to be studied
```
mkdir errorstudy/root/systematics/
```
and the sub-directory ```systematics``` contains spline root files and the correlation matrix.

The idea is that there is a single global folder with common parameters and multiple root folders, one for each systematic set.


## Creating the fitting space for atmospheric sample

**WARNING** *it uses ```GlobalOsc``` part of the original Osc3++ package. Skip if not interested in atmospheric stuff or do not have GlobalOsc installed.*

The first step is to create the fitting space, by running ```GlobalOsc```.
This can be done with the follwing script:
```
./global.sh [-i | -n] -f nFiles -g /path/to/global
```
where ```-i``` or ```-n``` select inverted or normal hiereachy for the oscillation.

If using batch jobs, the following script works with HTCondor.
```
./global_c.sh [-i | -n] -f nFiles -g /path/to/global
```


## Prepare the systematics
Creates correlation matrix of systematic parameters by combining matrices found in files of matrixn.root.
The script expects to find the systematics folder under root with the spline files to be processed (renaming of files and histograms).
It uses ```app/purifysystematics.cpp``` and ```app/addmatrix.cpp```.
```
./prepare_systematics.sh -r errorstudy/root matrix1.root [matrix2.root ...]
```

# Run the fit

The syntax is
```
./bin/fitter id all $output/this_sensitivity.card$1 $2 $3
```
The fitter is meant to be used with some parallel computing, like a batch system.
For this reason, the first two inputs are telling the process which set of files (or points) to study ```id``` and how many processes are running ```all```.
This maximises the effort from the different CPU running. 
The last parameter is the configuration file which specifies oscillation, samples, fitter, ecc.

The easiest way to run the fit is using the following script, which handles parallel computing, folder structure, and configuration files.
It work with HTCondor, but can be easily adapated to any other manager.
```
./trisens_c.sh -g errorstudy/global -r errorstudy/root/asim -1 [NH | IH] -2 [NH | IH] [-f] [-s]
```
where options 1 and 2 specify the mass hiererachies of respectively the observed and expected event samples.
The script uses the point specified in ```global/asim/point.info``` for the fit.
The option -s does a stats only fit (no systematics). The option -f performs the scan of multiple true points, using the ones specified in ```global/asim/scan.info```


The study multiple sets at the same time, the scheduler ```launch_trisens_c.sh``` working with HTCondor was devised.
Manually modify the model array and the point variable.



## Extract results

Once the fit is finished, chi-square penalty terms can be optionally added and the profiles created.
In the case of a scanning fit, the CP exclusion is also computed.

To penalise for instance the result of an asim NH vs NH fit, run the script
```
./penalise.sh -r errorstudy/root/asim/NH_NH/
```
The options are specified in the cards/penalty_asim.card.
The executable ```app/addpenalty.cpp``` is used.


To build the contour for point 49773 of an asim NH vs NH fit, run the script
```
./contours.sh -p -r errorstudy/root/asim/NH_NH point_49773
```
the name ```point_49773``` is the name of a folder (automatically generated), it is not a string.
The -p option specifies to use penalised files.

To create the sensitivity plot
```
./excludeall.sh -p -r errorstudy/root/asim/NH_NH
```

These operations could be done more easily and on muliple sets with the scheduler ```launch_contours_c.sh```.
It uses HTCondor and it works similarly to the fitter scheduler.
