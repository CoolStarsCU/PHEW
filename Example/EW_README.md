## EW Wrapper Documentation ##

The code “EW_WRAPPER_AUTOMATED.py” is an example of an automated code that’s purpose is to efficiently fit and get equivalent width measurements for a large quantity of spectrum. The code operates on the fact that many spectra can be properly fitted under highly similar parameters. Thus, the code fits spectra in rounds. In each round the code applies a fixed set of parameters for all unfit spectra in the directory, and removes those that were fit properly, leaving the rest to be fit under the new parameters fixed in the next round.

## Setup Steps ##
* Download the readspec.py file included in the example folder.
* Download desired .fits files (A sample set is included in “TEST_PHEW_SAMP”).
* Download "EW_WRAPPER_AUTOMATED.py".
* Open a new script in the same directory and import "EW_WRAPPER_AUTOMATED.py" as a package, from which the functions (below) "fit_runner" and "mcmc_run" can be called.
* Once the functions are called and proper path arguments are provided, the new (user-created) script can be executed via the command line or in the user's IDE of choice.

## Outputs ##

* EQW_FITS.csv
  * Contains the parameters used (or not used) for all spectrum files in the inputted directory. Also contains ‘types’, identifying whether a spectrum could not be read, was deemed to be flat, or remained undetermined.
* EQW_VALS.csv
  * Contains all spectra files in the inputted directory with their equivalent width value and error.

## Running the Code ##

*__All paths are absolute paths (eg /Users/stanislavdelaurentiis/desktop/spectra)__*

## fit_runner(_path_) ##
* Your input “path” here will be the directory containing the spectra.
* The code will prompt you to enter in general parameters that would be used to fit a curve to the spectral feature you are looking for. These general parameters will then be applied to every unfit spectra file in the directory during that round.
* You will then be prompted to make a series of decisions about the quality of the fit. If you are happy with the fit, the program will record the parameters you selected for this spectrum, save the output pdf figure into the same directory, ask for comments, and then show you the next file.
* If you are not happy with the fit you will be asked whether you would like to reapply the parameters (params) for this spectrum only (testing a fit with different parameters than those you chose at the beginning of the round). 
  *  If the amount of files remaining in the directory is still large (greater than 20) I discourage reapplying params individually. It would be more time efficient to wait until the next pass through the directory and adjust the general params accordingly.
  *  If this is not the case, you would apply new params to this individual file (as many times as you wish) and then make a decision on its fit.
* If you do not save the file, the program will then ask if the spectrum is effectively zero (meaning that the spectrum itself seems to lack the spectral feature you are looking for entirely).
  * If you answer ‘y’, the file will be filed as having no discernible feature, and you will be shown the next file.
* If you do not decide that the spectrum has a fit you are happy with or has no feature, the program will ask you whether you want to continue attempting to fit this spectrum.
  * If you do not wish to continue attempting to fit this spectrum, the code will file it as being undetermined, and will not refit it in the next round.
  * If you do wish to continue attempting to fit this spectrum, it will hold on to this spectrum as being undecided for the time being, including it in the next round.
* The program will then continue doing this for each spectrum for each pass through the loop, until there are no spectra left undecided in the directory you inputted. 
* In the original directory with the spectra you will find that there is an ‘EQW_FITS.csv’ file that has recorded each spectrum and its parameters.

## mcmc_run(_path_) ##

*__This function is fully automatic__*

* The inputted directory is the same as before.
* The program will then run through all the spectra that were fit with an emission, absorption, or labeled as having zero equivalent width, and apply an MCMC routine of 1000 steps to get a precise equivalent width and error.
* The code will then output the according pdf figures into the directory with the spectra, as well as recording the ew values in the “EQW_VALS_ALL.csv” file.



