To run Brian's BF RV analysis code in IDL:

1.  The idlastro libraries must be installed on your system and accessible by IDL (e.g., via the IDL_PATH).
2.  In IDL, do the following compiles:
.comp BFall_IDL.pro
.comp read_spectrum.pro
.comp my_bf_script
resolve_all
.comp make_and_fit_broadeningfuncs
3.  Run the main program.
make_and_fit_broadeningfuncs

Adjustments you must make to the code:

On about line 21, you'll need to adjust the directory names (the
observation_dirnames array) to point to where
you've stored the HET spectrum files and the ARCES spectrum files.  Inside
these directories, IDL will search for file name formats in the
observation_filebasenames array.  You will need to adjust these too if the
files are called something else.  The program will also be searching for the spectrum file names with specific order numbers;  the order numbers to be looped through are controlled via the variables observation_starting_orders and observation_n_orders. 

Optional adjustments:

In section 4.2 (about lines 337-392), the starting guesses for the model fitting parameters are specified.  MPFIT seems like it doesn't vary the instrumental LSF widths very much, so you may have to adjust the starting guesses for instrumental_lsf_sigma to get a more reasonable fit.  Each epoch gets a separate instrumental_lsf_sigma, so you may need to adjust up to six starting guesses here.