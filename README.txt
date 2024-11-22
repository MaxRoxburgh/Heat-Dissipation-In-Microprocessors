Hi,

Hope you are having a lovley day. I would first like to apologise for my variable names as I am dyslexic and have probably made a lot of spelling mistakes in the names.
This code was created as coursework for a computational physics course at Imperial College London. 
The shortened report submitted with the code can also be found here; A huge amount of design, implementation testing, validation, criticism of the model and references have been cut due to a word limit on the work.
Please message me if you have any questions.

In terms of using the code:

	I could not attach all of the inverse matricies that I had calculated with my code as the folder was over 115 GB...
	For this reason, if you wish to run my code, the inverse matricies will need to be recalculated.
	In the EXAMPLE_OF_PARRALEL_RUN.py file, it automatically calculates the necisary matricies which was why I made the file.
	If you wish to do ones yoursef, go to the _M_inverse_generator.py file and run the code, or take that format and make your own inverses


To get all of the plots and analysis used in the report, run the aptly named RUN_ME_FOR_REPORT_AND_ANALYSIS.py as you probably guessed by the title
The IPYNBS folder contains all early work from the project that was done on jupiter notebook

All functions utilised in the code are stored in the py files that DO NOT have an _ at the begining of the name

The primary file that investigations were run from was _parrallel_runs_until_convergence.py but that file is a mess and I would not recomend using it
all code from EXAMPLE_OF_PARRALEL_RUN.py is taken from this origonal file ^^

The data saved from the runs that were analysed for the majority of the report are in whole_system_DFO4
for each run there is a seperate data file with the temperature and difference between each run
and on the final run, the entire system and its segents are saved as seperate csv files.

All temperatures used to simulate are scaled by a factor of 1/293 or the ambient temperature modeled to be outside of the system.
This implies if you see a value of T = 21, the the temperature in kelvin is 21*293 or 8500

The primary function used to run the simulations is in Updates and is called update_system_n_times()
The subscripts in the variable names:
	T_0 implies it is the boundary values and has 0's everywhere eles
	T_k implies the value of temperature inside the block being simulated, this inclues the boudnary values
	_pr/_cs/_hs/_fn are processor, ceramic case, heat sink body and the list of fins respectivley

The list of fins should be half the size of the total number of fins as arguments of symetry were used to shorten compute time dramatically 
Finally the data from the last run with the optimised shape for minimum dimensions is found in:
	Heat_Dissipation_In_Microprocessors\whole_system_FDO4\fin_length_test_forced\data_s6_fin_length_39
	after increasing the processor temperature achieved due to scaling it got to 81 degrees kelvin

