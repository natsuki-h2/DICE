
# GenerateMockData

+ This code outputs the comparison of the ways to inject artificial source into the real/simulated imagecube or visibility data taken with ALMA.


+ How to run: in CASA, just type
		 execfile('GenerateMockdata.py')


+ The directory structure:
	GenerateMockdata
		- README (we're here)
	  - GenerateMockdata.py
		- Template
				- configurationfiles
				- ptgfiles
				- MesurementSet


+ Notes:
	- configuration files are created by using executable 'aU.buildConfigurationFile'
	- pointing data is referred from QA2 report. QA2 report is downloaded from ALMA archive query. 
	- datailed description and structure of the code is written in the code (GenerateMockdata.py)
	- there is inconsistency of RMS value of noise between the simulated data and original data if the data is combined. 


+ HISTORY:
	22 March 2018 share the code with project member by Natsuki Hayatsu
	
========

BliSS-ALMA grew out of an earlier repository, GenerateMockData, which was originally developed in 2018 for ALMA mock-data generation and artificial-source injection tests.
