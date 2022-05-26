# CSF_SSFP

% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA

The codes are solely for demonstration purpose and have not been optimized
for performance. 

Contents: 
a_Dictionary: To generate a dictionary
	* DictGen_main.m: The main function to genearte a new dictionary and
	prepare it for lookup
	* FlowVelocity.m: A function to generate a sinusoid velocity waveform
	with configurable properties
	* SSFP_Flow_Simu.m: A function to perform Bloch simulations. A single
	period of the velocity waveform is repeatedly used until convergence
	has been achieved.
	
b_MC_Simulation: To perform Monte-Carlo simulations
	* MC_main.m: The main function to create 3-sinusoid waveforms,
	simulate the perturbation patterns, and match to the dictionary.
	* SSFP_Flow_Simu_Aperiodic.m: A function to perform Bloch simulations
	using continuous velocity waveform.
	
c_InvivoValidation: To perform quantification in vivo at 3 CSF locations.

