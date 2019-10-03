# WESN Tools
This project holds useful functions in emulation and analysis of Wireless EEG Sensor Networks [1].


## Features
1. Channel Selection: Selects the best N channels of A  which minimizes the following LS problem 

	min_w ||Aw - b||^2
	
	where, A is (T X M) matrix. w is a (M X 1) filter. b is the desired (T X 1) signal which we are
	looking to reconstruct using the solution of the above problem.
	
	Currently supported methods: 
		1. Utility-based channel selection

a. MATLAB version: channel_select.m

b. Python version: channel_select.py

	


## References
[1] A. Mundanad Narayanan and A. Bertrand "Analysis of miniaturization effects and channel selection strategies for EEG sensor networks with application to auditory attention detection" IEEE Transactions on Biomedical Enginnering, 2019

[2] Bertrand, A. (2018). Utility Metrics for Assessment and Subset Selection of Input Variables for Linear Estimation [Tips & Tricks]. IEEE Signal Processing Magazine, 35(6), 93–99.
