# Channel Selection in a least-squares (LS) problem 
This project has code implementing channel-selection in an LS problem. Originally implemented for channel-selection in auditory attention decoding (AAD) based on EEG [1], the functions can be used in selecting relevant channels in an LS problem using any multi-channel signal.


## Features
Group-utility based channel Selection: Selects the _best N_ channels of A in the following LS minimization problem 

	min_w ||Aw - b||^2
	
	where, A is (T X M) matrix. w is a (M X 1) filter. b is the desired (T X 1) signal which we are
	looking to reconstruct using the solution of the above problem.
	
The best channels are selected based on group-utility, where group-utility is defined as the increase in mean-squared error (MSE) when a group of channels (i.e. columns of _A_) are removed from the problem. 
a. MATLAB version: channel_select.m

b. Python version: channel_select.py

	


## References
[1] A. Mundanad Narayanan and A. Bertrand "Analysis of miniaturization effects and channel selection strategies for EEG sensor networks with application to auditory attention detection" IEEE Transactions on Biomedical Enginnering, 2019

[2] Bertrand, A. (2018). Utility Metrics for Assessment and Subset Selection of Input Variables for Linear Estimation [Tips & Tricks]. IEEE Signal Processing Magazine, 35(6), 93–99.
