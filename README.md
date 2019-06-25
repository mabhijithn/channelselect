---
layout: post
mathjax: true
title: Readme
---
# WESN Tools
This project holds useful functions in emulation and analysis of Wireless EEG Sensor Networks [1].

a given wire happens to be carrying "$$\lvert 0\rangle$$."
By that we mean that it's carrying the linear combination
$$\begin{psmallmatrix} 1 \\ 0 \end{psmallmatrix}$$

## Features
1. Channel Selection: Selects the best N channels of A for the LS problem 

	$$ \min_{\mathbf{w}}\frac{1}{2}||\mathbf{A}_P\mathbf{w} - \mathbf{s}_a||_{2}^{2} $$

a. MATLAB version: channel_select.m

b. Python version: channel_select.py

	Currently supported methods: 
		1. Utility based method.


## References
[1] A. Mundanad Narayanan and A. Bertrand "Analysis of miniaturization effects and channel selection strategies for EEG sensor networks with application to auditory attention detection" IEEE Transactions on Biomedical Enginnering, 2019

[2] Bertrand, A. (2018). Utility Metrics for Assessment and Subset Selection of Input Variables for Linear Estimation [Tips & Tricks]. IEEE Signal Processing Magazine, 35(6), 93–99.
