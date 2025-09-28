# PDE-IF-Estimation-MulticomponentSignals
(Partial Differential Equation-Based Instantaneous Frequency Estimation for Multicomponent Signals)

## Overview
This repository contains the MATLAB implementation of the PDE-based instantaneous frequency (IF) estimation method described in our paper *[ A Partial Differential Equation-Based Method for Instantaneous Frequency Curves Estimation of Multicomponent Signals in the Time–Frequency Spectrogram]*.  
The method is designed to estimate instantaneous frequency for multi-component signals under **mode-crossing** and **spectral-mixing** conditions. By optimizing spectrogram gradients, adaptive frequency point selection, and multi-scale regularization, the method ensures robust IF trajectory estimation.

## Repository Contents
- `Code/` : Core MATLAB functions  
  - `Phi1Estimator.m`  
  - `Phi1Estimator4.m`  
  - `selectPointsBySlope1Peak.m`  
  - `selectPointsBySlope1PeakAdaptive.m`  
  - `tvdenoise.m`  
- `Demo/` : MATLAB Live Script demonstrating the full workflow  
  - `untitled23_3.mlx` (generates synthetic signals, runs the IF estimation, and produces figures)  
- `LICENSE` : Apache-2.0 license  
- `README.md` : This file  

## Requirements
- MATLAB R2024a (or later)  
- Signal Processing Toolbox   

## Usage
1. Open MATLAB and navigate to the repository root.  
2. Run the demo script:  
   ```matlab
   cd Demo
   untitled23_3
