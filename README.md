# XD_TDS

### What & Why
 In the paper they suggest to use a polynomial expression to counter the effect of TDS of the form alpha = a*x^2 + b*x^3 as a function of sin(theta)/lambda. The program will set- up a 12 scale factor refinement for a given set of .mas, .inp and .hkl files and refine three cycles. After every cycle the polynomial is fitted to the refined scale factors, optimal parameters for a and b are determined, the polynomial with those parameters is used to write a new ‘corrected’ .hkl file. This nested-intervals approach is necessary as the model needs to relax to the new more correct data but three cycles seemed sufficient, the problem is however that the overall correction factor is not directly obtainable. Therefore, the overall sum of the iteratively applied a, b parameters are then used to write the final .hkl file and a final refinement is performed. The iterative correction differs only marginally from the combined correction, however with the combined approach we can report the actual parameters used to correct the data.

### How to
 All you need to do is put your latest .inp, .hkl and .mas files into a folder together with the TDS.py program and run the program using python3 (rename the attached TDS.py_bak to TDS.py).
 The Masterfile should be set-up to at least refine the non-Hydrogen Uij and monopoles of all atoms, multipoles and Kappa would help in determining better parameters but take a lot more time to refine, please perform different TDS refinements and check for consistency! 
 The parameters a and b used to correct an hkl file are written to the .hkl and the line: !TDS CORRECTION FACTOR: a=x.xxx, b=x.xxx
is added and should remain there for later reference! XD will correctly ignore the line.

### The big note
 It seems that there is an issue with the very high resolution data as for that data the Fo/Fc vs. resolution distribution does not follow the expected TDS slope, e.g. the curve  falls off towards unity again instead of keeping to increase continuously. I had to cut the data (in the .mas file: SKIP   obs  0.0E+00 0.1E+11 *sigobs  3. 0.1E+07 *sinthl  0.000  1.200) and I am not sure if this is needed for your remaining data as well. Please check the DRK plot for the slope and set the sintl cut-off to where (if) the slope starts to  decrease.

### The smaller note
 The program needs the initial .hkl file to have 6 columns, that should be default anyways.

### Output (after the program completed!)
 xd00.* - Initial 1-scale refinement.
 
 xd01.* - Initial 12-scale refinement.
 
 xd02.* - 12-scale cycles.
 
 xd0X.* - more 12-scale cycles.
 
 xd.* - final corrected .hkl and 1-scale .mas and .inp files.
