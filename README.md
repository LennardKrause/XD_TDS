# XD_TDS

### What & Why
 In the paper [Empirical correction for resolution- and temperature-dependent errors caused by factors such as thermal diffuse scattering](https://scripts.iucr.org/cgi-bin/paper?ks5474) the Authors suggest to use a polynomial expression to counter the effect of TDS of the form *alpha* = *a*\*x&sup2; + *b*\*x&sup3; as a function of sin(&theta;)/&lambda;. The program will set-up a 12 scale-factor refinement for a given set of *.mas*, *.inp* and *.hkl* files and refine x [2] cycles. After every cycle the polynomial is fitted to the refined scale factors, optimal parameters for *a* and *b* are determined, the polynomial with those parameters is used to write a new ‘corrected’ *.hkl* file. This nested-intervals approach is necessary as the model needs to relax to the new more correct data but two cycles seemed sufficient, the problem  following this approach is, however, that the overall correction factor is not directly obtainable. Therefore, the overall sum of the iteratively applied *a*, *b* parameters are then used to write a final *.hkl* file and a refinement is performed on those data. The iterative correction differs only marginally from the combined correction, however, with the combined approach we can report the actual parameters used to correct the data.

### How to
 All that is needed is to put the latest *.inp*, *.hkl* and *.mas* files into a folder and run it.
 The Masterfile should be set-up to at least refine the non-Hydrogen *Uij* and *monopoles* of **all** atoms, *multipoles* and *Kappa* would help in determining better parameters but take a lot more time to refine, please perform different TDS refinements and check for consistency! 
 The parameters *a* and *b* used to correct an hkl file are written to the *.hkl* and the line: **!TDS CORRECTION FACTOR: a=x.xxx, b=x.xxx**
is added and should remain there for later reference! XD will correctly ignore the line.

### The small note
 The program needs the initial *.hkl* file to have 6 columns, that should be default anyways.

### Output (after the program completed!)
 xd00.* - Initial 1-scale refinement.
 
 xd01.* - Initial 12-scale refinement.
 
 xd02.* - 12-scale cycles.
 
 xd0X.* - some more 12-scale cycles.
 
 xd.* - final corrected .hkl and 1-scale .mas and .inp files.
