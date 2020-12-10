# XD_TDS

### What & Why
In the paper [Empirical correction for resolution- and temperature-dependent errors caused by factors such as thermal diffuse scattering](https://scripts.iucr.org/cgi-bin/paper?ks5474) the Authors suggest to use a polynomial expression to counter the effect of TDS of the form *alpha* = *a*\*x&sup2; + *b*\*x&sup3; as a function of sin(*&theta;*)/*&lambda;*. The program will set-up a X [12] scale-factor refinement for a given set of *.mas*, *.inp* and *.hkl* files and refine X [2] cycles. After every cycle the polynomial is fitted to the refined scale factors, optimal parameters for *a* and *b* are determined and a new ‘corrected’ *.hkl* file is written. This nested-intervals approach is necessary as the model needs to relax to this new, more correct data but two-three cycles seemed sufficient. The problem following this approach is that the overall correction factor is no longer directly obtainable. Therefore, the overall sum of the iteratively applied *a*, *b* parameters are finally used to write a *.hkl* file and a refinement is performed on those data. The iterative correction differs only marginally from the combined correction, however, with the combined approach we can report the actual parameters used to correct the data.

### How to
All that is needed is to put the latest *.inp*, *.hkl* and *.mas* files into a folder and run it.
The Masterfile should be set-up to at least refine the non-Hydrogen *Uij* and *monopoles* of **all** atoms, *multipoles* and *Kappa* would help in determining better parameters but take a lot more time to refine, please perform different TDS refinements and check for consistency! 
The parameters *a* and *b* used to correct an hkl file are written to the *.hkl* and the line: **!TDS CORRECTION FACTOR: a=x.xxx, b=x.xxx**
is added and should remain there for later reference! XD will correctly ignore the line.
The *XD path*, *number of scale factors* and *number of cycles* are of course hardcoded and have to be changed by editing the *.py* file directly.

### The small note
The program needs the initial *.hkl* file to have 6 columns, that should be default anyways.

### Output (after the program completed!)
*xd00.\** - Initial 1-scale refinement.
 
*xd01.\** - Initial X-scale refinement.
 
*xd02.\** - X-scale cycles.
 
*xd0X.\** - some more X-scale cycles.
 
*xd.\** - final corrected *.hkl* and 1-scale *.mas* and *.inp* files.

## Requirements

#### [Python](https://www.python.org/) 3.6 or later (f-strings!)

#### Libraries:
  - [numpy](https://www.numpy.org/)
  - [pandas](https://pandas.pydata.org/)
  - [scipy](https://www.scipy.org/)
