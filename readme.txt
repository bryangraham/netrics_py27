netrics: a Python 2.7 package for econometric analysis of networks
-----------------------------------------------------------------------------
by Bryan S. Graham, UC - Berkeley, e-mail: bgraham@econ.berkeley.edu


This package includes a Python 2.7 implementation of the two econometric
network formation models introduced in Graham (2014, NBER).

This package is offered "as is", without warranty, implicit or otherwise. While I would
appreciate bug reports, suggestions for improvements and so on, I am unable to provide any
meaningful user-support. Please e-mail me at bgraham@econ.berkeley.edu

Please cite both the code and the underlying source articles listed below when using this 
code in your research.

A simple example script to get started is::

	>>>> # Append location of netrics module root directory to systems path
	>>>> # NOTE: Only required if netrics not "permanently" installed
	>>>> import sys
	>>>> sys.path.append('/Users/bgraham/Dropbox/Sites/software/netrics/')

	>>>> # Load netrics package
	>>>> import netrics as netrics
	
	>>>> # View help file
	>>>> help(netrics.tetrad_logit)
	
	>>>> # Fit link formation model
	>>>> # See help for how to construct D and W in practice
	>>>> [beta_TL, vcov_beta_TL, tetrad_frac_TL, success] = \
            netrics.tetrad_logit(D, W, dtcon=None, silent=False, W_names=cov_names)


CODE CITATION
---------------
Graham, Bryan S. (2016). "netrics: a Python 2.7  package for econometric analysis of networks," (Version 0.0.1) 
	[Computer program]. Available at https://github.com/bryangraham/netrics (Accessed 04 September 2016) 
	
PAPER CITATIONS
---------------
Graham, Bryan S. (2014). "An econometric model of link formation with degree heterogeneity," NBER Working Paper No. w20341.	