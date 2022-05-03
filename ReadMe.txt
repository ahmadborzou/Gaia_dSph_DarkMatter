==========
How to run:

1. choose the name of dwarf spheroidal galaxy, its distance from the earth in kilo-parsec and the error of its distance from the earth in kilo-parsec. For example use table 1 of arXiv:0908.2995 [astro-ph.CO]. These will be needed below.  
An example would be: 
Fornax 147 3

2. 

=================
Code dependencies:

The needed libraries to be installed can be found in environment.yml 

One Gaia specific package not in this file is zero_point 
which can be installed for example using (see also: https://pypi.org/project/gaiadr3-zeropoint/)
$ pip install gaiadr3-zeropoint

File requirements.txt includes the full list of all the installed packages and their corresponding versions.  


To propagate errors, we need another package that can be installed using (see also: https://pythonhosted.org/uncertainties/)
$ pip install --upgrade uncertainties


Also, see runtime.txt for the python version.



