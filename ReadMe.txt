==========
How to run:

1. open a jupyter notebook and open the following file: run.ipynb

2. When asked insert the name of the dwarf spheroidal galaxy, its distance from the earth in kilo-parsec and the error of its distance from the earth in kilo-parsec. For example use table 1 of arXiv:0908.2995 [astro-ph.CO]. An example would be: 
Fornax 147 3

3. The code retrieves data from Gaia using GetGaia and store the files under a new folder with the name of the dSph galaxy. It will copy the final output into AnalyzeGaia folder for further processing. In this folder, a few informative histograms of the stars will be shown and stored. Also, the positions and velocities in the co-moving frame will be computed. The phase-space density estimation will be done.

4. The code will ask for an x and y of a position in the co-moving frame (from the center of the dSph) in parsec (pc). Examples are x = 150 & y = 200

5. The mass density of DM, and its static and systematic errors, at that (x,y) location will be returned. 


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


