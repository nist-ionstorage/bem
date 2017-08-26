Windows
=======

John Gaebler, 2015

These instructions have been tested on Windows 7 64bit

Note: It may be worth trying PythonXY, which comes with the vtk, mayavi, and
cvxopt packages, however pythonxy seems more apt to become confused with
multiple versions of python. Also, PythonXY is only available with 32 bit
Python.

Obtain BEM code:

  1. Download and install Git for Windows https://windows.github.com/

  2. Generate a .ssh key by following insturctions at:
       https://help.github.com/articles/generating-ssh-keys/ No password is
       needed (leave blank)

  3. Send the .ssh key to Robert

  4. Now you should be able to open a Git shell and type "git clone
       git@ions:bem.git" and give password

Obtain Electrode code:

  1. Go to https://github.com/nist-ionstorage/electrode

  2. Click the Clone in Desktop button

Install Anaconda version 2.7 (64 bit)

  1. http://continuum.io/downloads

Install Microsoft Visual Studio

  1. http://www.visualstudio.com/

Install Python Packages

  1. Open a Windows command prompt

  2. Install vtk - Type "conda install vtk"

  3. Install mayavi  - Type "conda install mayavi"

  4. Install cvxopt for electrode - Doesn't seem to work in conda64 - need to
       retest in conda32

Install bem package

  1. Open an ipython command prompt from anaconda

  2. Navigate to bem folder (use cd .. to move up dir, ls to see folders, and
       cd "foldername" to go to folder)

  2. type !python setup.py develop

Install Electrode Package

  1. In ipython command prompt Navigate to electrode folder

  2. type !python setup.py develop


Run Something:

  1. Open Spyder from Anaconda Startup Folder

  2. Run simple trap example from BEM code (Electrode part will fail)

  3. Note, the example needs to be changed so that you do not pass analyze
       static the u parameter


Linux
=====

After setting up git (get the ssh key to Robert) and installing numpy, scipy, myavi2 (for tvtk) and electrode if desiered , do the usual::

  git clone git@ions:bem.git
  cd bem
  python3 setup.py test
  pip3 install --user -e .
  examples/SimpleTrap/SimpleTrap.py
