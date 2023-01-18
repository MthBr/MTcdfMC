PDF_project
==============================

#https://stackoverflow.com/questions/501940/simple-simulations-for-physics-in-python

Multi-dimensional delta-Dirac Approximations for Uncertainty Quantification of Flows through Randomly Heterogeneous Porous Media



conda update -n base -c defaults conda


conda clean --packages --tarballs
conda clean --all
conda update --all
conda update conda
conda update conda-build


conda list --revisions

conda activate base
conda install --revision 0

conda info
conda config --show-sources
conda list --show-channel-urls




# how to create developing environment [env]

conda env create -f cdfARE-environment.yml 
NOTE: it may take some time to solve all the rependencies!

conda env remove -n  cdf-env
conda env list

conda activate cdf-env
cd CDF Code
pip install -e .

conda deactivate


# AMGX#bindings
cd AMGX-2.2.0
mkdir build
cd build
cmake ../
make -j16 all


#TODO
#nvcc warning : The 'compute_35', 'compute_37', 'compute_50', 'sm_35', 'sm_37' and 'sm_50' architectures are deprecated, and may be removed in a future release (Use -Wno-deprecated-gpu-targets to suppress warning).




# Building and installing pyamgx
conda activate cdf-env
export AMGX_DIR=/home/modal/Projets_local/AMGX-2.2.0
# export AMGX_BUILD_DIR=$AMGX_DIR/build
cd pyamgx-main
pip install .





# how to create developing environment [env]
pip install pysparse





# how to create developing environment [env]
conda env remove -n fenicsproject
conda activate fenicsproject
conda install -c conda-forge spyder 
spyder &


conda install opencv == 3.4.12
conda install tensorflow==1.14
conda install tensorflow-gpu==1.13.1

- "opencv>=4.5" 

channels:
  - defaults
  - conda-forge
  
"opencv-python>=4.4"

"opencv-contrib-python-headless>=4.4"


anaconda-navigator




# on server
conda activate kg-env ; spyder --new-instance &


# Extra
conda update conda
conda update anaconda
conda update python
conda update --all

conda clean --all





# Tips for developers
pip install pep8
pip install pylint


advice in 2020, to use Visual Studio Code

Install  Python support


VS Code Quick Open (Ctrl+P)

ext install ms-python.python



# Tips for VS code
Extension:
Python  (ms-python.python)
GitLens — Git supercharged
Code Spell Checker


# Configs for VS code
File > Preferences > keyboard shortcut
F5: file Python interactive
F9: line Python interactive

File > AutoSave


Press: CTRL + Shift + P

Click on "Preferences: Open Settings (JSON)"

Add this line into JSON : "python.linting.pylintArgs": ["--generate-members"]


#https://stackoverflow.com/questions/56844378/pylint-no-member-issue-but-code-still-works-vscode
  
  
  
  
  
  
  
  
name: cdf-env
channels:
  #- conda-forge
  - defaults
dependencies:
  - pathlib
  - "pip>=19.3"
  - "python>=3.7" 
  - pandas 
  - jupyter 
  - seaborn 
  - matplotlib
  - numpy
  - pep8
  - pylint
  - scipy
  #- gstools
  #- pysparse
  #- fenics
  #Pysparse, SciPy or Trilinos
  - pip:
      - gstools
      #- fipy
      - -e .
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
