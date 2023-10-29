# PYTOP

Python Yielding Tools for Optimization of Posterior distributions (PYTOP) is a python package for handling and displaying *posterior distributions* from gravitational-wave (GW) data analysis.

Given a set of files containing samples, PYTOP will collect the posterior distributions for the selected parameters and visualize them through three different plots: ***corner***, ***violin***, ***ridgeline***.

We show two possible applications as an example.
#### Comparing different pipelines on a set of GW events
<img src="https://github.com/vascogennari/PYTOP/assets/62053184/0753647d-03e0-4592-bfd6-5670b8614656" alt="drawing" width="750"/>

#### Model selection between two competing models at different starting times
<img src="https://github.com/vascogennari/PYTOP/assets/62053184/ac453ba8-150f-463e-a091-dbdc24218bd5" alt="drawing" width="750"/>

To run PYTOP, you only need a set of samples (input files) and a config file.
- The name of the input files needs to contain the values of the *fields* following the structure:<br> `event_pipeline_model_submodel_time_GR-tag.ext`
- Choose a set of parameters and one *field* to compare the different files.
- Select the types of plots and customize bounds, dimensions and colors of the figures.
- Enjoy the results.


## Installing and running the code
To install PYTOP, follow the instructions:

    git clone https://github.com/vascogennari/PYTOP.git
    cd PYTOP
    pip install -r requirements.txt
    python setup.py install
    
Once the package is installed, you can launch PYTOP as:

    PYTOP --config-file config_filename.ini
    
You can get familiar with PYTOP, try to replicate the above plots by launching the config files in the *config_files/examples* directory.
