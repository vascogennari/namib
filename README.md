# PYTOP

Python Yielding Tools for Optimization of Posterior distributions (PYTOP) is a python package for handling and displaying *posterior distributions* from gravitational-wave (GW) data analysis.

Given a set of files containing samples, PYTOP will display the posterior distributions for the selected parameters through three different plots: ***corner***, ***violin***, ***ridgeline***.

We show two possible applications as an example.
#### Comparing two pipelines on a set of GW events

#### Comparing two models at different starting times

To run PYTOP, you need a set of samples (input files) and a config file.
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
    
You can get familiar with PYTOP by replicating the above examples and an exercise. To do that, try to launch the config files *config_GWTC-3_LVK-TEOBPM.ini* and *config_ringdown_HMs.ini*
