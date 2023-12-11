# `namib`

`namib` is a python package for handling and visualizing multiple *posterior distributions*.

Given a set of input files with samples, `namib` filters the samples for the selected parameters and displays them through different plots: ***corner***, ***violin***, ***ridgeline***.

For example, these are possible applications.
#### Comparing multidimensional posterior distributions from different models
<img src="https://github.com/vascogennari/namib/assets/62053184/73deb9b2-a7a8-4ea0-adfd-df3f25d82c24" alt="drawing" width="400"/>

#### Comparing different pipelines on a set of GW events
<img src="https://github.com/vascogennari/namib/assets/62053184/0753647d-03e0-4592-bfd6-5670b8614656" alt="drawing" width="750"/>

#### Model selection between two competing models at different starting times
<img src="https://github.com/vascogennari/namib/assets/62053184/ac453ba8-150f-463e-a091-dbdc24218bd5" alt="drawing" width="750"/>

To run `namib`, you only need a set of samples (input files) and a config file.
- The name of the input files needs to contain the values of the *fields* following the structure:<br> `event_pipeline_model_submodel_time_GR-tag.ext`
- Choose a set of parameters and one *field* to compare the different files.
- Select the types of plots and customize bounds, dimensions and colors of the figures.
- Enjoy the results.

### Guide
You can find detailed instructions on `namib` usage in the [`namib` guide](https://github.com/vascogennari/namib/files/13635307/namib_guide.pdf).


## Installing and running the code
To install `namib`, follow the instructions:

    git clone https://github.com/vascogennari/namib.git
    cd namib
    conda install -c conda-forge hdf5 gwsurrogate
    pip install -r requirements.txt
    python setup.py install
    
Once the package is installed, you can easily launch `namib` as:

    namib --config-file config_path/config_filename.ini

To get familiar with `namib`, try to replicate the above plots by launching the config files in the `config_files/examples` directory. Remember to change the `~namib_path/` for the `samples` and `output` keys within the config files.
