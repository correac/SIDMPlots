SIDMPlots
=========

Python package that produces various plots to analyse SIDM impact in structure formation.

Requirements
------------

The morpholopy package requires:

+ `python3.6` or above
+ see requirements.txt

Installing
----------

To get started using the package you need to set up a python virtual environment. The steps are as follows:

Clone SIDMPlots
```
git clone https://github.com/correac/SIDMPlots.git

cd SIDMPlots

python3 -m venv sidmplots_env
```

Now activate the virtual environment.

```
source sidmplots_env/bin/activate
```

Update pip just in case
```
pip install pip --upgrade

pip install -r requirements.txt
```

How to use it
-------------

See run.sh as an example. There you can specify in the folder where your simulation is stored (-d), the simulation
snapshot (e.g. -s), subhalo catalogue (with -c), simulation name (with -n) and the output folder where plots will be stored (-o).

Then type 
```
bash run.sh
```


