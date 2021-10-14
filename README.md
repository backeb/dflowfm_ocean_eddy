# dflowfm_ocean_eddy

## Purpose: 
Scripts and Jupyter Notebooks to analyse output from idealised D-Flow FM simulations of a baroclinic eddy.

## Installation:
The scripts and notebooks make use of [dfm_tools](https://github.com/openearth/dfm_tools).

To use dfm_tools in Spyder, follow the [dfm_tools installation instructions](https://github.com/openearth/dfm_tools#installation). 

Follow the below steps to set up dfm_tools to use in a Jupyter Notebooks:
```
conda create --name dfm_tools_env -c conda-forge python=3.7 git spyder -y
conda activate dfm_tools_env
python -m pip install git+https://github.com/openearth/dfm_tools.git
conda install -c conda-forge "shapely>=1.7.0" -y
conda install -c conda-forge cartopy -y
conda install -c conda-forge geopandas -y
conda install -c conda-forge notebook
conda install -c anaconda ipykernel
python -m ipykernel install --user --name=dfm_tools_env
```

Start the Notebook from the activated dfm_tools_env, e.g.
```
conda activate dfm_tools_env
jupyter-notebook BaroclinicVortex.ipynb
```
