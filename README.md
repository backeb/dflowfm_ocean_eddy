# D-FlowFM_baroclinic_eddy_test

## Purpose: 
Notebook to analyse output from idealised D-Flow FM simulations of a baroclinic eddy.

## Installation:
The Notebook makes use of [dfm_tools](https://github.com/openearth/dfm_tools).

Follow the below steps to set up the Notebook:
`conda create --name dfm_tools_env -c conda-forge python=3.7 git spyder -y`
`conda activate dfm_tools_env`
`python -m pip install git+https://github.com/openearth/dfm_tools.git`
`conda install -c conda-forge "shapely>=1.7.0" -y`
`conda install -c conda-forge cartopy -y`
`conda install -c conda-forge geopandas -y`
`conda install -c conda-forge notebook`
`conda install -c anaconda ipykernel`
`python -m ipykernel install --user --name=dfm_tools_env`
