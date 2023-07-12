import iris
import numpy as np
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import os

# Merge years from 1980 to 2000 together
obs_1980_to_2000_list = iris.cube.CubeList()
for year in range(1980, 2001):
    cube_year = iris.load_cube(f'/project/ukcp18/ncic_observations/post_processed/badc/ukmo-hadobs/'
                               f'data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.0.0/1km/rainfall/mon/'
                               f'v20181126/rainfall_hadukgrid_uk_1km_mon_{year}01-{year}12.nc')

    obs_1980_to_2000_list.append(cube_year)
    iris.util.equalise_attributes(obs_1980_to_2000_list)

obs_1980_to_2000 = obs_1980_to_2000_list.concatenate_cube()

# Slice cube to be from Dec 1980 to November 2000
obs_198012_to_200011 = obs_1980_to_2000[11:-1]

# Regrid cube to 2.2km grid to be the same as the local model's grid
local_grid = iris.load_cube('/scratch/mgrant/UKCP/local/full_time_slice/TS3/'
                            'monthly/mi-bd784/pr/mi-bd784_206012-208011_pr.nc')

regrid_scheme = iris.analysis.Linear(extrapolation_mode='mask')
obs_198012_to_200011_LocalGrid = obs_198012_to_200011.regrid(local_grid, regrid_scheme)

filepath_1km = '/scratch/mgrant/UKCP/obs/1km/rainfall/'
filename_1km = 'obs_198012-200011_rainfall_1km.nc'

if not os.path.exists(filepath_1km):
    os.makedirs(filepath_1km)

filepath_2km = '/scratch/mgrant/UKCP/obs/2.2km/rainfall/'
filename_2km = 'obs_198012-200011_rainfall_2.2km.nc'

if not os.path.exists(filepath_2km):
    os.makedirs(filepath_2km)

iris.save(obs_198012_to_200011, filepath_1km + filename_1km)
iris.save(obs_198012_to_200011_LocalGrid, filepath_2km + filename_2km)
