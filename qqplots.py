import numpy as np
import statsmodels.api as sm
import pylab as py
import iris
import matplotlib.pyplot as plt


global_cube = iris.load_cube('/data/users/mgrant/ukcp/droughts/data/global/ACCESS1-3/TS1/monthly/dsi-3/ACCESS1-3_198102-200011_dsi-3.nc')
obs_60km_cube = iris.load_cube('/data/users/mgrant/ukcp/droughts/data/obs/60km/dsi-3/1980-2000/obs_60km_198102-200011_dsi-3.nc')

global_data_1d = global_cube.data.flatten()
non_mask_index = np.where(~global_data_1d.mask)
global_data_1d = global_data_1d[non_mask_index]
non_zero_index = np.where(global_data_1d != 0)
global_data_1d = global_data_1d[non_zero_index]

obs_60km_data_1d = obs_60km_cube.data.flatten()
non_mask_index = np.where(~obs_60km_data_1d.mask)
obs_60km_data_1d = obs_60km_data_1d[non_mask_index]
non_zero_index = np.where(obs_60km_data_1d != 0)
obs_60km_data_1d = obs_60km_data_1d[non_zero_index]

# Save the array into a file so we can plot it in R

global_data_1d.sort()
obs_60km_data_1d.sort()
max_global = max(global_data_1d)
min_global = min(global_data_1d)
max_obs = max(obs_60km_data_1d)
min_obs = min(obs_60km_data_1d)


plt.scatter(obs_60km_data_1d, global_data_1d)
plt.plot([min_obs, max_obs], [min_obs, max_obs], color='red')
plt.xlabel("Observation Quantiles")
plt.ylabel("Global Model's Quantiles")
plt.show()
# sm.qqplot(global_data_1d, dist=obs_60km_data_1d, line='45')
# py.show()
