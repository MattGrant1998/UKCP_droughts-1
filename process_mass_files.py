import os.path
import iris
import pre_processing

# runid = pre_processing.RUNID['MRI-CGCM3']['TS1']
# print(runid[3:])

def load_total_precip(model, time_slice, year, month):
    runid = pre_processing_model.LOCAL_RUNID[model][time_slice]

    filepath = f'/scratch/mgrant/UKCP/local/mass_pp_files/{time_slice}/{runid}/'
    filename = f'{runid[3:]}a.pm{year}{month}.pp'

    cube_rainfall = iris.load_cube(filepath + filename, 'stratiform_rainfall_flux')
    cube_snowfall = iris.load_cube(filepath + filename, 'stratiform_snowfall_flux')

    cube_precip = cube_rainfall + cube_snowfall

    return cube_precip

def merge_monthly_cubes_to_years(model, time_slice, year):
    runid = pre_processing_model.LOCAL_RUNID[model][time_slice]
    months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov']

    precip_cubes = iris.cube.CubeList()
    dec_cube = load_total_precip(model, time_slice, year - 1, 'dec')
    precip_cubes.append(dec_cube)

    for month in months:
        precip_cube = load_total_precip(model, time_slice, year, month)
        precip_cubes.append(precip_cube)

    year_cube = precip_cubes.merge_cube()
    filepath = f'/scratch/mgrant/UKCP/local/hires_rcm/UKCP18/cmip5_downscale' \
               f'/cpm_output/{time_slice}/monthly/{runid}/pr/'
    filename = f'{runid}_{year - 1}12-{year}11_pr.nc'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    # Convert units from 'kg m-2 s-1' to 'mm (day)-1'
    year_cube = year_cube * 86400

    year_cube.units = 'mm (day)-1'
    year_cube.var_name = 'precipitation_flux'

    iris.save(year_cube, filepath + filename)

def merge_all_monthly_cubes_to_years(time_slice):
    years = pre_processing_model.SLICES[time_slice]
    for model in pre_processing_model.CMIP_MODELS:
        for year in years:
            merge_monthly_cubes_to_years(model, time_slice, year)
