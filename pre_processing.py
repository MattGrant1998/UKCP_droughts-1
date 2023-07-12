import iris
import os
import numpy as np
import matplotlib.pyplot as plt
import iris.quickplot as qplt
# import process_mass_files

CMIP_MODELS = [
    'MRI-CGCM3',
    'MPI-ESM-LR',
    'ACCESS1-3',
    'IPSL-CM5A-MR'
]


MEANS = [
    'CMIP_multimodel_mean',
    'PPE_multimodel_mean',
    'full_multimodel_mean'
]

MEDIANS = [
    'CMIP_multimodel_median',
    'PPE_multimodel_median',
    'full_multimodel_median'
]

PPE_MEMBERS = [
    1,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    15
]

MEMBERS = {
    'MRI-CGCM3': 27,
    'MPI-ESM-LR': 26, # this is actually MPI-ESM-MR in UKCP Global - have to check if this is an issue/what difference is
    'ACCESS1-3': 23,
    'IPSL-CM5A-MR': 25,
}

SLICES = {
    'TS1': range(1981, 2001),
    # 'TS2': range(2021, 2041),
    'TS3': range(2061, 2081),
}

ALT_TIMESLICES = [
    [1980, 2000],
    [1970, 2000],
    [1980, 2010],
    [1970, 2010],
    [1960, 2010]
]

ENSEMBLE = {
    # 'full': PPE_MEMBERS,
    'full': PPE_MEMBERS + CMIP_MODELS,
    'PPE': PPE_MEMBERS,
    'CMIP': CMIP_MODELS
}

LOCAL_RUNID = {
    'MRI-CGCM3': {
        'TS1': 'mi-bd773',
        'TS2': 'mi-bd774',
        'TS3': 'mi-bd775',
    },
    'MPI-ESM-LR': {
        'TS1': 'mi-bd776',
        'TS2': 'mi-bd777',
        'TS3': 'mi-bd778',
    },
    'ACCESS1-3': {
        'TS1': 'mi-bd779',
        'TS2': 'mi-bd780',
        'TS3': 'mi-bd781',
    },
    'IPSL-CM5A-MR': {
        'TS1': 'mi-bd782',
        'TS2': 'mi-bd783',
        'TS3': 'mi-bd784',
    }
}

REGIONAL_RUNID = {
    'MRI-CGCM3': {
        'TS1': 'u-co311',
        'TS2': 'u-co312',
        'TS3': 'u-co314'
    },
    'MPI-ESM-LR': {
        'TS1': 'u-co315',
        'TS2': 'u-co316',
        'TS3': 'u-co317'
    },
    'ACCESS1-3': {
        'TS1': 'u-co318',
        'TS2': 'u-co319',
        'TS3': 'u-co320'
    },
    'IPSL-CM5A-MR': {
        'TS1': 'u-co321',
        'TS2': 'u-co322',
        'TS3': 'u-co323'
    }
}

COORDS = {
    'lat_long': ('grid_latitude', 'grid_longitude'),
    'x_y': ('projection_y_coordinate', 'projection_x_coordinate')
}

COUNTRIES = [
    'uk',
    'Scotland',
    'England',
    'Wales',
    'Northern_Ireland'
]


VARIABLES = ['pr', ]

CALENDAR = {
    'MRI-CGCM3': 'gregorian',
    'MPI-ESM-LR': 'gregorian',
    'ACCESS1-3': 'gregorian',
    'IPSL-CM5A-MR': '365day',
}

DOMAIN = [
    'global',
    'regional',
    'local'
]

RESOLUTION = {
    'global': '60km',
    'regional': '12km',
    'local': '2.2km'
}

DROUGHT_INDICES = [
    'dsi',
]

OBS_RES = [
    '1km',
    '12km',
    '60km'
]

DOM_OBS_RES = {
    'global': '60km',
    'regional': '12km',
    'local': '1km'
}

COORD_LIMS = {
    'y_x': [(-87500, 1220500), (-99500, 699500)],
    'lat_long': [(-3.5, 8.08), (355.8, 363.2)]
}

datadir_unprocessed = '/data/users/mgrant/ukcp/droughts/unprocessed_data'
data_files = '/data/users/mgrant/ukcp/droughts/data'
plots_files = '/scratch/mgrant/UKCP/droughts/plots'
stat_tests_files = '/data/users/mgrant/ukcp/droughts/stat_tests'
cmip5_ppe_data = ''
ppe_cpm_data = '/project/ukcp/land-cpm/uk/2.2km/rcp85'
ppe_rcm_data = '/project/ukcp/land-rcm/uk/12km/rcp85'


def check_leap_year(year):
    if (year % 4) == 0:
        if (year % 100) == 0:
            if (year % 400) == 0:
                return True
            else:
                return False
        else:
            return True
    else:
        return False


def apply_landsea_mask(cube, lsm, domain):
    """
    Applies the land sea mask to the data - designed to work with the UK so may not work with other
    masks.
    Args:
         cube (iris.cube.Cube): cube that needs to be masked
         lsm (iris.cube.Cube): cube with required mask on data
    Returns:
        masked_cube (iris.cube.Cube): input cube but with the sea masked out
    """
    if domain == 'global':
        # Convert the data in lsm from 1/0 to True/False, where True is over the sea
        lsm_boolean = np.where(lsm.data == 1, False, True)
        # Apply mask to each time step of the cube
        time_length = len(cube.coord('time').points)
        masked_cube = cube.copy()
        for t in range(time_length):
            masked_cube.data[t, ...] = iris.util.mask_cube(cube[t, ...], lsm_boolean).data

    elif domain == 'local':
        lats = cube.coord('grid_latitude').points
        longs = cube.coord('grid_longitude').points
        cube_lat_constraint = iris.Constraint(grid_latitude=lambda cell: lats[0] <= cell <= lats[-1])
        cube_long_constraint = iris.Constraint(grid_longitude=lambda cell: longs[0] <= cell <= longs[-1])
        lsm = lsm.extract(cube_lat_constraint & cube_long_constraint)
        lsm_data = lsm.data.mask
        masked_cube = iris.util.mask_cube(cube, lsm_data)

    elif domain == 'regional':
        lats = lsm.coord('projection_x_coordinate').points
        longs = lsm.coord('projection_y_coordinate').points
        lsm_lat_constraint = iris.Constraint(projection_x_coordinate=lambda cell: lats[0] <= cell <= lats[-1])
        lsm_long_constraint = iris.Constraint(projection_y_coordinate=lambda cell: longs[0] <= cell <= longs[-1])
        cube = cube.extract(lsm_lat_constraint & lsm_long_constraint)
        lsm_data = lsm.data.mask
        masked_cube = iris.util.mask_cube(cube, lsm_data)

    return masked_cube



def regrid_regional_to_OSGB(regional_cube):
    # Load obs cube to regrid the CMIP5 cube to
    ppe_regional_cube = iris.load_cube(f'{data_files}/uk/obs/12km/rainfall/'
                                       f'1980-2000/obs_12km_198012-200011_rainfall.nc')
    # Regrid cube to PPE grid
    regrid_scheme = iris.analysis.Linear()
    regional_cube_regrid = regional_cube.regrid(ppe_regional_cube, regrid_scheme)

    return regional_cube_regrid


def spatially_constrain_cube(cube):
    if 'projection_x_coordinate' in [coord.name() for coord in cube.dim_coords]:
        coord_lims = COORD_LIMS['y_x']

        y_min = coord_lims[0][0]
        y_max = coord_lims[0][1]
        x_min = coord_lims[1][0]
        x_max = coord_lims[1][1]

        x_constraint = iris.Constraint(projection_x_coordinate=lambda cell: x_min <= cell <= x_max)
        y_constraint = iris.Constraint(projection_y_coordinate=lambda cell: y_min <= cell <= y_max)

    else:
        coord_lims = COORD_LIMS['lat_long']

        lat_min = coord_lims[0][0]
        lat_max = coord_lims[0][1]
        long_min = coord_lims[1][0]
        long_max = coord_lims[1][1]

        x_constraint = iris.Constraint(grid_longitude=lambda cell: long_min <= cell <= long_max)
        y_constraint = iris.Constraint(grid_latitude=lambda cell: lat_min <= cell <= lat_max)

    constrained_cube = cube.extract(x_constraint & y_constraint)

    return constrained_cube


# def create_monthly_regional_TS2(model, variable):

# obs_cube = iris.load_cube('/data/users/mgrant/ukcp/droughts/data/obs/12km/dsi-3/1980-2000/obs_12km_198102-200011_dsi-3.nc')
# obs_const = spatially_constrain_cube(obs_cube)
# print(obs_const)
def load_local_or_regional_data(model, time_slice, variable, year, PPE=False, domain='local'):
    """
    Basic function to simply load the cube of interest for a 1-year timeslice - unless PPE=True and domain='regional'
    then it loads in the cube of the full timeslice.
    Args:
        model (str): name of the CMIP5 model wanting to be loaded - available ones in above MODELS dictionary
        time_slice (str): selection of timeslice - either TS1 (1980-2000), TS2 (2020-2040), or TS3 (2060-2080)
        variable (str): variable of interest
        year (int): specified year of interest (must be within chosen timeslice)
    """
    if domain == 'local':
        # filepath_in = f'/scratch/mgrant/UKCP/local/unprocessed_data/hires_rcm/UKCP18/cmip5_downscale/cpm_output'
        if PPE:
            model = str(model).zfill(2)
            filepath_in = f'{ppe_cpm_data}/{model}/{variable}/mon/v20210615/'
            filename_in = f'{variable}_rcp85_land-cpm_uk_2.2km_{model}_mon_{year - 1}12-{year}11.nc'
        else:
            runid = LOCAL_RUNID[model][time_slice]
            if time_slice == 'TS2':
                filepath_in = f'{datadir_unprocessed}/local/hires_rcm/UKCP18/cmip5_downscale' \
                              f'/cpm_output/{time_slice}/monthly/{runid}/{variable}/'
            else:
                filepath_in = f'/project/hires_rcm/UKCP18/cmip5_downscale' \
                              f'/cpm_output/{time_slice}/monthly/{runid}/{variable}/'

            filename_in = f'{runid}_{year - 1}12-{year}11_{variable}.nc'
    elif domain == 'regional':
        runid = REGIONAL_RUNID[model][time_slice]
        filepath_in = f'/project/hires_rcm/cshort/cmip5/{time_slice}/{runid}/monthly/r001i1p00000/{variable}/'
        filename_in = f'r001i1p00000_{year - 1}12-{year}11_{variable}.nc'

    # if not os.path.exists(filepath_in):
    #     process_mass_files.merge_monthly_cubes_to_years(model, time_slice, variable, year)

    # Load cube of mm/day precipitation averaged over the month
    cube_mmperday = iris.load_cube(filepath_in + filename_in)

    # Mask the cube
    # lsm = iris.load_cube('/project/ukcp/extra/lsm_land-cpm_BI_2.2km.nc')
    # masked_cube_mmperday = apply_landsea_mask(cube_mmperday, lsm, 'local')

    # Multiply mm/day cube by number of days in month to get mm/month data
    cube_mmpermonth = cube_mmperday.copy()
    if PPE:
        cube_mmpermonth.data = cube_mmperday.data * 30
        cube_mmpermonth = cube_mmpermonth[0, ...]
    else:
        days_in_month = np.array([31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30])
        if (CALENDAR[model] == 'gregorian') & (check_leap_year(year) == True):
            days_in_month[2] = 29

        for i in range(12):
            cube_mmpermonth.data[i, ...] = cube_mmperday.data[i, ...] * days_in_month[i]

    cube_mmpermonth.units = 'mm (month)-1'

    return cube_mmpermonth


def create_year_coord(year_range):
    number_of_years = len(year_range)
    start_year = year_range[0]

    years = np.zeros(number_of_years * 12)
    for i in range(number_of_years):
        years[0] = start_year - 1
        years[1 + (12 * i):13 + (12 * i)] = year_range[i]

    years_coord = iris.coords.AuxCoord(years, standard_name=None, long_name='year', var_name='year',
                                       units='1', bounds=None, attributes=None, coord_system=None)

    return years_coord


def create_month_number_coord(year_range):
    number_of_years = len(year_range)
    month_number_one_year = np.roll(np.arange(1, 13), 1)
    month_number = np.tile(month_number_one_year, number_of_years)

    month_number_coord = iris.coords.AuxCoord(month_number, standard_name=None, units='1',
                                              long_name='month_number', var_name='month_number')

    return month_number_coord


def create_full_timeslice_of_local_and_CMIP5_regional_cube(model, time_slice, variable, PPE=False, domain='local'):
    # Initialise CubeList
    yearly_cubes = iris.cube.CubeList()

    if PPE:
        model = str(model).zfill(2)

    # Load cube from each year in time slice and append onto yearly_cubes CubeList
    for year in SLICES[time_slice]:
        cube_year_i = load_local_or_regional_data(model, time_slice, variable, year, PPE, domain)
        yearly_cubes.append(cube_year_i)

    # Merge all yearly cubes together
    iris.util.equalise_attributes(yearly_cubes)
    full_timeslice_cube = yearly_cubes.concatenate_cube()

    if domain == 'regional':
        full_timeslice_cube = regrid_regional_to_OSGB(full_timeslice_cube)

    # if UK_mask:
    res = RESOLUTION[domain]
    UK_mask = iris.load_cube(f'{data_files}/uk/obs/{res}/rainfall/'
                             f'1980-2000/obs_{res}_198012-200011_rainfall.nc')

    full_timeslice_cube = spatially_constrain_cube(full_timeslice_cube)
    UK_mask = spatially_constrain_cube(UK_mask)
    # lats = UK_mask.coord('grid_latitude').points
    # longs = UK_mask.coord('grid_longitude').points
    # lsm_lat_constraint = iris.Constraint(grid_latitude=lambda cell: lats[0] <= cell <= lats[-1])
    # lsm_long_constraint = iris.Constraint(grid_longitude=lambda cell: longs[0] <= cell <= longs[-1])
    # full_timeslice_cube = full_timeslice_cube.extract(lsm_lat_constraint & lsm_long_constraint)
    # UK_mask = UK_mask_cube.data.mask
    # qplt.pcolormesh(UK_mask[0, ...])
    # plt.show()
    if PPE and domain == 'local':
        full_timeslice_cube = apply_landsea_mask(full_timeslice_cube, UK_mask, 'local')
    elif domain == 'regional':
        full_timeslice_cube = apply_landsea_mask(full_timeslice_cube, UK_mask, 'regional')
    else:
        mask = UK_mask.data.mask
        full_timeslice_cube = iris.util.mask_cube(full_timeslice_cube, mask)

    full_timeslice_cube = spatially_constrain_cube(full_timeslice_cube)

    filepath = f'{data_files}/uk/{domain}/{model}/{time_slice}/monthly/{variable}/'
    # qplt.pcolormesh(full_timeslice_cube[0,...])
    # plt.show()
    # else:
    #     filepath = f'/scratch/mgrant/UKCP/local/data/{time_slice}/monthly/{runid}/{variable}/LandSeaMask/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    year_range = SLICES[time_slice]
    start_year = year_range[0]
    end_year = year_range[-1]
    filename = f'{model}_{start_year - 1}12-{end_year}11_{variable}.nc'

    if not PPE:
        years_coord = create_year_coord(year_range)
        full_timeslice_cube.add_aux_coord(years_coord, 0)

        month_number_coord = create_month_number_coord(year_range)
        full_timeslice_cube.add_aux_coord(month_number_coord, 0)

    iris.save(full_timeslice_cube, filepath + filename)


# create_full_timeslice_of_local_cube(1, 'TS1', 'pr', PPE=True)


def pre_process_all_global_and_PPE_regional_timeslices(model, timeslice, variable, PPE=False, domain='global'):
    if PPE or domain == 'regional':
        member = str(model).zfill(2)
    else:
        member = MEMBERS[model]

    years = SLICES[timeslice]
    firstyear = years[0] - 1
    lastyear = years[-1]

    if domain == 'regional':
        filepath_in = f'{ppe_rcm_data}/{member}/{variable}/mon/v20190731/'
        filename_in = f'{variable}_rcp85_land-rcm_uk_12km_{member}_mon_198012-208011.nc'
        cube = iris.load_cube(filepath_in + filename_in)
        cube = cube[0, ...]  # Remove ensemble_member coord
    else:
        filepath_in = f'/project/ukcp18/data/land-gcm/uk/60km/rcp85/{member}/{variable}/mon/latest/'
        filename_in = f'{variable}_rcp85_land-gcm_uk_60km_{member}_mon_189912-209911.nc'
        cube = iris.load_cube(filepath_in + filename_in)
        cube = cube[0, ...]  # Remove ensemble_member coord

    # Slice the cube to be from Dec first year in timeslice to Nov of last year in timeslice
    year_constraint = iris.Constraint(year=lambda cell: firstyear <= cell <= lastyear)
    cube_yearsliced = cube.extract(year_constraint)
    if domain == 'regional' and timeslice == 'TS1':
        cube_timesliced = cube_yearsliced[:-1, ...]
    elif domain == 'regional' and timeslice == 'TS3':
        cube_timesliced = cube_yearsliced[11:, ...]
    else:
        cube_timesliced = cube_yearsliced[11:-1, ...]

    # Apply mask to cube
    # if UK_mask:
    if domain == 'regional':
        cube_timesliced = regrid_regional_to_OSGB(cube_timesliced)
        UK_mask = iris.load_cube(f'{data_files}/uk/obs/12km/rainfall/1980-2000/obs_12km_198012-200011_rainfall.nc')
        masked_cube = apply_landsea_mask(cube_timesliced, UK_mask, 'regional')
    else:
        UK_mask = iris.load_cube('/project/ukcp/extra/lsm_land-gcm_uk_60km.nc')
        masked_cube = apply_landsea_mask(cube_timesliced, UK_mask, 'global')
    # UK_mask = UK_mask_cube.data.mask

    masked_cube = spatially_constrain_cube(masked_cube)
    if PPE or domain == 'regional':
        filepath = f'{data_files}/uk/{domain}/{member}/{timeslice}/monthly/{variable}/'
        filename = f'{member}_{firstyear}12-{lastyear}11_{variable}.nc'
    else:
        filepath = f'{data_files}/uk/{domain}/{model}/{timeslice}/monthly/{variable}/'
        filename = f'{model}_{firstyear}12-{lastyear}11_{variable}.nc'

    # else:
    #     lsm = iris.load_cube('/project/ukcp/extra/lsm_land-gcm_uk_60km.nc')
    #     masked_cube = apply_landsea_mask(cube_timesliced, lsm, 'global')
    #     filepath = f'/scratch/mgrant/UKCP/global/data/{timeslice}/monthly/{member}/{variable}/LandSeaMask/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    masked_cube.remove_coord('yyyymm')

    iris.save(masked_cube, filepath + filename)

# pre_process_all_global_and_PPE_regional_timeslices(1, 'TS1', 'pr', PPE=True, domain='global')
def pre_process_obs_full_timeslice(variable, resolution):
    # variable name for precipitation differs between obs and model so set to correct for obs
    if variable == 'pr':
        obs_variable = 'rainfall'

    obs_cubelist = iris.load(f'/project/ukcp18/ncic_observations/post_processed/badc/ukmo-hadobs/'
                             f'data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.0.0/{resolution}/{obs_variable}/mon/'
                             f'v20181126/*')

    iris.util.equalise_attributes(obs_cubelist)
    obs_cube = obs_cubelist.concatenate_cube()

    obs_cube = spatially_constrain_cube(obs_cube)

    year_range = range(1862, 2018)
    years_coord = create_year_coord(year_range)
    obs_cube.add_aux_coord(years_coord, 0)

    filepath_1km = f'{data_files}/uk/obs/{resolution}/{obs_variable}/full_timeslice/'
    filename_1km = f'obs_{resolution}_194001-201712_{obs_variable}.nc'

    if not os.path.exists(filepath_1km):
        os.makedirs(filepath_1km)

    iris.save(obs_cube, filepath_1km + filename_1km)

# pre_process_obs_full_timeslice('pr', '1km')
def pre_process_timesliced_obs(resolution, firstyear, lastyear):
    obs_full_timeslice = iris.load_cube(f'{data_files}/uk/obs/{resolution}/rainfall/full_timeslice/'
                                        f'obs_{resolution}_194001-201712_rainfall.nc')

    timeslice_constraint = iris.Constraint(year=lambda cell: firstyear <= cell <= lastyear)
    cube_yearsliced = obs_full_timeslice.extract(timeslice_constraint)
    cube_timesliced = cube_yearsliced[11:-1, ...]

    filepath = f'{data_files}/uk/obs/{resolution}/rainfall/{firstyear}-{lastyear}/'
    filename = f'obs_{resolution}_{firstyear}12-{lastyear}11_rainfall.nc'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    iris.save(cube_timesliced, filepath + filename)


def create_2p2km_obs(variable, firstyear, lastyear):
    if variable[:-2] == 'dsi':
        firstmonth = (str(int(variable[-1]) - 1)).zfill(2)
        file_firstyear = firstyear + 1
    else:
        firstmonth = 12
        file_firstyear = firstyear

    obs_1km = iris.load_cube(f'{data_files}/uk/obs/1km/{variable}/{firstyear}-{lastyear}/'
                             f'obs_1km_{file_firstyear}{firstmonth}-{lastyear}11_{variable}.nc')

    local_grid = iris.load_cube(f'/project/hires_rcm/UKCP18/cmip5_downscale/cpm_output/'
                                f'TS1/monthly/mi-bd773/pr/mi-bd773_198012-198111_pr.nc')

    regrid_scheme = iris.analysis.Linear(extrapolation_mode='mask')
    obs_2p2km = obs_1km.regrid(local_grid, regrid_scheme)
    filepath = f'{data_files}/uk/obs/2.2km/{variable}/{firstyear}-{lastyear}/'
    filename = f'obs_2.2km_{firstyear}{firstmonth}-{lastyear}11_{variable}.nc'

    obs_2p2km = spatially_constrain_cube(obs_2p2km)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    iris.save(obs_2p2km, filepath + filename)


# def pre_process_obs(variable, resolution):
#     # variable name for precipitation differs between obs and model so set to correct for obs
#     if variable == 'pr':
#         obs_variable = 'rainfall'
#
#     # Merge years from 1980 to 2000 together
#     obs_1980_to_2000_list = iris.cube.CubeList()
#     for year in range(1980, 2001):
#         cube_year = iris.load_cube(f'/project/ukcp18/ncic_observations/post_processed/badc/ukmo-hadobs/'
#                                    f'data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.0.0/{resolution}/{obs_variable}/mon/'
#                                    f'v20181126/rainfall_hadukgrid_uk_{resolution}_mon_{year}01-{year}12.nc')
#
#         obs_1980_to_2000_list.append(cube_year)
#         iris.util.equalise_attributes(obs_1980_to_2000_list)
#
#     obs_1980_to_2000 = obs_1980_to_2000_list.concatenate_cube()
#
#     # Slice cube to be from Dec 1980 to November 2000
#     obs_198012_to_200011 = obs_1980_to_2000[11:-1]
#
#     # Regrid cube to 2.2km grid to be the same as the local model's grid
#     local_grid = iris.load_cube(f'/project/hires_rcm/UKCP18/cmip5_downscale/cpm_output/TS1/monthly/'
#                                 f'mi-bd773/{variable}/mi-bd773_198012-198111_{variable}.nc')
#
#     global_grid = iris.load_cube(f'/project/ukcp18/data/land-gcm/uk/60km/rcp85/27/{variable}/mon/'
#                                  f'latest/{variable}_rcp85_land-gcm_uk_60km_27_mon_189912-209911.nc')
#
#     regrid_scheme = iris.analysis.Linear(extrapolation_mode='mask')
#     obs_198012_to_200011_LocalGrid = obs_198012_to_200011.regrid(local_grid, regrid_scheme)
#     obs_198012_to_200011_GlobalGrid = obs_198012_to_200011.regrid(global_grid, regrid_scheme)
#
#     filepath_1km = f'/scratch/mgrant/UKCP/data/obs/1km/{obs_variable}/'
#     filename_1km = f'obs_1km_198012-200011_{obs_variable}.nc'
#
#     if not os.path.exists(filepath_1km):
#         os.makedirs(filepath_1km)
#
#     filepath_2km = f'/scratch/mgrant/UKCP/data/obs/2.2km/{obs_variable}/'
#     filename_2km = f'obs_2.2km_198012-200011_{obs_variable}.nc'
#
#     if not os.path.exists(filepath_2km):
#         os.makedirs(filepath_2km)
#
#     filepath_60km = f'/scratch/mgrant/UKCP/data/obs/60km/{obs_variable}/'
#     filename_60km = f'obs_60km_198012-200011_{obs_variable}.nc'
#
#     if not os.path.exists(filepath_60km):
#         os.makedirs(filepath_60km)
#
#     year_range = range(1981, 2001)
#     years_coord = create_year_coord(year_range)
#     obs_198012_to_200011.add_aux_coord(years_coord, 0)
#     obs_198012_to_200011_LocalGrid.add_aux_coord(years_coord, 0)
#
#     iris.save(obs_198012_to_200011, filepath_1km + filename_1km)
#     iris.save(obs_198012_to_200011_LocalGrid, filepath_2km + filename_2km)
#     iris.save(obs_198012_to_200011_GlobalGrid, filepath_60km + filename_60km)

# pre_process_all_global_and_PPE_regional_timeslices('MRI-CGCM3', 'TS1', 'pr', PPE=False)
def main():
    print('Processing...')
    for var in VARIABLES:
        print('Variable: ', var)
        for timeslice in SLICES:
            print('Time Slice: ', timeslice)
            for model in CMIP_MODELS:
                print('Model: ', model)
                # create_full_timeslice_of_local_and_CMIP5_regional_cube(model, timeslice, var)
                # pre_process_all_global_and_PPE_regional_timeslices(model, timeslice, var, PPE=False)
                # create_full_timeslice_of_local_and_CMIP5_regional_cube(model, timeslice, var, domain='regional')
            for member in PPE_MEMBERS:
                print('Member: ', member)
                # create_full_timeslice_of_local_and_CMIP5_regional_cube(member, timeslice, var, PPE=True)
                # pre_process_all_global_and_PPE_regional_timeslices(member, timeslice, var, PPE=True)
                # pre_process_all_global_and_PPE_regional_timeslices(member, timeslice, var, domain='regional')
    for var in VARIABLES:
        for res in OBS_RES:
            pre_process_obs_full_timeslice(var, res)
    #
    #
    for res in OBS_RES:
        for slice in ALT_TIMESLICES:
            pre_process_timesliced_obs(res, slice[0], slice[1])
    #
    # for slice in ALT_TIMESLICES:
    # create_2p2km_obs('rainfall', 1980, 2000)

# create_full_timeslice_of_local_and_regional_cube(1, 'TS1', 'pr', PPE=True, domain='regional')
# pre_process_all_global_and_PPE_regional_timeslices(1, 'TS3', 'pr', domain='regional')

# pre_process_obs('pr')
# create_full_timeslice_of_local_cube('ACCESS1-3', 'TS1', 'pr')
# pre_process_global('MRI-CGCM3', 'TS1', 'pr')
#
if __name__ == '__main__':
    main()


# cube_mmperday = iris.load_cube(
#         f'/project/hires_rcm/UKCP18/cmip5_downscale/cpm_output/TS1/monthly'
#         f'/mi-bd773/pr/mi-bd773_198312-198411_pr.nc'
#     )
# days_in_month = np.array([31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30])
# print(days_in_month[2])
# for i in range(12):
#     cube_mmpermonth = cube_mmperday[i, ...] * days_in_month[i, ...]
