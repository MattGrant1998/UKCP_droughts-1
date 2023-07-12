import iris
import os
import dask.array as da
import pre_processing


def spi(model, time_slice, spi_number):
    """
    Calculates the SPI-x for the given model and timeslice, where x is the spi_number.

    Args:
        model (str): name of the model of interest (model names defined in pre_processing.MODELS)
        time_slice (str): time slice of interest (timeslice names stored in pre_processing.SLICES)
        spi_number (int): number of months for which SPI is calculated over
    """
    runid = pre_processing_model.LOCAL_RUNID[model][time_slice]
    filepath_in = f'/scratch/mgrant/UKCP/local/full_time_slice/{time_slice}/monthly/{runid}/pr/'

    if not os.path.exists(filepath_in):
        pre_processing_model.create_full_timeslice_cube(model, time_slice, 'pr')

    years = pre_processing_model.SLICES[time_slice]
    start_year = years[0]
    end_year = years[-1]
    filename_in = f'{runid}_{start_year - 1}12-{end_year}11_pr.nc'

    precip_cube = iris.load_cube(filepath_in + filename_in)

    rolling_average_precip = precip_cube.rolling_window(['time'], iris.analysis.MEAN, spi_number)
    std = rolling_average_precip.collapsed('time', iris.analysis.STD_DEV)
    mean = rolling_average_precip.collapsed('time', iris.analysis.MEAN)

    # Calculate SPI
    spi = (rolling_average_precip - mean) / std # should mean always be 1980-2000 average or vary with each timeslice??
    spi.data = da.ma.masked_invalid(spi.core_data())

    # Set some attributes
    spi.var_name = f'spi_{spi_number}month'
    spi.long_name = f'Standardised Precipitation Index ({spi_number} month)'
    spi.attributes['variable_id'] = 'spi'

    filepath_out = f'/scratch/mgrant/UKCP/local/full_time_slice/{time_slice}/monthly/{runid}/SPI-{spi_number}/'

    if not os.path.exists(filepath_out):
        os.makedirs(filepath_out)

    filename_out = f'{runid}_{start_year - 1}12-{end_year}11_SPI{spi_number}.nc'

    iris.save(spi, filepath_out + filename_out)

def calculate_all_spi(spi_number):
    for model in pre_processing_model.CMIP_MODELS:
        for time_slice in pre_processing_model.SLICES:
            spi(model, time_slice, spi_number)


def main():
    calculate_all_spi(3)

if __name__ == '__main__':
    main()
