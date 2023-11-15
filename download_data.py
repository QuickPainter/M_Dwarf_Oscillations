import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy
from astropy.io import fits
import os
from astropy.timeseries import LombScargle
import glob
import pathlib
from astroquery.mast import Observations
from astropy import units as u
import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = 'white'
from astropy.coordinates import SkyCoord
from tqdm import tqdm
import pickle
import lightkurve as lk
import scipy.stats as stats
import traceback
from scipy.signal import find_peaks
import lightkurve



def download_data():

    with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_names.pkl','rb') as f:
        coords = pickle.load(f)


    paths = []
    for i in tqdm(range(1,5)):
        try:
            star_paths = []
            try:
                print(coords[i])
                search_result = lightkurve.search_lightcurve(coords[i],exptime='fast',mission='TESS')
                if len(search_result) > 0:
                    print('Data Found')
                    data = search_result.download_all(quality_bitmask="default", flux_column='pdcsap_flux')
                    for obs in range(len(data)):
                        dilution = data[obs].meta["CROWDSAP"]
                        data[obs].to_fits(path=f'lightkurve_data/heavy_data/{i}_{obs}_lc.fits', overwrite=True,CROWDSAP=dilution)
                        star_paths.append(f'{i}_{obs}_lc.fits')
                else:
                    star_paths.append('')
            except:
                search_result = lightkurve.search_lightcurve(coords[i],exptime='fast',mission='TESS')
                search_result.table["dataURL"] = search_result.table["dataURI"]  # workaround MAST issue
                data = search_result.download_all(quality_bitmask="default", flux_column='pdcsap_flux')
                if len(data)>0:
                    print('Data Found')
                    for obs in range(len(data)):
                        dilution = data[obs].meta["CROWDSAP"]
                        data[obs].to_fits(path=f'lightkurve_data/heavy_data/{i}_{obs}_lc.fits', overwrite=True,CROWDSAP=dilution)
                        star_paths.append(f'{i}_{obs}_lc.fits')
                else:
                    star_paths.append('')
            paths.append(star_paths)
        except:
            print(f"Error on star {i}")

    
    # with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data_paths.pkl','wb') as f:
    #     pickle.dump(paths,f)

if __name__ == '__main__':
    download_data()
