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

def split_obs(lcs):
    result = []
    current_subarray = [lcs[0]]

    for i in range(1, len(lcs)):
        if abs(lcs.sector[i]-lcs.sector[i - 1]) < 3:
            current_subarray.append(lcs[i])
        else:
            result.append(current_subarray)
            current_subarray = [lcs[i]]

    result.append(current_subarray)
    result = [lightkurve.LightCurveCollection(i) for i in result]
    return result


def get_paths():
    cwd = os.getcwd()
    current_dir = pathlib.Path(cwd)
    all_paths = []
    for num in range(998):
        file_paths = [str(i) for i in list(current_dir.glob(f"lightkurve_data/heavy_data/{num}_*_lc.fits"))]
        all_paths.append(file_paths)
    return all_paths

def process(lc): 
    # flatten and remove flares
    clean_lc = lc.remove_nans().remove_outliers(3)
    ppm_unit = clean_lc.flux[0].unit
    flat_lc, trend_lc = clean_lc.flatten(window_length=101, polyorder=2, return_trend=True, break_tolerance=5, niters=3, sigma=3)
    # residuals = np.array(clean_lc.flux)-np.array(trend_lc.flux)
    # flat_lc['flux'] = residuals
    flat_lc = flat_lc.remove_outliers(3)

    flat_lc['flux'] = (flat_lc['flux']-1)*ppm_unit*10**6
    # get perio
    # frequency = np.arange(1, 4000,0.008)

    # fig, axs = plt.subplots(2,1,figsize=(20,10))
    # flat_lc.plot(ax=axs[0])
    # plt.show()    


    perio = flat_lc.to_periodogram(method='lombscargle',frequency=np.arange(2, 2160,0.008)).smooth(filter_width=2)
    threshold =  4*np.mean(np.array(perio.power))
    threshold_freq = [(i*10**6)/(24*60*60) for i in np.array(perio.smooth(filter_width=20).frequency)]
    # perio = perio.smooth(filter_width=2)
    perio_psd = flat_lc.to_periodogram(method='lombscargle',frequency=np.arange(20, 20000,0.008),normalization='psd').smooth(filter_width=2)

    # fig, axs = plt.subplots(2,1,figsize=(20,10))
    # perio_psd.plot(scale='log',ax=axs[0])
    # plt.show()

    # get thresholds

    frequency_hz = [(i*10**6)/(24*60*60) for i in np.array(perio.frequency)]
    # plt.plot(frequency_hz,np.array(perio.power),color='black')
    # plt.xlabel("Frequency [uHz]")
    # plt.ylabel("Amp. [ppm]")
    # plt.axhline(y=threshold,color='red')
    # plt.xlim(0,1000)
    # plt.show()
    # print('max',frequency_hz[np.argmax(np.array(perio.power))])
    return perio, perio_psd, threshold

def get_perio(lc_data):
    consecutive_lcs = split_obs(lc_data)

    perios = []
    perio_psds = []
    flats = []


    thresholds = []


    for lcs in consecutive_lcs:
        if len(lcs) > 1:
            lc = lcs.stitch().normalize(unit='ppm')
            perio, perio_psd,threshold = process(lc)
            ppm_unit = perio.power[0].unit
            psd_unit = perio_psd.power[0].unit
            perios.append(np.array(perio.power))
            perio_psds.append(np.array(perio_psd.power))
            thresholds.append(threshold)
        else:
            for lc in lcs:
                lc = lc.normalize(unit='ppm')
                perio, perio_psd,threshold = process(lc)
                ppm_unit = perio.power[0].unit
                psd_unit = perio_psd.power[0].unit
                thresholds.append(threshold)
                perios.append(np.array(perio.power))
                perio_psds.append(np.array(perio_psd.power))

    
    summed_perio = sum(perios)
    summed_perio_psd = sum(perio_psds)

    perio.power = np.array(summed_perio)*ppm_unit
    perio_psd.power = np.array(summed_perio_psd)*psd_unit
    return perio, perio_psd, np.min(thresholds)

def run_analysis():


    m_dwarf_sample = pd.read_csv('/Users/caleb/research/Astro_98/all_stars_mags_mass.csv',delimiter=',',index_col=0)
    # ra = np.array(m_dwarf_sample[m_dwarf_sample.columns[0]])
    # dec = np.array(m_dwarf_sample[m_dwarf_sample.columns[1]])
    names = np.array(m_dwarf_sample['SimbadName'])
    masses = [float(i) for i in np.array(m_dwarf_sample[["Mass"]])]
    comps = [i[0].strip() for i in np.array(m_dwarf_sample[["Comp"]])]

    # ras = [i.replace(' ',':') for i in ra]
    # decs = [i.replace(' ',':') for i in dec]

    # coords = []
    # for i in range(len(decs)):
    #     coords.append(ras[i]+" "+decs[i])

    all_paths = get_paths()

    powers = []
    power_psds = []
    flats = []
    thresholds = []
    dilutions = []
    stars_analyzed = []

    for num in tqdm(range(10,len(names))):
        # print("Star: ",num)
        # try:
            if len(all_paths[num]) != 0:
                if len(all_paths[num]) > 0:
                    lcs_list = []
                    for i in all_paths[num]:
                        lc = lightkurve.io.read(f'{i}')
                        lcs_list.append(lc)
                    lcs = lightkurve.LightCurveCollection(lcs_list)
                else:
                    lc = lightkurve.io.read(f'{all_paths[num][0]}')
                    lcs = lightkurve.LightCurveCollection(lc)



                # print(lcs)
                # get contamination
                dilution = lc.meta["CROWDSAP"]
                fig, axs = plt.subplots(2,1,figsize=(20,20))
                lc.remove_outliers().plot(ax=axs[0],color='black',alpha=.8)

                # sort lightkurves by sector
                sorted_indexes = sorted(range(len(lcs.sector)), key=lambda x: lcs.sector[x])
                lcs = lcs[sorted_indexes]
                perio, perio_psd,threshold = get_perio(lcs)

                # perio.plot(scale='log')
                # perio_psd.plot(scale='log')


                flat = perio_psd.flatten()
                flat.plot(ax=axs[1],color='black',alpha=.8)
                axs[0].tick_params(labelsize=20)
                axs[1].tick_params(labelsize=20)
                plt.rcParams.update({'font.size': 30})
                axs[0].set_xlabel("Time [d]",fontsize=25)
                axs[0].set_ylabel("Flux [e/s]",fontsize=25)
                axs[1].set_xlabel("Frequency [uHz]",fontsize=25)
                axs[1].set_ylabel("SNR",fontsize=25)
                axs[0].set_title(f'Star {num} -- {names[num]} {comps[num]} -- (Mass: {np.round(masses[num],2)}M)',fontsize=35)
                plt.tight_layout()
                plt.savefig(f'lightkurve_data/heavy_data/power_fig_plots/{num}_flat_power_fig.png')
                plt.close()

                # flat_perio = perio.smooth(filter_width=2)
                # flat_perio.plot(color='red')
                powers.append(np.array(perio.power))
                power_psds.append(np.array(perio_psd.power))
                thresholds.append(threshold)
                flats.append(np.array(flat.power))
                dilutions.append(dilution)
                stars_analyzed.append(num)


                with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/perio/{num}_perio.pkl','wb') as f:
                    pickle.dump(np.array(perio.power),f)

                with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/psd/{num}_psd.pkl','wb') as f:
                    pickle.dump(np.array(perio_psd.power),f)

                with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/flat/{num}_flat.pkl','wb') as f:
                    pickle.dump(np.array(flat.power),f)

                with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/thresholds/{num}_threshold.pkl','wb') as f:
                    pickle.dump(threshold,f)

            else:
                powers.append(np.array([]))
                power_psds.append(np.array([]))
                thresholds.append(-1)
                flats.append(np.array([]))
                dilutions.append(-1)
                stars_analyzed.append(num)
        # except:
        #     print("error on star: ",num)
        #     powers.append(np.array([]))
        #     power_psds.append(np.array([]))
        #     thresholds.append(0)
        #     flats.append(np.array([]))
        #     dilutions.append(0)
        #     stars_analyzed.append(f'X{num}X')




    with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/thresholds.pkl','wb') as f:
        pickle.dump(thresholds,f)

    with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/power_psds.pkl','wb') as f:
        pickle.dump(power_psds,f)

    with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/flats.pkl','wb') as f:
        pickle.dump(flats,f)

    with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/powers.pkl','wb') as f:
        pickle.dump(powers,f)

    with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/dilutions.pkl','wb') as f:
        pickle.dump(dilutions,f)

    with open(f'/Users/caleb/research/Astro_98/lightkurve_data/heavy_data/perio_data/stars_analyzed.pkl','wb') as f:
        pickle.dump(stars_analyzed,f)


if __name__ == '__main__':
    run_analysis()
