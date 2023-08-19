from astropy.io import fits
from glob import glob
from gammapy.irf import Background2D
from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file
import sys
import numpy as np
import os
from tqdm import tqdm

if __name__ == "__main__" :
    '''
        Main function:
        - Search for background files in target dir
        - Search for DL3 files in target dir
        - Append DL3 files with backgrounds when within range
        - Write files to target dir
        - Create index files from target dir
    '''

    # Get the relevent directories
    in_dir = sys.argv[1]
    back_dir = sys.argv[2]
    out_dir = sys.argv[3]

    # Read in the list of backgrounds
    back_files = glob(back_dir + "/background_*_az.fits")
    back_fmt = "background_{ze}_ze_{nsb}_nsb_{az}_az.fits"
    # print (len(back_files))


    # Read in the list of input files
    in_files = glob(in_dir + "/*anasum.fits")
    # print (len(in_files))

    # Check if the directory exists
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


    zen_str = []
    nsb_str = []
    az_str = []
    
    for f in back_files:

        # ['background', '20-30', 'ze', '7.0-8.0', 'nsb', '90.0-180.0', 'az.fits']

        _, ze, _, nsb, _, az, _ = f.split("/")[-1].split("_")
        zen_str.append(ze)
        nsb_str.append(nsb)
        az_str.append(az)
    
    # Get the unique entries
    zen_str = np.unique(zen_str)
    nsb_str = np.unique(nsb_str)
    az_str = np.unique(az_str)

    # Get the mean values
    zen = np.mean([ np.array(ze.split("-"), dtype = float) for ze in zen_str ], axis = 1)
    nsb = np.mean([ np.array(ns.split("-"), dtype = float) for ns in nsb_str ], axis = 1)
    az = np.mean([ np.array(az.split("-"), dtype = float) for az in az_str ], axis = 1)


    # Loop over the runs
    fits_list = []
    for f in tqdm(in_files):
        
        hdul = fits.open(f)
        zen_obs = 90 - float(hdul[1].header["ALT_PNT"])
        nsb_obs = float(hdul[1].header["NSBLEVEL"])
        az_obs = float(hdul[1].header["AZ_PNT"])

        # Get the closest
        # ToDo: Interpolation, Circular means/ Interpolation in cos(angle)
        zen_min = np.argmin(
            np.abs(zen_obs - zen)
        )
        nsb_min = np.argmin(
            np.abs(nsb_obs - nsb)
        )
        az_min = np.argmin(
            np.abs(az_obs - az)
        )

        # Check if background exists
        back_file = back_dir + "/" +  back_fmt.format(
            ze = zen_str[zen_min],
            nsb = nsb_str[nsb_min],
            az = az_str[az_min],
            
        )
        
        if os.path.isfile(back_file):


            # Append the hdul with the background and save to files

            out_file = out_dir + "/" + f.split("/")[-1]
            hdul_back = fits.open(back_file)
            background = Background2D.from_hdulist(hdul_back)
            hdul.append(background.to_table_hdu())
            hdul.writeto(out_file, overwrite = True)

            hdul_back.close()
            fits_list.append(out_file)

        hdul.close()

    # Write the index files
    create_obs_hdu_index_file(
        fits_list,
        out_dir, 
    )



