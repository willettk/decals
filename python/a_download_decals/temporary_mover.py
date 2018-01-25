import download_decals_settings
import subprocess
from multiprocessing.dummy import Pool as ThreadPool
import os

from tqdm import tqdm
from astropy.table import Table

from get_images.download_images_threaded import get_loc

joint_catalog = Table.read(download_decals_settings.joint_catalog_loc)


def get_old_fits_loc(fits_dir, galaxy):
    return '{}/{}.fits'.format(fits_dir, galaxy['iauname'])

# new_root_fits_dir = '/data/galaxy_zoo/decals/fits/dr5'
new_root_fits_dir = 'Volumes/new_hdd/dr5'


def move_galaxy(galaxy, current_root_fits_dir=download_decals_settings.fits_dir, new_root_fits_dir=new_root_fits_dir):
    current_fits_loc = get_old_fits_loc(current_root_fits_dir, galaxy)
    target_fits_loc = get_loc(new_root_fits_dir, galaxy, 'fits')
    if os.path.exists(current_fits_loc) and not os.path.exists(target_fits_loc):
        command = 'mv {} {}'.format(current_fits_loc, target_fits_loc)
        result = subprocess.Popen(command, shell=True)
        name = result.pid
        try:
            result.wait(timeout=10)
        except:
            print('timeout doing {}, {}'.format(command, name))

n_threads = 100
pool = ThreadPool(n_threads)
list(tqdm(pool.imap(move_galaxy, joint_catalog), total=len(joint_catalog), unit=' moved'))
pool.close()
pool.join()
