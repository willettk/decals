import urllib.request


def download_image(download_loc, type='jpg', ra=114.5970, dec=21.5681, pixscale=0.262, size=512):
    '''
    Retrieve image from DECALS server and save to disk

    Args:
        download_loc (str): location to save file, excluding type e.g. /data/fits/test_image
        type (str): filetype to retrieve, 'jpg' or 'fits
        ra (float): right ascension (corner? center?)
        dec (float): declination (corner? center?)
        pixscale (float): proportional to decals pixels vs. image pixels. 0.262 for 1-1 map.
        size (int): image edge length in pixels

    Returns:
        None
    '''

    base_url = 'http://legacysurvey.org/viewer/'
    if type == 'jpg':
        base_url += 'jpeg-cutout/?'
    elif type == 'fits':
        base_url += 'fits-cutout/?'
    else:
        raise UserWarning('File type "{}" not recognised'.format(type))

    params = {'ra': ra, 'dec': dec, 'pixscale': pixscale, 'size': size}
    param_string = 'ra={ra}&dec={dec}&layer=decals-dr3&pixscale={pixscale}&bands=grz&size={size}'.format(**params)

    # query and save the response
    print('target loc: ' + r'{}.{}'.format(download_loc, type))
    print('target dir: ' + download_loc)
    urllib.request.urlretrieve(base_url + param_string, r'{}.{}'.format(download_loc, type))

# matrix_rolled = np.moveaxis(matrix, 0, 2)
# print(matrix_rolled.shape)
# plt.imshow(matrix_rolled)
# plt.show()
