import numpy as np


def dstn_rgb(imgs, bands, mnmx=None, arcsinh=None, scales=None, desaturate=False):
    '''
    Given a list of image arrays in the given bands, returns a scaled RGB image.

    Args:
        imgs (list): numpy arrays, all the same size, in nanomaggies
        bands (list): strings, eg, ['g','r','z']
        mnmx (min,max), values that will become black/white *after* scaling. Default is (-3,10)):
        arcsinh (bool): if True, use nonlinear scaling (as in SDSS)
        scales (str): Override preset band scaling. Dict of form {band: (plane index, scale divider)}
        desaturate (bool): If [default=False] desaturate pixels dominated by a single colour

    Returns:
        (np.array) of shape (H, W, 3) with values between 0 and 1
    '''

    bands = ''.join(bands)  # stick list of bands into single string

    # TODO refactor grzscales
    # first number is index of that band
    # second number is scale divisor - divide pixel values by scale divisor for rgb pixel value
    grzscales = dict(
        g=(2, 0.0066),
        r=(1, 0.01385),
        z=(0, 0.025),
    )

    if scales is None:
        if bands == 'grz':
            scales = grzscales
        elif bands == 'urz':
            scales = dict(
                u=(2, 0.0066),
                r=(1, 0.01),
                z=(0, 0.025),
            )
        elif bands == 'gri':
            scales = dict(
                g=(2, 0.002),
                r=(1, 0.004),
                i=(0, 0.005),
            )
        else:
            scales = grzscales

    #  create blank matrix to work with
    h, w = imgs[0].shape
    rgb = np.zeros((h, w, 3), np.float32)

    # Copy each band matrix into the rgb image, dividing by band scale divisor to increase pixel values
    for im, band in zip(imgs, bands):
        plane, scale = scales[band]
        rgb[:, :, plane] = (im / scale).astype(np.float32)

    # TODO mnmx -> (min, max)
    # cut-off values for non-linear arcsinh map
    if mnmx is None:
        mn, mx = -3, 10
    else:
        mn, mx = mnmx

    if arcsinh is not None:
        # TODO move out of if statement?
        def nonlinear_map(x):
            # arcsinh (as opposed to np.arcsinh) is parameter (1. from downloader), inverse of softening factor B
            # applies both arcsinh and inverse sqrt rescaling
            return np.arcsinh(x * arcsinh) / np.sqrt(arcsinh)
        # TODO image rescaled by single-pixel not image-pixel, which means colours depend on brightness at this point??
        rgb = nonlinear_map(rgb)
        mn = nonlinear_map(mn)
        mx = nonlinear_map(mx)

    # lastly, rescale image to be between min and max
    rgb = (rgb - mn) / (mx - mn)

    # default False, but downloader sets True
    if desaturate:
        # optionally desaturate pixels that are dominated by a single
        # colour to avoid colourful speckled sky

        # reshape rgb from (h, w, 3) to (3, h, w)
        RGBim = np.array([rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]])
        a = RGBim.mean(axis=0)  # a is mean pixel value across all bands, (h, w) shape
        # putmask: given array and mask, set all mask=True values of array to new value
        np.putmask(a, a == 0.0, 1.0)  # set pixels with 0 mean value to mean of 1. Inplace?
        acube = np.resize(a, (3, h, w))  # copy mean value array (h,w) into 3 bands (3, h, w)
        bcube = (RGBim / acube) / 2.5  # bcube: divide image by mean-across-bands pixel value, and again by 2.5 (why?)
        mask = np.array(bcube)  # isn't bcube already an array?
        wt = np.max(mask, axis=0)  # maximum per pixel across bands of mean-band-normalised rescaled image
        # i.e largest relative deviation from mean
        np.putmask(wt, wt > 1.0, 1.0)  # clip largest allowed relative deviation to one (inplace?)
        wt = 1 - wt  # invert relative deviations
        wt = np.sin(wt*np.pi/2.0)  # non-linear rescaling of relative deviations
        temp = RGBim * wt + a*(1-wt) + a*(1-wt)**2 * RGBim  # multiply by weights in complicated fashion
        rgb = np.zeros((h, w, 3), np.float32)  # reset rgb to be blank
        for idx, im in enumerate((temp[0, :, :], temp[1, :, :], temp[2, :, :])):  # fill rgb with weight-rescaled rgb
            rgb[:, :, idx] = im

    clipped = np.clip(rgb, 0., 1.)  # set max/min to 0 and 1

    return clipped


def get_rgb(imgs, bands, mnmx=None, arcsinh=None, scales=None,
            clip=True):
    '''
    Given a list of images in the given bands, returns a scaled RGB
    image.
    *imgs*  a list of numpy arrays, all the same size, in nanomaggies
    *bands* a list of strings, eg, ['g','r','z']
    *mnmx*  = (min,max), values that will become black/white *after* scaling.
        Default is (-3,10)
    *arcsinh* use nonlinear scaling as in SDSS
    *scales*
    Returns a (H,W,3) numpy array with values between 0 and 1.
    '''
    bands = ''.join(bands)

    grzscales = dict(g=(2, 0.0066),
                     r=(1, 0.01),
                     z=(0, 0.025),
                     )

    # print('get_rgb: bands', bands)

    if scales is None:
        if bands == 'grz':
            scales = grzscales
        elif bands == 'urz':
            scales = dict(u=(2, 0.0066),
                          r=(1, 0.01),
                          z=(0, 0.025),
                          )
        elif bands == 'gri':
            # scales = dict(g = (2, 0.004),
            #               r = (1, 0.0066),
            #               i = (0, 0.01),
            #               )
            scales = dict(g=(2, 0.002),
                          r=(1, 0.004),
                          i=(0, 0.005),
                          )
        elif bands == 'ugriz':
            scales = dict(g=(2, 0.0066),
                          r=(1, 0.01),
                          i=(0, 0.05),
                          )
        else:
            scales = grzscales

    # print('Using scales:', scales)
    h, w = imgs[0].shape
    rgb = np.zeros((h, w, 3), np.float32)
    for im, band in zip(imgs, bands):
        if not band in scales:
            print('Warning: band', band, 'not used in creating RGB image')
            continue
        plane, scale = scales.get(band, (0, 1.))
        # print('RGB: band', band, 'in plane', plane, 'scaled by', scale)
        rgb[:, :, plane] = (im / scale).astype(np.float32)

    if mnmx is None:
        mn, mx = -3, 10
    else:
        mn, mx = mnmx

    if arcsinh is not None:
        def nlmap(x):
            return np.arcsinh(x * arcsinh) / np.sqrt(arcsinh)

        rgb = nlmap(rgb)
        mn = nlmap(mn)
        mx = nlmap(mx)

    rgb = (rgb - mn) / (mx - mn)

    # print(rgb.min(), rgb.max())
    if clip:
        return np.clip(rgb, 0., 1.)
    return rgb


def f_lupton(x, arcsinh=.1):
    return np.arcsinh(x * arcsinh) / np.sqrt(arcsinh)


def lupton_rgb(imgs, bands='grz'):

    grzscales = dict(g=(2, 0.00526),
                     r=(1, 0.008),
                     z=(0, 0.0135),
                     )

    h, w = imgs[0].shape
    img = np.zeros((h, w, 3), np.float32)
    for im, band in zip(imgs, bands):
        plane, scale = grzscales.get(band, (0, 1.))
        img[:, :, plane] = (im / scale).astype(np.float32)

    I = img.sum(axis=2).reshape(424, 424, 1)
    # print(I.shape)

    rescaling = f_lupton(I)/I

    rescaled_img = img * rescaling

    max_rgb = rescaled_img.max()

    # if max_rgb > 1.:
    #     return rescaled_img / max_rgb
    # else:
    #     return rescaled_img

    mn, mx = 0.5, 100.
    mn = f_lupton(mn)
    mx = f_lupton(mx)
    print(mn, mx)
    # mn, mx = rescaled_img.min(), rescaled_img.max()


    # return rescaled_img / rescaled_img.max()

    rescaled_img = (rescaled_img - mn) / (mx - mn)
    # return rescaled_img


    clip = True
    if clip:
        return np.clip(rescaled_img, 0., 1.)

    # TODO unit tests creating rgb image from artificial colour-path fits file,
    # then verifying the colours produced are unique per g-r, r-z colour in input file