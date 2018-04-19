# only used to make DR1 color images - DR2 uses Dustin's function

##############################################################
# Trilogy - Color / grayscale image maker (from FITS input)
# Dan Coe
# http://www.stsci.edu/~dcoe/trilogy/
##############################################################

# Trilogy produces nice color / grayscale images based on your input FITS files
# Just tell the input file (e.g., trilogy.in)
#   which images you'd like applied to which channel (RGB)

# When assigning multiple filters to a channel (R,G,B)
#   Trilogy currently adds the data from those filters

# If you do not specify channels, a grayscale image will be made

# To determine scaling,
#   Samples a (samplesize,samplesize) section of the image core (center)

#################################
# Usage:

# trilogy
# - Assumes an input file trilogy.in

# trilogy acsir.in
# - sets the input file; output name is acsir.png by default

# trilogy acsir.in -indir ../images -outname a383 -showwith PIL -sampledx 300
# - set the input file and override some options from the command line

#################################
# Example input file (everything between the quotes)
# and default values (beginning with outdir):
"""
B
a383_uvis_f225w.fits
a383_uvis_f275w.fits
a383_uvis_f336w.fits
a383_uvis_f390w.fits

G
a383_acs_f435w.fits
a383_acs_f475w.fits
a383_acs_f606w.fits
a383_acs_f625w.fits
a383_acs_f775w.fits
a383_acs_f814w.fits
a383_acs_f850lp.fits

R
a383_ir_f105w.fits
a383_ir_f110w.fits
a383_ir_f125w.fits
a383_ir_f140w.fits
a383_ir_f160w.fits

indir  /Users/dcoe/CLASH/a383/
outname  a383
outdir   .
samplesize 1000
sampledx  0
sampledy  0
stampsize  1000
showstamps  0
satpercent  0.001
noiselum    0.10
colorsatfac  1
deletetests  0
testfirst    1
show  1
showstamps  0
showwith  open
deletetests  0
scaling  None
verbose 1
legend  1
thumbnail  None
maxstampsize  6000
"""

# Can also set noiselum individually for each filter with:
# noiselums 0.10,0.15,0.20

# Some explanations:
"""
    'noiselum':0.10,  # Make higher to dig into noise more (between 0 - 1)
    'satpercent':0.001,  # *Percentage* of pixels which will be saturated
    # (satpercent = 0.001 means 1 / 100,000 pixels will be saturated)
    'samplesize':1000,  # to determine levels
    'sampledx':0,  # offset
    'sampledy':0,  # offset
    'stampsize': 1000,  # for making final color image one section at a time (just a memory issue)
    'maxstampsize':6000,   # My memory can't handle an array larger than 6000x6000
    'showwith':'open',  # Command to display images; set to 0 to display with PIL (as lossy jpeg)
    'scaling':None,  # Use an input scaling levels file
    'legend':1,  # Adds legend to top-left corner indicating which filters were used
"""

#################################

# TO DO:
# Images larger than 6000
# More robust input? if just input 3 fits files, have them be RGB?
# Change temperature to make image redder/bluer, if desired
# I should allow for input of a weight image and make use of it, but I don't currently
# better way to combine images within channel?
# allow one to add a constant to the image before scaling

#################################
# Requirements (python libraries):
# PIL - Python Image Library
# pyfits - FITS handler
# numpy - handles data arrays
# scipy - "golden" root finder (maybe there's another way that doesn't require scipy?)

#################################
# Log scaling constrained at 3 data points: "tri-log-y"

# Inputs:
# % of pixels that saturate
# output brightness of noise
# color saturation boost

# y = log10( k * (x - xo) + 1 ) / r

# input data (x) -> output scaling (y) from 0-1 (in gray / color image)
# x0 yields 0
# x1 yields y1
# x2 yields 1

# Current settings:
# x0: 0 (0 in the input yields black in the output)
# x1: mean + std (1-sigma above the noise)
# x2: set so only some small fraction of pixels saturate (with output = 1)

# DERIVATION

# log10( k * (x - xo) + 1 ) / r
# x0, x1, x2  YIELD  0, 0.5, 1,  RESPECTIVELY

# (0): log10( k (x0 - xo) + 1 ) / r = 0
# (1): log10( k (x1 - xo) + 1 ) / r = 0.5
# (2): log10( k (x2 - xo) + 1 ) / r = 1

# (0)        gives xo = x0
# (2) - (0)  gives r = log( k (x2-x0) + 1 )
# (2) = 2(1) gives k = (x2 - 2*x1 + x0) / (x2 - x0)**2

# This is not easily generalized to output values other than (0, 0.5, 1)
# The first two require y0 = 0
# The last step is possible because y2 / y1 = 1 / 0.5 = 2

# Of course one could always solve numerically...

#################################
# More resources:

# Robert Lupton
# http://www.astro.princeton.edu/~rhl/PrettyPictures/

# Robert Hurt
# http://hea-www.harvard.edu/~ascpub/viztalks/talks/Chandra%20Dynamic%20Range.ppt

#################################

try:
    import pyfits
except ImportError:
    from astropy.io import fits as pyfits
import os
import string
import sys
from glob import glob
from os.path import exists, join

from PIL import Image, ImageDraw
from numpy import *
from scipy.optimize import golden

defaultvalues = {
    'indir':'',
    #'outname':'trilogy',
    'outname':None,
    'outdir':'',
    #'noiselum':0.15,  # Make higher to dig into noise more (between 0 - 1)
    'noiselum':0.10,  # Make higher to dig into noise more (between 0 - 1)
    'noiselums':{},
    'satpercent':0.001,  # *Percentage* of pixels which will be saturated
    # (satpercent = 0.001 means 1 / 100,000 pixels will be saturated)
    'colorsatfac':3,
    'thumbnail':None,
    'samplesize':1000,  # to determine levels
    'sampledx':0,  # offset
    'sampledy':0,  # offset
    'stampsize': 1000,  # for making final color image (just a memory issue)
    'testfirst':1,
    'show':1,
    'showstamps':0,
    'showwith':'open',  # Command to display images; set to 0 to display with PIL (as lossy jpeg)
    'deletetests':0,
    'verbose':1,
    'scaling':None,
    'maxstampsize':6000,   # My memory can't handle an array larger than 6000x6000
    'legend':1,  # Adds legend to top-left corner indicating which filters were used
    }

#################################
# A few general tools

class stat_robust:
    #Generates robust statistics using a sigma clipping
    #algorithm. It is controlled by the parameters n_sigma
    #and n, the number of iterations
    # -from Narciso Benitez
    def __init__(self,x,n_sigma=3,n=5,reject_fraction=None):
        self.x=x
        self.n_sigma=n_sigma
        self.n=n
        self.reject_fraction=reject_fraction

    def run(self):
        good=ones(len(self.x))
        nx=sum(good)
        if self.reject_fraction==None:
            for i in range(self.n):
                if i>0: xs=compress(good,self.x)
                else: xs=self.x
                #            aver=mean(xs)
                aver=median(xs)
                std1=std(xs)
                good=good*less_equal(abs(self.x-aver),self.n_sigma*std1)
                nnx=sum(good)
                if nnx==nx: break
                else: nx=nnx
        else:
            np=float(len(self.x))
            nmin=int((0.5*self.reject_fraction)*np)
            nmax=int((1.-0.5*self.reject_fraction)*np)
            orden=argsort(self.x)
            connect(arange(len(self.x)),sort(self.x))
            good=greater(orden,nmin)*less(orden,nmax)

        self.remaining=compress(good,self.x)
        self.max=max(self.remaining)
        self.min=min(self.remaining)
        self.mean=mean(self.remaining)
        self.rms=std(self.remaining)
        #self.rms0=rms(self.remaining)  # --DC
        self.median=median(self.remaining)
        self.outliers=compress(logical_not(good),self.x)
        self.n_remaining=len(self.remaining)
        self.n_outliers=len(self.outliers)
        self.fraction=1.-(float(self.n_remaining)/float(len(self.x)))


def rms(x):
    return sqrt(mean(x**2))

class meanstd_robust:
    #Generates robust statistics using a sigma clipping
    #algorithm. It is controlled by the parameters n_sigma
    #and n, the number of iterations
    # ADAPTED from Txitxo's stat_robust
    # Now much quicker for large arrays
    def __init__(self,x,n_sigma=3,n=5,sortedalready=False):
        self.x=x
        self.n_sigma=n_sigma
        self.n=n
        self.sortedalready = sortedalready

    def run(self):
        ihi = nx = len(self.x)
        ilo = 0
        #self.x[isnan(self.x)]=0  # set all nan values to zero
        if not self.sortedalready:
            print 'sorting...'
            self.xsort = sort(self.x)
        else:
            self.xsort = self.x
        #print self.xsort
        for i in range(self.n):
            #print i
            xs = self.xsort[ilo:ihi]
            #print 'median'
            #aver = median(xs)
            #print xs
            #print xs[-1]
            #print xs[-2]
            #print len(xs)
            imed = (ilo+ihi) / 2
            #print imed
            aver = xs[imed]
            #print 'std'
            std1 = std(xs)
            std1 = rms(xs - aver)
            #print 'lohi'
            lo = aver - self.n_sigma * std1
            hi = aver + self.n_sigma * std1
            #print 'searching...'
            ilo = searchsorted(self.xsort, lo)
            ihi = searchsorted(self.xsort, hi, side='right')
            nnx = ihi - ilo
            #print ilo, ihi, nnx, nx, lo, hi
            if nnx==nx: break
            else: nx=nnx

        self.remaining = xrem = xs[ilo:ihi]
        self.mean = mean(xrem)
        self.std  = rms(xrem - self.mean)


def strend(str, phr):
    return str[-len(phr):] == phr

def decapfile(name, ext=''):
    """REMOVE EXTENSION FROM FILENAME IF PRESENT
    IF ext LEFT BLANK, THEN ANY EXTENSION WILL BE REMOVED"""
    if ext:
        if ext[0] <> '.':
            ext = '.' + ext
        n = len(ext)
        if name[-n:] == ext:
            name = name[:-n]
    else:
        if strend(name, '.gz'):
            name = name[-3:]
        i = string.rfind(name, '.')
        if i > -1:
            name = name[:i]
    return name

def stringsplitatof(str, separator=''):
    """Splits a string into floats"""
    if separator:
        words = string.split(str, separator)
    else:
        words = string.split(str)
    vals = []
    for word in words:
        vals.append(string.atof(word))
    return vals

def str2num(str, rf=0):
    """CONVERTS A STRING TO A NUMBER (INT OR FLOAT) IF POSSIBLE
    ALSO RETURNS FORMAT IF rf=1"""
    try:
        num = string.atoi(str)
        format = 'd'
    except:
        try:
            num = string.atof(str)
            format = 'f'
        except:
            if not string.strip(str):
                num = None
                format = ''
            else:
                words = string.split(str)
                if len(words) > 1:
                    num = map(str2num, tuple(words))
                    format = 'l'
                else:
                    num = str
                    format = 's'
    if rf:
        return (num, format)
    else:
        return num

def clip2(m, m_min=None, m_max=None):
    if m_min == None:
        m_min = min(m)
    if m_max == None:
        m_max = max(m)
    return clip(m, m_min, m_max)

def striskey(str):
    """IS str AN OPTION LIKE -C or -ker
    (IT'S NOT IF IT'S -2 or -.9)"""
    iskey = 0
    if str:
        if str[0] == '-':
            iskey = 1
            if len(str) > 1:
                iskey = str[1] not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.']
    return iskey

def params_cl(converttonumbers=True):
    """RETURNS PARAMETERS FROM COMMAND LINE ('cl') AS DICTIONARY:
    KEYS ARE OPTIONS BEGINNING WITH '-'
    VALUES ARE WHATEVER FOLLOWS KEYS: EITHER NOTHING (''), A VALUE, OR A LIST OF VALUES
    ALL VALUES ARE CONVERTED TO INT / FLOAT WHEN APPROPRIATE"""
    list = sys.argv[:]
    i = 0
    dict = {}
    oldkey = ""
    key = ""
    list.append('')  # EXTRA ELEMENT SO WE COME BACK AND ASSIGN THE LAST VALUE
    while i < len(list):
        if striskey(list[i]) or not list[i]:  # (or LAST VALUE)
            if key:  # ASSIGN VALUES TO OLD KEY
                if value:
                    if len(value) == 1:  # LIST OF 1 ELEMENT
                        value = value[0]  # JUST ELEMENT
                dict[key] = value
            if list[i]:
                key = list[i][1:] # REMOVE LEADING '-'
                value = None
                dict[key] = value  # IN CASE THERE IS NO VALUE!
        else: # VALUE (OR HAVEN'T GOTTEN TO KEYS)
            if key: # (HAVE GOTTEN TO KEYS)
                if value:
                    if converttonumbers:
                        value.append(str2num(list[i]))
                    else:
                        value = value + ' ' + list[i]
                else:
                    if converttonumbers:
                        value = [str2num(list[i])]
                    else:
                        value = list[i]
        i += 1

    return dict

#################################
# TRILOGY-specific tools

def determinescaling(data, unsatpercent):
    """Determines data values (x0,x1,x2) which will be scaled to (0,noiselum,1)"""
    datasorted = sort(data.flat)
    datasorted[isnan(datasorted)]=0  # set all nan values to zero
    if datasorted[0] == datasorted[-1]:
        levels = 0, 1, 100  # whatever
    else:
        s = meanstd_robust(datasorted,sortedalready=True)
        s.run()
        m = s.mean
        r = s.std

        x0 = 0
        x1 = m+r
        x2 = setlevels(datasorted, array([unsatpercent]), sortedalready=True)[0]
        levels = x0, x1, x2
    return levels
    
# PREVIOUSLY in colorimage.py
def setlevels(data, pp, stripneg=False, sortedalready=False):
    if sortedalready:
        vs = data
    else:
        vs = sort(data.flat)
    if stripneg:  # Get rid of negative values altogether!
        # This is the way I was doing it for a while
        # Now that I'm not, resulting images should change (get lighter)
        i = searchsorted(vs, 0)
        vs = vs[i+1:]
    else:  # Clip negative values to zero
        vs = clip2(vs, 0, None)
    ii = array(pp) * len(vs)
    ii = ii.astype(int)
    ii = clip(ii, 0, len(vs)-1)
    levels = vs.take(ii)
    return levels

def imscale1(data, levels):
    # x0, x1, x2  YIELD  0, 0.5, 1,  RESPECTIVELY
    x0, x1, x2 = levels  
    k = (x2 - 2 * x1 + x0) / float(x1 - x0) ** 2
    r1 = log10( k * (x2 - x0) + 1)
    v = ravel(data)
    v = clip2(v, 0, None)
    d = k * (v - x0) + 1
    d = clip2(d, 1e-30, None)
    z = log10(d) / r1
    z = clip(z, 0, 1)
    z.shape = data.shape
    z = z * 255
    z = z.astype(int)
    return z

def da(k):
    a1 = k * (x1 - x0) + 1
    a2 = k * (x2 - x0) + 1
    a1n = a1**n
    a1n = abs(a1n)  # Don't want the solutions where a1 & a2 are both negative!
    da1 = a1n - a2
    k = abs(k)
    if k == 0:
        return da(1e-10)
    else:
        da1 = da1 / k  # To avoid solution k = 0!
    return abs(da1)

# Fixed: (eliminated the negative solution in da above)
# For some reason, setting noiselum = 0.2 (exactly) was making an all yellow image
# it alters k for some of the channels
# levels stay the same
def imscale2(data, levels, y1):
    # x0, x1, x2  YIELD  0, y1, 1,  RESPECTIVELY
    global n, x0, x1, x2  # So that golden can use them
    x0, x1, x2 = levels  
    if y1 == 0.5:
        k = (x2 - 2 * x1 + x0) / float(x1 - x0) ** 2
    else:
        n = 1 / y1
        k = abs(golden(da))
    r1 = log10( k * (x2 - x0) + 1)
    v = ravel(data)
    v = clip2(v, 0, None)
    d = k * (v - x0) + 1
    d = clip2(d, 1e-30, None)
    z = log10(d) / r1
    z = clip(z, 0, 1)
    z.shape = data.shape
    z = z * 255
    z = z.astype(uint8)
    return z

#########

def satK2m(K):
    m00 = rw * (1-K) + K
    m01 = gw * (1-K)
    m02 = bw * (1-K)
    
    m10 = rw * (1-K)
    m11 = gw * (1-K) + K
    m12 = bw * (1-K)
    
    m20 = rw * (1-K)
    m21 = gw * (1-K)
    m22 = bw * (1-K) + K
    
    m = array([[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]])
    return m

# Luminance vector
# All pretty similar; yellow galaxy glow extended a bit more in NTSC
rw, gw, bw = 0.299,  0.587,  0.114  # NTSC (also used by PIL in "convert")
rw, gw, bw = 0.3086, 0.6094, 0.0820  # linear
rw, gw, bw = 0.212671, 0.715160, 0.072169  # D65: red boosted, blue muted a bit, I like it

# also see PIL's ImageEnhance.Contrast
def adjsat(RGB, K):
    """Adjust the color saturation of an image.  K > 1 boosts it."""
    m = satK2m(K)
    three, nx, ny = RGB.shape
    RGB.shape = three, nx*ny
    RGB = dot(m, RGB)
    RGB.shape = three, nx, ny
    return RGB

# Now using the coeim.py version rather than the colorimage.py version
def RGB2im(RGB):
    """r, g, b = data  (3, ny, nx)
    Converts to an Image"""
    data = RGB
    data = transpose(data, (1,2,0))  # (3, ny, nx) -> (ny, nx, 3)
    data = clip(data, 0, 255)
    data = data.astype(uint8)
    three = data.shape[-1]  # 3 if RGB, 1 if L
    if three == 3:
        im = Image.fromarray(data)
    elif three == 1:
        im = Image.fromarray(data[:,:,0], 'L')
    else:
        print 'Data shape not understood: expect last number to be 3 for RGB, 1 for L', data.shape
        raise Exception  # Raise generic exception and exit
    
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    return im

def RGBscale2im(RGB, levdict, noiselums, colorsatfac, mode='RGB'):
    three, nx, ny = RGB.shape  # if 'L', then three = 1 !
    if nx * ny > 2000 * 2000:
        print 'Warning: You should probably feed smaller stamps into RGBscale2im.'
        print "This may take a while..."

    scaled = zeros(RGB.shape, float)
    for i in range(three):
        channel = mode[i]  # 'RGB' or 'L'
        levels = levdict[channel]
        noiselum = noiselums[channel]
        scaled[i] = imscale2(RGB[i], levels, noiselum)

    if (colorsatfac <> 1) and (mode == 'RGB'):
        scaled = adjsat(scaled, colorsatfac)
    
    im = RGB2im(scaled)
    return im


def grayimage(scaled):
    ny, nx = scaled.shape
    im = Image.new('L', (nx,ny))
    im.putdata(scaled.ravel())
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    return im

def grayscaledimage(stamp, levels, noiselum):
    global scaled
    scaled = imscale2(stamp, levels, noiselum)
    im = grayimage(scaled)
    return im


def savelevels(levdict, outfile='levels.txt', outdir=''):
    outfile = join(outdir, outfile)
    fout = open(outfile, 'w')
    for filt in levdict.keys():
        levels = x0, x1, x2 = levdict[filt]

        line = '%s  %g  %g  %g' % (filt, x0, x1, x2)
        fout.write(line+'\n')

    fout.close()



def loadfile(filename, dir="", silent=1, keepnewlines=0):
    infile = join(dir, filename)
    if not silent:
        print "Loading ", infile, "...\n"
    fin = open(infile, 'r')
    sin = fin.readlines()
    fin.close()
    if not keepnewlines:
        for i in range(len(sin)):
            sin[i] = sin[i][:-1]
    return sin

def loaddict(filename, dir="", silent=0):
    lines = loadfile(filename, dir, silent)
    dict = {}
    for line in lines:
        if line[0] <> '#':
            words = string.split(line)
            key = str2num(words[0])
            val = ''  # if nothing there
            valstr = string.join(words[1:], ' ')
            valtuple = False
            valarray = True
            if valstr[0] in '[(' and valstr[-1] in '])':  # LIST / TUPLE!
                valtuple = valstr[0] == '('
                valstr = valstr[1:-1].replace(',', '')
                words[1:] = string.split(valstr)
            if len(words) == 2:
                val = str2num(words[1])
            elif len(words) > 2:
                val = []
                for word in words[1:]:
                    val.append(str2num(word))
                if valtuple:
                    val = tuple(val)
                if valarray:
                    val = array(val)
                
            dict[key] = val
    return dict


#################################
# Apply offsets
# Implement below??

offsets = {}

def offsetarray(data, offset):
    new = zeros(data.shape)
    dx, dy = offset
    if dy >= 0:
        if dx >= 0:
            new[dy:,dx:] = data[:-dy,:-dx]
        else:
            new[dy:,:-dx] = data[:-dy,dx:]
    else:
        if dx >= 0:
            new[:-dy,dx:] = data[dy:,:-dx]
        else:
            new[:-dy,:-dx] = data[dy:,dx:]
    return new

for channel in offsets.keys():
    dataRGB[channel] = offsetarray(dataRGB[channel], offsets[channel])

#################################


#################################

def extractfilter(header):
    """Extracts filter from a FITS header"""
    filt = header.get('FILTER', '')
    if filt == '':
        filt1 = header.get('FILTER1', '')
        if filt1 <> '':
            if filt1[:5] == 'CLEAR':
                filt2 = header.get('FILTER2', '')
                filt = filt2
            else:
                filt = filt1
    return filt

def loadfitsimagedata(image, indir='', silent=1):
    global imfilt, imfilts
    if image[-1] == ']':
        iext = int(image[-2])
        image = image[:-3]  # Remove [0]
    else:
        iext = 0

    image0 = decapfile(image)

    image = join(indir, image)
    hdulist = pyfits.open(image, memmap=1)
    hdu = hdulist[iext]
    data = hdu.data
    if image0 in imfilts.keys():
        imfilt = imfilts[image0]
    else:
        imfilt = extractfilter(hdu.header)
        imfilts[image0] = imfilt
    if not silent:
        print image+'[%d]' % iext, data.shape
    
    return data

imfilts = {}

def processimagename(image):
    global imfilts
    if image[-1] == ')':
        i = image.find('(')
        filt = image[i+1:-1]
        image = image[:i]
        imfilts[image] = filt
    if image[-1] == ']':
        ext = image[-3:]
        image = image[:-3]
    else:
        ext = ''
    if not strend(image, 'fits') and not strend(image, 'fits.gz'):
        image += '.fits'
    image = image + ext
    return image

class Trilogy:
    def __init__(self, infile=None, images=None, imagesorder='BGR', **inparams):
        self.nx = None  # image size
        self.ny = None  # image size
        self.xlo = 0
        self.ylo = 0
        self.xhi = None
        self.yhi = None
        self.imagesRGB = {'R':[], 'G':[], 'B':[], 'L':[]}  # File names
        self.inkeys = []
        self.mode = 'L'  # reset below if color
        self.weightext = None  # No weighting unless weight images are declared
        # Can use either:
        # weightext drz wht
        # weightext drz -> wht

        self.infile = infile
        if infile:
            self.loadinputs()

        self.images = images
        if images:
            self.setimages()

        self.inparams = inparams
        if inparams:
            self.setinparams()

        #if infile or inparams:
        self.setdefaults()

        #self.setoutfile()
        
    def setinparams(self):
        if self.verbose:
            print 'From input parameters:'
        for key in self.inparams.keys():
            val = self.inparams[key]
            cmd = 'self.%s = val' % key
            exec(cmd)
            if self.verbose:
                print key, '=', val
            self.inkeys.append(key)

    def setdefaults(self):
        if self.verbose:
            print 'Default:'
        for key in defaultvalues.keys():
            if key not in self.inkeys:
                val = defaultvalues[key]
                cmd = 'self.%s = val' % key
                exec(cmd)
                if self.verbose:
                    print key, '=', val, '(default)'

    def setimages(self, images=None):
        images = images or self.images
        if images <> None:
            if type(images) == str:  # Single image
                images = processimagename(images)
                self.imagesRGB['L'] = [images]
                self.mode = 'L'
            elif type(images[0]) == str:  # List of images
                images = map(processimagename, images)
                self.imagesRGB['L'] = images
                self.mode = 'L'
            else:  # List of 3 lists of images, one for each channel
                self.mode = 'RGB'
                for i in range(3):
                    channel = imagesorder[i]
                    channelimages = map(processimagename, images[i])
                    self.imagesRGB[channel] = channelimages

    def setnoiselums(self):
        for channel in self.mode:
            self.noiselums[channel] = self.noiselum
            
    def loadinputs(self):
        """Load R,G,B filenames and options"""

        f = open(self.infile)
        prevline = ''
        channel = 'L'  # if no channel declared, then it's grayscale!
        self.noiselums = {}
        for line in f:
            if line[0] == '#':
                continue
            
            word = string.strip(line)
            if len(word):
                words = string.split(word)
                if len(words) == 1:  # Channel or image name
                    if (word in 'RGB') and (prevline == ''):
                        channel = word
                        self.mode = 'RGB'
                    else:
                        image = word
                        image = processimagename(image)
                        self.imagesRGB[channel].append(image)
                else:  # parameter and value(s)
                    key = words[0]
                    val = str2num(string.join(words[1:]))
                    if key == 'weightimages':
                        if len(words[1:]) == 2:  # drz wht
                            keyagain, self.imext, self.weightext = words
                        else:  # drz -> wht
                            self.imext = words[1]
                            self.weightext = words[3]
                    elif key == 'noiselums':
                        if ',' in val:
                            val = val.split(',')
                            val = map(float, val)
                        for i, channel in enumerate(self.mode[::-1]):
                            self.noiselums[channel] = val[i]
                    else:
                        cmd = 'self.%s = val' % key
                        exec(cmd)
                    self.inkeys.append(key)
            prevline = word
        
        if self.noiselums == {}:
            if 'noiselum' in self.inkeys:
                for channel in self.mode:
                    self.noiselums[channel] = self.noiselum
                    self.inkeys.append('noiselums')

        f.close()

    def setoutfile(self, outname=None):
        self.outname = outname or self.outname or self.infile
        if self.outname == None:
            self.outname = self.images
            if self.outname:
                self.outname = os.path.basename(self.outname)
                self.outname = decapfile(self.outname)
        
        if (len(self.outname) > 4) and (self.outname[-4] == '.'):
            # Has extension
            self.outfile = self.outname  # Use whatever extension they picked
            self.outname = self.outname[:-4]  # Remove extension
        else:  # Just root
            self.outfile = self.outname + '.png'

    def outfilterfile(self):
        outfile = 'trilogy_filters.txt'
        outfile = join(self.outdir, outfile)
        return outfile

    def loadimagesize(self):
        global filters
        if self.verbose:
            print "Loading image data.",
            print "If multiple filters per channel, adding data."
        filters = {'B':[], 'G':[], 'R':[], 'L':[]}
        for ichannel, channel in enumerate(self.mode):
            for image in self.imagesRGB[channel]:
                if self.verbose:
                    print channel,
                data = loadfitsimagedata(image, self.indir, silent=logical_not(self.verbose))
                filters[channel].append(imfilt)
                ny, nx = data.shape
                if self.ny == None:
                    self.ny = ny
                    self.nx = nx
                    self.yc = ny / 2
                    self.xc = nx / 2
                else:
                    if (self.ny <> ny) or (self.nx <> nx):
                        print "Input FAIL.  Your images are not all the same size as (%d,%d)." % (self.ny, self.nx)
                        for channel in self.mode[::-1]:  # 'BGR'
                            for image in imagesRGB[channel]:
                                data = loadfitsimagedata(image, self.indir, silent=logical_not(self.verbose))
                        raise  # Raise Exception (error) and quit


        fout = open(self.outfilterfile(), 'w')
        for channel in self.mode[::-1]:
            filtstr = string.join(filters[channel], ' + ')
            if self.verbose:
                print channel, '=',
                print filtstr
            fout.write('%s = %s\n' % (channel, filtstr))

        fout.close()

        if 0:
            fout = open(self.outfilterfile(), 'w')
            for channel in 'BGR':
                filtstr = string.join(filters[channel], ' + ')
                fout.write('%s = %s\n' % (channel, filtstr))

            fout.close()

        # if hi is negative, then trims from that edge (as 1000:-1000)
        if self.xhi == None:
            self.xhi = self.nx
        elif self.xhi < 0:
            self.xhi = self.nx + self.xhi

        if self.yhi == None:
            self.yhi = self.ny
        elif self.yhi < 0:
            self.yhi = self.ny + self.yhi


    def loadstamps(self, limits,silent=1):
        ylo, yhi, xlo, xhi = limits

        ylo = clip(ylo, 0, self.ny)
        yhi = clip(yhi, 0, self.ny)
        xlo = clip(xlo, 0, self.nx)
        xhi = clip(xhi, 0, self.nx)

        ny = yhi - ylo
        nx = xhi - xlo
        
        three = len(self.mode)
        stampRGB = zeros((three, ny, nx), float)
        if self.weightext:
            weightstampRGB = zeros((three, ny, nx), float)
        for ichannel, channel in enumerate(self.mode):
            for image in self.imagesRGB[channel]:
                if self.verbose:
                    print channel,
                data = loadfitsimagedata(image, self.indir, silent=silent)
                stamp = data[ylo:yhi,xlo:xhi]

                # weight image?
                if self.weightext <> None:
                    weightimage = image.replace(self.imext, self.weightext)
                    weightfile = join(self.indir, weightimage)
                    if exists(weightfile):
                        weight = loadfitsimagedata(weightimage, self.indir, silent=silent)
                        weightstamp = weight[ylo:yhi,xlo:xhi]
                        weightstamp = greater(weightstamp, 0)  # FLAG IMAGE!!  EITHER 1 or 0
                        weightstampRGB[ichannel] = weightstampRGB[ichannel] + weightstamp
                        stamp = stamp * weightstamp
                    else:
                        print weightfile, 'DOES NOT EXIST'

                stampRGB[ichannel] = stampRGB[ichannel] + stamp

        if self.weightext <> None:
            for ichannel, channel in enumerate(self.mode):
                stampRGB[ichannel] = where(weightstampRGB[ichannel],
                                           stampRGB[ichannel] / weightstampRGB[ichannel], 0)

        return stampRGB

    def determinescalings(self):
        """Determine data scalings
        will sample a (samplesize x samplesize) region of the (centered) core
        make color image of the core as a test if desired"""

        self.testimages = []
        redo = True
        while redo:  # Until user is happy with test image of core
            dx = dy = self.samplesize
            unsatpercent = 1 - 0.01 * self.satpercent
            self.levdict = {}

            if dx * dy == 0:
                print 'By setting samplesize = 0, you have asked to sample the entire image to determine the scalings.'
                print '(Note this will be clipped to a maximum of %dx%d.)' % (self.maxstampsize, self.maxstampsize)
                dx = dy = self.maxstampsize  # Maximum size possible
            
            ylo = clip(self.yc-dy/2 + self.sampledy, 0, self.ny)
            yhi = clip(self.yc+dy/2 + self.sampledy, 0, self.ny)
            xlo = clip(self.xc-dx/2 + self.sampledx, 0, self.nx)
            xhi = clip(self.xc+dx/2 + self.sampledx, 0, self.nx)
            dy = yhi - ylo
            dx = xhi - xlo 
            if self.verbose:
                print "Determining image scaling based on %dx%d core sample" % (dx, dy),
            if self.sampledx or self.sampledy:
                if self.verbose:
                    print 'offset by (%d,%d)' % (self.sampledx, self.sampledy),
            
            limits = ylo, yhi, xlo, xhi
            stampRGB = self.loadstamps(limits)
            for ichannel, channel in enumerate(self.mode):
                self.levdict[channel] = determinescaling(stampRGB[ichannel], unsatpercent)
            
            savelevels(self.levdict, outdir=self.outdir)
            
            redo = False
            if self.testfirst:

                im = RGBscale2im(stampRGB, self.levdict, self.noiselums, self.colorsatfac, self.mode)
                
                outfile = '%s_test_%g_%g_%g.png' % (self.outname, self.satpercent, self.noiselum, self.colorsatfac)
                outfile = join(self.outdir, outfile)
                self.testimages.append(outfile)

                print "Creating test image", outfile
                im.save(outfile)

                # NOTE I use "open" instead of im.show()
                # because the latter converts the image to a jpg for display
                # which degrades it slightly
                if self.show:
                    self.showimage(outfile, Image)

                print 'Like what you see?'
                print 'If so, press <Enter> a few times'

                print 'Otherwise, enter new values:'

                line = '  noise yields brightness: %g? ' % self.noiselum
                inp = raw_input(line)
                if string.strip(inp) <> '':
                    self.noiselum = float(inp)
                    for channel in self.mode:
                        self.noiselums[channel] = self.noiselum
                    redo = True

                line = '  %% of pixels that saturate: %g? ' % self.satpercent
                inp = raw_input(line)
                if string.strip(inp) <> '':
                    self.satpercent = float(inp)
                    redo = True

                if self.mode == 'RGB':
                    line = '  color saturation factor: %g? ' % self.colorsatfac
                    inp = raw_input(line)
                    if string.strip(inp) <> '':
                        self.colorsatfac = float(inp)
                        redo = True

                line = '  Sample size: %d? ' % self.samplesize
                inp = raw_input(line)
                if string.strip(inp) <> '':
                    self.samplesize = int(inp)
                    redo = True

                line = '  Sample offset x: %d? ' % self.sampledx
                inp = raw_input(line)
                if string.strip(inp) <> '':
                    self.sampledx = int(inp)
                    redo = True

                line = '  Sample offset y: %d? ' % self.sampledy
                inp = raw_input(line)
                if string.strip(inp) <> '':
                    self.sampledy = int(inp)
                    redo = True

    def makecolorimage(self):
        """Make color image (in sections)"""
        if (self.stampsize == self.samplesize == 0) and self.testfirst:
            # Already did the full image!
            print 'Full size image already made.'
            imfile = self.testimages[-1]
            outfile = join(self.outdir, self.outfile)
            if self.deletetests:
                print 'Renaming to', outfile
                os.rename(imfile, outfile)
            else:
                print 'Copying to', outfile
                os.copy(imfile, outfile)
            imfull = Image.open(outfile)
            return imfull
        
        # Clean up: Delete test images
        if self.deletetests:
            for testimage in self.testimages:
                if exists(testimage):
                    os.remove(testimage)

        dx = dy = self.stampsize
        if dx * dy == 0:
            dx = dy = self.maxstampsize
        
        imfull = Image.new(self.mode, (self.nx, self.ny))

        if self.mode == 'RGB':
            if self.verbose:
                print 'Making full color image, one stamp (section) at a time...'
        elif self.mode == 'L':
            if self.verbose:
                print 'Making full grayscale image, one stamp (section) at a time...'
        for yo in range(self.ylo,self.yhi,dy):
            dy1 = min([dy, self.yhi-yo])
            for xo in range(self.xlo,self.xhi,dx):
                dx1 = min([dx, self.xhi-xo])
                if self.verbose:
                    print '%5d, %5d  /  (%d x %d)' % (xo, yo, self.nx, self.ny)
                limits = yo, yo+dy, xo, xo+dx
                stamps = self.loadstamps(limits)
                im = RGBscale2im(stamps, self.levdict, self.noiselums, self.colorsatfac, self.mode)
                if self.show and self.showstamps:
                    im.show()

                imfull.paste(im, (xo,self.ny-yo-dy1,xo+dx1,self.ny-yo))

        outfile = join(self.outdir, self.outfile)
        if self.legend:
            self.addlegend(im=imfull)
        else:
            if self.verbose:
                print 'Saving', outfile, '...'
            imfull.save(outfile)
        
        if self.show:
            self.showimage(outfile, Image)
        
        return imfull

    def makethumbnail1(self, outroot, width, fmt='jpg'):
        nx = width
        ny = 1000 * (ny / float(nx))
        im = open(outroot+'.png')
        im2 = im.resize((nx,ny))
        im2.save(self.outname+'_%d.%s' % (width, fmt))
        return im2

    def makethumbnail(self):
        if self.thumbnail not in [None, 'None']:
            outname = self.thumbnail
            if outname[-4] == '.':
                outname = outname[:-4]
                fmt = outname[-3:]
            width = int(outname)
            self.makethumbnail1(self.outname, width, fmt)

    def showsample(self, outfile):
        dx = dy = self.samplesize
        if dx * dy == 0:
            if self.verbose:
                print 'By setting samplesize = 0, you have asked to sample the entire image to determine the scalings.'
                print '(Note this will be clipped to a maximum of %dx%d.)' % (self.maxstampsize, self.maxstampsize)
            dx = dy = self.maxstampsize  # Maximum size possible
        
        ylo = clip(self.yc-dy/2 + self.sampledy, 0, self.ny)
        yhi = clip(self.yc+dy/2 + self.sampledy, 0, self.ny)
        xlo = clip(self.xc-dx/2 + self.sampledx, 0, self.nx)
        xhi = clip(self.xc+dx/2 + self.sampledx, 0, self.nx)
        dy = yhi - ylo
        dx = xhi - xlo 
        if self.verbose:
            print "Showing %dx%d core sample" % (dx, dy),
        if (self.sampledx or self.sampledy) and self.verbose:
            print 'offset by (%d,%d)' % (self.sampledx, self.sampledy),

        limits = ylo, yhi, xlo, xhi
        stampRGB = self.loadstamps(limits)
        im = RGBscale2im(stampRGB, self.levdict, self.noiselums, self.colorsatfac, self.mode)
        
        outfile = join(self.outdir, outfile)
        
        if self.verbose:
            print "Creating test image", outfile
        im.save(outfile)
        
        # NOTE I use "open" instead of im.show()
        # because the latter converts the image to a jpg for display
        # which degrades it slightly
        if self.show:
            self.showimage(outfile, Image)

        if self.verbose:
            print 'Like what you see?'
        inp = raw_input()

    def addlegend(self, outfile=None, im=None):

        if im == None:
            outfile1 = join(self.outdir, self.outfile)
            if self.verbose:
                print 'Adding legend to', outfile1, '...'
            im = Image.open(outfile1)
        else:
            if self.verbose:
                print 'Adding legend...'
        
        nx, ny = im.size
        draw = ImageDraw.Draw(im)

        x = 20
        y0 = 20
        dy = 15

        txt = loadfile(self.outfilterfile(), silent=1)
        
        if self.mode == 'L':
            white = 255
            line = txt[0][4:]  # get rid of leading "L = "
            draw.text((x, y0), line, fill=white)
        else:
            blue   = tuple(255 * array([0,0.5,1]))
            green  = tuple(255 * array([0,1,0]))
            red    = tuple(255 * array([1,0,0]))

            colors = blue, green, red
            colors = array(colors).astype(int)

            #print 'ImageDraw will complain about float color but handles it fine:'
            for i, line in enumerate(txt):
                y = y0 + dy*i
                color = tuple(colors[i])
                draw.text((x, y), ''.join(line.split("=")[1:]).strip(), fill=color)

        if outfile == None:
            outfile = join(self.outdir, self.outfile)
        
        if self.verbose:
            print 'Saving', outfile, '...'
        im.save(outfile)


    def showimage(self, outfile, Image):
        cmd = self.showwith
        if (not cmd) or (cmd.upper() == 'PIL'):
            Image.open(outfile).show()
        else:
            try:
                os.system(cmd + ' ' + outfile)
            except:  # In case "open" doesn't work on their system (not a Mac)
                # Although it may not work but not raise an error either!
                # Should do better error handling here
                Image.open(outfile).show()

    def addtofilterlog(self):
        fout = open('trilogyfilterlog.txt', 'a')
        fout.write(self.outfile+'\n')
        txt = loadfile(self.outfilterfile())
        for line in txt:
            fout.write(line+'\n')
        
        fout.write('\n')
        fout.close()

    def run(self):
        self.setimages()  # not needed from command line
        self.setoutfile()  # adds .png if necessary to outname
        self.loadimagesize()
        self.addtofilterlog()
        if 'justaddlegend' in self.inkeys:
            self.addlegend()
            quit()

        if self.noiselums == {}:
            self.setnoiselums()
        if self.scaling == None:
            self.determinescalings()
        else:
            if self.verbose:
                print 'Loading scaling saved in', self.scaling
            self.levdict = loaddict(self.scaling)
            self.showsample(self.outname+'_'+self.scaling[:-4]+'.png')
        if self.verbose:
            print "Scalings:"
        for channel in self.mode:
            if self.verbose:
                print channel, self.levdict[channel]
        self.makecolorimage()
        
        self.makethumbnail()

def pause(text=''):
    inp = raw_input(text)

if __name__ == '__main__':    
    if len(sys.argv) == 1:
        infile = 'trilogy.in'
        images = None
        Trilogy(infile, images=images, **params_cl()).run()
    else: # > 1
        input1 = sys.argv[1]
        if ('*' in input1) or ('?' in input1):
            indir = params_cl().get('indir', '')
            input1 = join(indir, input1)
            images = glob(input1)
            for image in images:
                Trilogy(images=image, **params_cl()).run()
        else:
            images = None
            if strend(input1, '.fits') or strend(input1, '.fits.gz'):
                images = input1
                infile = None
            else:
                infile = input1
            
            Trilogy(infile, images=images, **params_cl()).run()
