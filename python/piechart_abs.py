from matplotlib import pyplot as plt
from gz_class import plurality
import numpy as np
from collections import Counter
import re
import pandas as pd
from astropy.io import fits

gzpath = '/Users/willettk/Astronomy/Research/GalaxyZoo'

"""

piechart
=========

Generate pie chart of the plurality classifications for the DECaLS data.

Kyle Willett (UMN) - 10 Dec 2015

"""

def survey_dict():

    # Information about the specific group settings in the project

    d = {u'candels':        {'name':u'CANDELS','retire_limit':80},
        u'candels_2epoch':  {'name':u'CANDELS 2-epoch','retire_limit':80},
        u'decals':          {'name':u'DECaLS','retire_limit':40},
        u'ferengi':         {'name':u'FERENGI','retire_limit':40},
        u'goods_full':      {'name':u'GOODS full-depth','retire_limit':40},
        u'illustris':       {'name':u'Illustris','retire_limit':40},
        u'sloan_singleband':{'name':u'SDSS single-band','retire_limit':40},
        u'ukidss':          {'name':u'UKIDSS','retire_limit':40},
        #u'sloan':          {'name':u'SDSS DR8','retire_limit':60},
        u'stripe82':        {'name':u'Stripe 82','retire_limit':40},
        u'gz2':             {'name':u'SDSS DR7','retire_limit':40}}

    return d
    
def is_number(s):

    # Is a string a representation of a number?

    try:
        int(s)
        return True
    except ValueError:
        return False

def morphology_distribution(survey='decals'):

    # What's the plurality distribution of morphologies?

    # Get weights
    try:
        collation_file = "{0}/gz_reduction_sandbox/data/decals_unweighted_classifications_00.csv".format(gzpath)
        collated = pd.read_csv(collation_file)
    except IOError:
        print "Collation file for {0:} does not exist. Aborting.".format(survey)
        return None

    columns = collated.columns

    fraccols,colnames = [],[]
    for c in columns:
        if c[-4:] == 'frac':
            fraccols.append(c)
        if c[0] == 't' and is_number(c[1:3]):
            colnames.append(c[:3])

    collist = list(set(colnames))
    collist.sort()

    # Plot distribution of vote fractions for each task

    ntasks = len(collist)
    ncols = 4 if ntasks > 9 else int(np.sqrt(ntasks))
    nrows = int(ntasks / ncols) if ntasks % ncols == 0 else int(ntasks / ncols) + 1

    sd = survey_dict()[survey]
    survey_name = sd['name']

    def f7(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))] 

    tasklabels = f7([re.split("[ax][0-9]",f)[0][11:-1] for f in fraccols])
    labels = [re.split("[ax][0-9]",f)[-1][1:-5] for f in fraccols]

    # Make pie charts of the plurality votes

    votearr = np.array(collated[fraccols])
    class_arr,task_arr,task_ans = [],[],[]
    for v in votearr:
        e,a = plurality(v,survey) 
        task_arr.append(e)
        task_ans.append(a)

    task_arr = np.array(task_arr)
    task_ans = np.array(task_ans)

    fig,axarr = plt.subplots(nrows=nrows,ncols=ncols,figsize=(15,12))

    colors=[u'#377EB8', u'#E41A1C', u'#4DAF4A', u'#984EA3', u'#FF7F00',u'#A6761D',u'#1B9E77']

    n = (task_arr.shape)[1]
    for i in range(n):
        ax = axarr.ravel()[i]
        c = Counter(task_ans[:,i][task_arr[:,i] == True])
        pv,pl = [],[]
        for k in c:
            pv.append(c[k])
            pl.append(labels[k])
        ax.pie(pv,labels=pl,colors=colors,autopct=lambda(p): '{:.0f}'.format(p * sum(pv) / 100))
        title = '{0:} - t{1:02} {2:}'.format(survey_name,i,tasklabels[i]) if i == 0 else 't{0:02} {1:}'.format(i,tasklabels[i])
        ax.set_title(title)
        ax.set_aspect('equal')

    # Remove empty axes from subplots
    if axarr.size > ntasks:
        for i in range(axarr.size - ntasks):
            ax = axarr.ravel()[axarr.size-(i+1)]
            ax.set_axis_off()

    fig.set_tight_layout(True)
    plt.savefig('{1}/decals/plots/pie_{0:}.eps'.format(survey,gzpath))
    plt.close()

    return None

def morph_table_gz2():

    overlap = True
    survey = 'decals'

    # Get weights
    try:
        fitsfile = "{0}/dr10/dr10_gz2_main_specz.fits".format(gzpath)
        hdr = fits.getheader(fitsfile,1)
        colnames = []
        for i in range(hdr['TFIELDS']):
            colnames.append(hdr['TTYPE{0}'.format(i+1)])

        if overlap:
            if survey == 'gz2':
                collation_file = "{0}/decals/csv/decals_gz2_main.csv".format(gzpath)
            elif survey == 'stripe82':
                collation_file = "{0}/decals/csv/decals_gz2_stripe82c1.csv".format(gzpath)
            elif survey == 'decals':
                collation_file = "{0}/decals/csv/decals_gz2_union.csv".format(gzpath)
            collated = pd.read_csv(collation_file)
        else:
            if survey == 'gz2':
                collation_file = "{0}/dr10/dr10_gz2_main_specz.csv".format(gzpath)
            elif survey == 'stripe82':
                collation_file = "{0}/dr10/dr10_gz2_stripe82_coadd1.csv".format(gzpath)
            collated = pd.read_csv(collation_file,names=colnames)
    except IOError:
        print "Collation file for {0:} does not exist. Aborting.".format(survey)
        return None

    columns = collated.columns

    fraccols,colnames = [],[]
    if survey == 'decals':
        for c in columns:
            if len(c) > 10:
                if c[-4:] == 'frac' and c[:6] == 'decals':
                    fraccols.append(c)
                if c[7] == 't' and is_number(c[8:10]):
                    colnames.append(c[7:10])
    else:
        for c in columns:
            if c[-17:] == 'weighted_fraction':
                fraccols.append(c)
            if c[0] == 't' and is_number(c[1:3]):
                colnames.append(c[:3])

    collist = list(set(colnames))
    collist.sort()

    # Plot distribution of vote fractions for each task

    ntasks = len(collist)
    ncols = 4 if ntasks > 9 else int(np.sqrt(ntasks))
    nrows = int(ntasks / ncols) if ntasks % ncols == 0 else int(ntasks / ncols) + 1

    sd = survey_dict()[survey]
    survey_name = sd['name']

    def f7(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))] 

    if survey == 'decals':
        tasklabels = f7([re.split("[ax][0-9]",f)[0][11:-1] for f in fraccols])
        labels = [re.split("[ax][0-9]",f)[-1][1:-5] for f in fraccols]
    else:
        tasklabels = f7([re.split("[ax][0-9]",f)[0][4:-1] for f in fraccols])
        labels = [re.split("[ax][0-9]",f[4:-18])[-1][2:] for f in fraccols]

    # Make pie charts of the plurality votes

    votearr = np.array(collated[fraccols])
    class_arr,task_arr,task_ans = [],[],[]
    for v in votearr:
        e,a = plurality(v,survey) 
        task_arr.append(e)
        task_ans.append(a)

    task_arr = np.array(task_arr)
    task_ans = np.array(task_ans)

    fig,axarr = plt.subplots(nrows=nrows,ncols=ncols,figsize=(15,12))

    colors=[u'#377EB8', u'#E41A1C', u'#4DAF4A', u'#984EA3', u'#FF7F00',u'#A6761D',u'#1B9E77']

    n = (task_arr.shape)[1]
    for i in range(n):
        ax = axarr.ravel()[i]
        c = Counter(task_ans[:,i][task_arr[:,i] == True])
        pv,pl = [],[]
        task_total = sum(c.values())
        for k in c:
            pv.append(c[k])
            pl.append(labels[k])

            # Print to screen in LaTeX format
            print "{0:20} & {1:6} & {3:.2f} & {2:.2f}".format(labels[k],c[k],c[k] * 1./task_total,c[k] * 1./len(collated))
        print ""
        ax.pie(pv,labels=pl,colors=colors,autopct=lambda(p): '{:.0f}'.format(p * sum(pv) / 100))
        title = '{0:} - t{1:02} {2:}'.format(survey_name,i,tasklabels[i]) if i == 0 else 't{0:02} {1:}'.format(i,tasklabels[i])
        ax.set_title(title)
        ax.set_aspect('equal')

    # Remove empty axes from subplots
    if axarr.size > ntasks:
        for i in range(axarr.size - ntasks):
            ax = axarr.ravel()[axarr.size-(i+1)]
            ax.set_axis_off()

    fig.set_tight_layout(True)
    suffix = '_overlap' if overlap else ''
    plt.savefig('{1}/decals/plots/pie_{0}{2}.eps'.format(survey,gzpath,suffix))
    plt.close()

    return None


if __name__ == "__main__":

    #morph_table_gz2()
    morphology_distribution()
