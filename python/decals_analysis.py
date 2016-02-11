# Do some preliminary analysis on the results of the DECaLS-Galaxy Zoo data. 

from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
from astropy.io import fits
from astropy.cosmology import WMAP9

decals_path = '/Users/willettk/Astronomy/Research/GalaxyZoo/decals'
plot_path = '{0}/plots'.format(decals_path)

def load_data():

    mgs = fits.getdata('{0}/matched/gz2_main.fits'.format(decals_path),1)
    s82 = fits.getdata('{0}/matched/gz2_s82_coadd1.fits'.format(decals_path),1)
    decals = fits.getdata('{0}/matched/decals_dr1.fits'.format(decals_path),1)

    return mgs,s82,decals

def color_mag_plots(mgs,s82,decals,savefig=False):

    # Make paneled histograms of the color distribution for several magnitude bins of three sets of data:
    #   SDSS main sample (GZ2)
    #   Stripe 82 coadded data (GZ2)
    #   DECaLS


    redshifts = (0.12,0.08,0.05)
    appmag_lim = 17.0
    # Work out the magnitude limit from cosmology

    fig,axarr = plt.subplots(num=1,nrows=3,ncols=3,figsize=(12,10))

    for z,ax in zip(redshifts,axarr.ravel()):
        absmag_lim = appmag_lim - WMAP9.distmod(z).value
        maglim = (mgs['PETROMAG_MR'] < absmag_lim) & (mgs['REDSHIFT'] <= z)
        spiral = mgs['t01_smooth_or_features_a02_features_or_disk_weighted_fraction'] >= 0.8
        elliptical = mgs['t01_smooth_or_features_a01_smooth_weighted_fraction'] >= 0.8
        ax.hist(mgs[maglim & spiral]['PETROMAG_U'] - mgs[maglim & spiral]['PETROMAG_R'],range=(0,4),bins=25,color='blue',histtype='step',label='spiral')
        ax.hist(mgs[maglim & elliptical]['PETROMAG_U'] - mgs[maglim & elliptical]['PETROMAG_R'],range=(0,4),bins=25,color='red',histtype='step',label='elliptical')

        ax.set_xlabel(r'$(u-r)$',fontsize=16)
        ax.set_title(r'$M_r<{0:.2f}, z<{1:.2f}$'.format(absmag_lim,z),fontsize=16)
        ax.text(0.95,0.95,'MGS',ha='right',va='top',transform=ax.transAxes)

        if ax == axarr.ravel()[0]:
            ax.legend(loc='upper left',fontsize=10)

    s82_lim = 17.77
    for z,ax in zip(redshifts,axarr.ravel()[3:6]):
        absmag_lim = s82_lim - WMAP9.distmod(z).value
        maglim = (s82['PETROMAG_MR'] < absmag_lim) & (s82['REDSHIFT'] <= z)
        spiral = s82['t01_smooth_or_features_a02_features_or_disk_weighted_fraction'] >= 0.8
        elliptical = s82['t01_smooth_or_features_a01_smooth_weighted_fraction'] >= 0.8
        ax.hist(s82[maglim & spiral]['PETROMAG_U'] - s82[maglim & spiral]['PETROMAG_R'],range=(0,4),bins=25,color='blue',histtype='step',label='spiral')
        ax.hist(s82[maglim & elliptical]['PETROMAG_U'] - s82[maglim & elliptical]['PETROMAG_R'],range=(0,4),bins=25,color='red',histtype='step',label='elliptical')

        ax.set_xlabel(r'$(u-r)$',fontsize=16)
        ax.set_title(r'$M_r<{0:.2f}, z<{1:.2f}$'.format(absmag_lim,z),fontsize=16)
        ax.text(0.95,0.95,'Stripe 82',ha='right',va='top',transform=ax.transAxes)

    decals_lim = 17.77
    for z,ax in zip(redshifts,axarr.ravel()[6:]):
        absmag_lim = decals_lim - WMAP9.distmod(z).value
        maglim = (decals['metadata.mag.abs_r'] < absmag_lim) & (decals['metadata.redshift'] <= z)
        spiral = decals['t00_smooth_or_features_a1_features_frac'] >= 0.8
        elliptical = decals['t00_smooth_or_features_a0_smooth_frac'] >= 0.8
        ax.hist(decals[maglim & spiral]['metadata.mag.u'] - decals[maglim & spiral]['metadata.mag.r'],range=(0,4),bins=25,color='blue',histtype='step',label='spiral')
        ax.hist(decals[maglim & elliptical]['metadata.mag.u'] - decals[maglim & elliptical]['metadata.mag.r'],range=(0,4),bins=25,color='red',histtype='step',label='elliptical')

        ax.set_xlabel(r'$(u-r)$',fontsize=16)
        ax.set_title(r'$M_r<{0:.2f}, z<{1:.2f}$'.format(absmag_lim,z),fontsize=16)
        ax.text(0.95,0.95,'DECaLS',ha='right',va='top',transform=ax.transAxes)

    fig.tight_layout()
    if savefig:
        plt.savefig('{0}/color_hist.pdf'.format(plot_path))
    else:
        plt.show()

    return None

def color_mag_ratio(mgs,s82,decal,savefig=False):

    # Plot the spiral to elliptical ratio as a function of optical color.

    redshifts = (0.12,0.08,0.05)
    linestyles = ('solid','dashed','dashdot')

    datasets = ({'data':mgs,
                 'title':'MGS',
                 'appmag':17.0,
                 'sp':'t01_smooth_or_features_a02_features_or_disk_weighted_fraction',
                 'el':'t01_smooth_or_features_a01_smooth_weighted_fraction',
                 'umag':'PETROMAG_U',
                 'rmag':'PETROMAG_R',
                 'absr':'PETROMAG_MR',
                 'redshift':'REDSHIFT'},
                {'data':s82,
                 'title':'Stripe 82',
                 'appmag':17.77,
                 'sp':'t01_smooth_or_features_a02_features_or_disk_weighted_fraction',
                 'el':'t01_smooth_or_features_a01_smooth_weighted_fraction',
                 'umag':'PETROMAG_U',
                 'rmag':'PETROMAG_R',
                 'absr':'PETROMAG_MR',
                 'redshift':'REDSHIFT'},
                {'data':decals,
                 'title':'DECaLS',
                 'appmag':17.77,
                 'sp':'t00_smooth_or_features_a1_features_frac',
                 'el':'t00_smooth_or_features_a0_smooth_frac',
                 'umag':'metadata.mag.u',
                 'rmag':'metadata.mag.r',
                 'absr':'metadata.mag.abs_r',
                 'redshift':'metadata.redshift'})

    # Work out the magnitude limit from cosmology

    fig,axarr = plt.subplots(num=2,nrows=1,ncols=3,figsize=(12,5))

    for ax,d in zip(axarr.ravel(),datasets):
        for z,ls in zip(redshifts,linestyles):
            absmag_lim = d['appmag'] - WMAP9.distmod(z).value
            maglim = (d['data'][d['absr']] < absmag_lim) & (d['data'][d['redshift']] <= z)
            spiral = d['data'][d['sp']] >= 0.8
            elliptical = d['data'][d['el']] >= 0.8
            n_sp,bins_sp = np.histogram(d['data'][maglim & spiral][d['umag']] - d['data'][maglim & spiral][d['rmag']],range=(0,4),bins=25)
            n_el,bins_el = np.histogram(d['data'][maglim & elliptical][d['umag']] - d['data'][maglim & elliptical][d['rmag']],range=(0,4),bins=25)

            plotval = np.log10(n_sp * 1./n_el)
            ax.plot(bins_sp[1:],plotval,linestyle=ls,label=r'$M_r<{0:.2f}, z<{1:.2f}$'.format(absmag_lim,z))

        ax.set_xlabel(r'$(u-r)$',fontsize=16)
        ax.set_ylabel(r'$\log(n_{sp}/n_{el})$',fontsize=16)
        ax.set_ylim(-1.5,1.5)
        ax.set_title(d['title'],fontsize=16)
        if ax == axarr.ravel()[0]:
            ax.legend(loc='upper left',fontsize=8)

    fig.tight_layout()
    
    if savefig:
        plt.savefig('{0}/feature_ratio.pdf'.format(plot_path))
    else:
        plt.show()

    return None

def feature_comparison():

    # Plot the difference in vote fractions for the matched galaxies

    filename = '{0}/fits/decals_gz2_union.fits'.format(decals_path)

    data = fits.getdata(filename,1)

    # Map the columns
    columns = [{'gz2':"gz2_t01_smooth_or_features_a01_smooth_fraction",  "decals":"decals_t00_smooth_or_features_a0_smooth_frac"},
                {'gz2':"gz2_t01_smooth_or_features_a02_features_or_disk_fraction",   "decals":"decals_t00_smooth_or_features_a1_features_frac"},
                {'gz2':"gz2_t01_smooth_or_features_a03_star_or_artifact_fraction",   "decals":"decals_t00_smooth_or_features_a2_artifact_frac"},
                {'gz2':"gz2_t02_edgeon_a04_yes_fraction",                            "decals":"decals_t01_disk_edge_on_a0_yes_frac"},
                {'gz2':"gz2_t02_edgeon_a05_no_fraction",                             "decals":"decals_t01_disk_edge_on_a1_no_frac"},
                {'gz2':"gz2_t03_bar_a06_bar_fraction",                               "decals":"decals_t02_bar_a0_bar_frac"},
                {'gz2':"gz2_t03_bar_a07_no_bar_fraction",                            "decals":"decals_t02_bar_a1_no_bar_frac"},
                {'gz2':"gz2_t04_spiral_a08_spiral_fraction",                         "decals":"decals_t03_spiral_a0_spiral_frac"},
                {'gz2':"gz2_t04_spiral_a09_no_spiral_fraction",                      "decals":"decals_t03_spiral_a1_no_spiral_frac"},
                {'gz2':"gz2_t05_bulge_prominence_a10_no_bulge_fraction",             "decals":"decals_t04_bulge_prominence_a0_no_bulge_frac"},
                {'gz2':"gz2_t05_bulge_prominence_a11_just_noticeable_fraction",      "decals":"decals_t04_bulge_prominence_a1_obvious_frac"},
                {'gz2':"gz2_t05_bulge_prominence_a12_obvious_fraction",              "decals":"decals_t04_bulge_prominence_a2_dominant_frac"},
                {'gz2':"gz2_t07_rounded_a16_completely_round_fraction",              "decals":"decals_t08_rounded_a0_completely_round_frac"},
                {'gz2':"gz2_t07_rounded_a17_in_between_fraction",                    "decals":"decals_t08_rounded_a1_in_between_frac"},
                {'gz2':"gz2_t07_rounded_a18_cigar_shaped_fraction",                  "decals":"decals_t08_rounded_a2_cigar_shaped_frac"},
                {'gz2':"gz2_t08_odd_feature_a19_ring_fraction",                      "decals":"decals_t10_odd_feature_x1_ring_frac"},
                {'gz2':"gz2_t08_odd_feature_a20_lens_or_arc_fraction",               "decals":"decals_t10_odd_feature_x2_lens_frac"},
                {'gz2':"gz2_t08_odd_feature_a22_irregular_fraction",                 "decals":"decals_t10_odd_feature_x4_irregular_frac"},
                {'gz2':"gz2_t08_odd_feature_a23_other_fraction",                     "decals":"decals_t10_odd_feature_x5_other_frac"},
                {'gz2':"gz2_t08_odd_feature_a38_dust_lane_fraction",                 "decals":"decals_t10_odd_feature_x3_dustlane_frac"},
                {'gz2':"gz2_t09_bulge_shape_a25_rounded_fraction",                   "decals":"decals_t07_bulge_shape_a0_rounded_frac"},
                {'gz2':"gz2_t09_bulge_shape_a26_boxy_fraction",                      "decals":"decals_t07_bulge_shape_a1_boxy_frac"},
                {'gz2':"gz2_t09_bulge_shape_a27_no_bulge_fraction",                  "decals":"decals_t07_bulge_shape_a2_no_bulge_frac"},
                {'gz2':"gz2_t10_arms_winding_a28_tight_fraction",                    "decals":"decals_t05_arms_winding_a0_tight_frac"},
                {'gz2':"gz2_t10_arms_winding_a29_medium_fraction",                   "decals":"decals_t05_arms_winding_a1_medium_frac"},
                {'gz2':"gz2_t10_arms_winding_a30_loose_fraction",                    "decals":"decals_t05_arms_winding_a2_loose_frac"},
                {'gz2':"gz2_t11_arms_number_a31_1_fraction",                         "decals":"decals_t06_arms_number_a0_1_frac"},
                {'gz2':"gz2_t11_arms_number_a32_2_fraction",                         "decals":"decals_t06_arms_number_a1_2_frac"},
                {'gz2':"gz2_t11_arms_number_a33_3_fraction",                         "decals":"decals_t06_arms_number_a2_3_frac"},
                {'gz2':"gz2_t11_arms_number_a34_4_fraction",                         "decals":"decals_t06_arms_number_a3_4_frac"},
                {'gz2':"gz2_t11_arms_number_a36_more_than_4_fraction",               "decals":"decals_t06_arms_number_a4_more_than_4_frac"}]

    # Working, but still needs to sort for questions that are ACTUALLY ANSWERED. Lots of pileup at 0,0.

    fig,axarr = plt.subplots(num=1,nrows=4,ncols=8,figsize=(12,8))

    for i,ax in enumerate(axarr.ravel()):
        try:
            col = columns[i]
            ax.hist2d(data[col['gz2']],data[col['decals']],bins=(20,20),range=[[0,1],[0,1]],normed=True,cmap = cm.hot)
            #ax.set_xlabel(r"$f_{GZ2}$",fontsize=20)
            #ax.set_ylabel(r"$f_{DECaLS}$",fontsize=20)
        except IndexError:
            break

                
    fig.tight_layout()

    plt.show()

    return None


if __name__ == "__main__":

    #mgs,s82,decals = load_data()
    #color_mag_plots(mgs,s82,decals,savefig=True)
    #color_mag_ratio(mgs,s82,decals,savefig=True)
    feature_comparison()
