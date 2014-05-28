

    import os
    import glob
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from astropy.io import fits


    #-- iPython notebook magic command to make plots appear within the 
    #-- browser instead of opening a separate plotting window.
    %matplotlib inline


    mpl.rcParams['font.size'] = 16


    fig, (ax1, ax2) = plt.subplots(2, figsize=(7, 7), sharex=True)
    
    hdu = fits.open('airglow/lbry01ilq_corrtag_a.fits')
    ax1.hist(hdu['events'].data['PHA'], bins=range(0, 32))
    ax1.set_ylabel('FUVA Distribution')
    ax1.axvspan(xmin=2, xmax=3, ymin=0, ymax=1, color='grey', alpha=.5)
    ax1.axvspan(xmin=13, xmax=23, ymin=0, ymax=1, color='grey', alpha=.5)
    
    
    hdu = fits.open('airglow/lbry01ilq_corrtag_b.fits')
    ax2.hist(hdu['events'].data['PHA'], bins=range(0, 32))
    ax2.set_xlabel('PHA')
    ax2.set_ylabel('FUVB Distribution')
    ax2.axvspan(xmin=12, xmax=23, ymin=0, ymax=1, color='grey', alpha=.5)
    
    for axis in [ax1, ax2]:
        axis.axvline(x=2, color='r', ls='--', lw=3, label='Lower Bound')
        axis.axvline(x=23, color='r', ls='--', lw=3, label='Upper Bound')
    
    fig.subplots_adjust(hspace=0)
    fig.savefig('filtering_by_segment.pdf', bbox_inches='tight')


![png](filtering_tutorial_files/filtering_tutorial_3_0.png)



    
