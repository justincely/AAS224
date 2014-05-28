

    import os
    import glob
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from astropy.io import fits
    
    import lightcurve


    #-- iPython notebook magic command to make plots appear within the 
    #-- browser instead of opening a separate plotting window.
    %matplotlib inline


    mpl.rcParams['font.size'] = 16


    os.environ['lref'] = '/grp/hst/cdbs/lref/'


    fig = plt.figure(figsize=(13, 7))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    
    hdu = fits.open('lightcurves/lc1va0zgq_x1d.fits')
    ax1.plot(hdu[1].data['wavelength'].ravel(), hdu[1].data['flux'].ravel())
    ax1.set_ylim(0, 3e-14)
    ax1.set_xlim(1000, 2200)
    ax1.set_xlabel('Wavelength $(\AA)$')
    ax1.set_ylabel('Flux')
    
    lc = lightcurve.open(filename='lightcurves/lc1va0zgq_corrtag_a.fits')
    ax2.plot(lc.times, lc.counts, 'o')
    ax2.set_ylabel('Counts')
    ax2.set_xlabel('Time (s)')
    
    fig.tight_layout()  
    fig.savefig('spectrum_vs_lightcurve.pdf', bbox_inches='tight')

    Extracting from: ['lightcurves/lc1va0zgq_corrtag_a.fits']



![png](lightcurve_tutorial_files/lightcurve_tutorial_4_1.png)



    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(1, 1, 1)
    
    datasets = glob.glob('lightcurves/*x1d.fits')
    datasets.sort()
    
    for item in datasets:
        hdu = fits.open(item)
        ax.plot(hdu[1].data['wavelength'].ravel(), hdu[1].data['flux'].ravel())
        
    ax.set_ylim(0, 3e-14)
    ax.set_xlim(1000, 2200)




    (1000, 2200)




![png](lightcurve_tutorial_files/lightcurve_tutorial_5_1.png)



    fig, (ax1, ax2) = plt.subplots(2, figsize=(15, 7), sharex=True)
    
    datasets = glob.glob('lightcurves/*corrtag_a*')
    datasets.sort()
    
    for item in datasets:
        #-- Minimal wavelength screening
        lc = lightcurve.open(filename=item, step=10)
        ax1.plot(lc.mjd, lc.counts, 'o')
        
        #-- Specifically an area with no airglow
        lc.extract(wlim=(1600, 1800), step=10)
        ax2.plot(lc.mjd, lc.counts, 'o')
        
    ax1.set_xlabel('MJD')
    ax1.set_ylabel('Counts')
    
    ax2.set_xlabel('MJD')
    ax2.set_ylabel('Counts')
        
    fig.subplots_adjust(hspace=0)
    fig.savefig('eclipse.pdf', bbox_inches='tight')

    Extracting from: ['lightcurves/lc1va0yeq_corrtag_a.fits']
    Extracting from: ['lightcurves/lc1va0yjq_corrtag_a.fits']
    Extracting from: ['lightcurves/lc1va0zcq_corrtag_a.fits']
    Extracting from: ['lightcurves/lc1va0zgq_corrtag_a.fits']



![png](lightcurve_tutorial_files/lightcurve_tutorial_6_1.png)



    
    dataset = 'lightcurves/lc1va0zgq_corrtag_a.fits'
    #dataset = 'lightcurves/lc1va0yeq_corrtag_a.fits'
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(10,15), sharex=True)
    
    lc = lightcurve.open(filename=dataset, step=1)
    ax1.plot(lc.times, lc.counts, 'o', label='Stepsize=1s')
    ax1.legend(shadow=True, numpoints=1, loc='best', fontsize=14)
    ax1.set_ylabel('Counts')
    
    
    #-- Re-extract with different stepsize
    lc.extract(step=15)
    ax2.plot(lc.times, lc.counts, 'o', label='Stepsize=15s')   
    ax2.legend(shadow=True, numpoints=1, loc='best', fontsize=14)
    ax2.set_ylabel('Counts')
    
    #-- A select wavelength or two   
    lc.extract(step=15, wlim=(1400, 1500))
    ax3.plot(lc.times, lc.counts, 'o', label='($1400\AA - 1500\AA)$')
    
    lc.extract(step=15, wlim=(1500, 1600))
    ax3.plot(lc.times, lc.counts, 'o', label='($1500\AA - 1600\AA)$')
    
    ax3.legend(shadow=True, numpoints=1, loc='best', fontsize=14)
    ax3.set_ylabel('Counts')
    ax3.set_xlabel('Time (s)')
    
    fig.subplots_adjust(hspace=0)
    #fig.tight_layout() 
    fig.savefig('eclipse_things.pdf', bbox_inches='tight')

    Extracting from: ['lightcurves/lc1va0zgq_corrtag_a.fits']



![png](lightcurve_tutorial_files/lightcurve_tutorial_7_1.png)



    
