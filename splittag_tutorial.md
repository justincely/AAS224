

    import os
    import glob
    import numpy as np
    from scipy.ndimage.filters import convolve
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from astropy.io import fits
    
    import lightcurve
    
    from calcos import calcos
    from costools import splittag


    #-- iPython notebook magic command to make plots appear within the 
    #-- browser instead of opening a separate plotting window.
    %matplotlib inline


    mpl.rcParams['font.size'] = 16

### Transit spectroscopy


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



![png](splittag_tutorial_files/splittag_tutorial_4_1.png)



    dataset = 'lightcurves/lc1va0zgq_corrtag_a.fits'
    file_path, file_name = os.path.split(dataset)
    splittag.splittag(dataset, 
                      os.path.join('lightcurves', file_name[:9]),
                      time_list="0,400,900,1300")

    lightcurves/lc1va0zgq_1_corrtag_a.fits written
    lightcurves/lc1va0zgq_2_corrtag_a.fits written
    lightcurves/lc1va0zgq_3_corrtag_a.fits written



    dataset_list = glob.glob('lightcurves/*_?_corrtag_a.fits')
    
    if not os.path.exists('lightcurves/out/'):
        os.mkdir('lightcurves/out/')
        
    for item in dataset_list:
        calcos(item, outdir='lightcurves/out/')

    CALCOS version 2.21 (2014-02-18)
    numpy version 1.8.0
    astropy version 0.3
    Begin 28-May-2014 00:59:12 EDT
    Input file = lightcurves/lc1va0zgq_1_corrtag_a.fits
    
    TIME-TAG calibration -- 28-May-2014 00:59:14 EDT
    Input     lightcurves/lc1va0zgq_1_corrtag_a.fits
    OutTag    lightcurves/out/lc1va0zgq_1_corrtag_a.fits
    OutFlt    lightcurves/out/lc1va0zgq_1_flt_a.fits
    OutCounts lightcurves/out/lc1va0zgq_1_counts_a.fits
    OutFlash  lightcurves/out/lc1va0zgq_1_lampflash_a.fits
    DETECTOR  FUV, segment A
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G140L, CENWAVE 1105, FPOFFSET 1
    APERTURE  PSA
    
    HVTAB   = lref$wa41551al_hv.fits
    BADTCORR  OMIT (already complete)
    RANDCORR  OMIT (already complete)
    TEMPCORR  OMIT (already complete)
    GEOCORR   OMIT (already complete)
    WALKCORR  OMIT (already complete)
    DEADCORR  OMIT (already complete)
    PHACORR   OMIT (already complete)
    DOPPCORR  OMIT (already complete)
    XTRACTAB= lref$x1v17418l_1dx.fits
    DISPTAB = lref$x1v17415l_disp.fits
    BRFTAB  = lref$x1u1459il_brf.fits
    FLATCORR  OMIT
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$x1u1459hl_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Airglow region for LY_ALPHA is X: 2283 to 2468, Y: 509 to 570
    Airglow region for OI_1304 is X: 3370 to 3603, Y: 509 to 570
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 28-May-2014 01:00:10 EDT
    Input     lightcurves/out/lc1va0zgq_1_flt_a.fits
    Incounts  lightcurves/out/lc1va0zgq_1_counts_a.fits
    Output    lightcurves/out/lc1va0zgq_1_x1d_a.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$x1v17418l_1dx.fits
    DISPTAB = lref$x1v17415l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$x1v17416l_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVA spectrum was found at y = 536.43 vs. nominal y = 540.31
        error estimate for y location = 5.16, FWHM = 5.51
    Spectrum will be extracted at y = 540.31
    SPWCSTAB= lref$x1v17417l_spwcs.fits
    End   28-May-2014 01:00:21 EDT
    CALCOS version 2.21 (2014-02-18)
    numpy version 1.8.0
    astropy version 0.3
    Begin 28-May-2014 01:00:21 EDT
    Input file = lightcurves/lc1va0zgq_2_corrtag_a.fits
    
    TIME-TAG calibration -- 28-May-2014 01:00:23 EDT
    Input     lightcurves/lc1va0zgq_2_corrtag_a.fits
    OutTag    lightcurves/out/lc1va0zgq_2_corrtag_a.fits
    OutFlt    lightcurves/out/lc1va0zgq_2_flt_a.fits
    OutCounts lightcurves/out/lc1va0zgq_2_counts_a.fits
    OutFlash  lightcurves/out/lc1va0zgq_2_lampflash_a.fits
    DETECTOR  FUV, segment A
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G140L, CENWAVE 1105, FPOFFSET 1
    APERTURE  PSA
    
    HVTAB   = lref$wa41551al_hv.fits
    BADTCORR  OMIT (already complete)
    RANDCORR  OMIT (already complete)
    TEMPCORR  OMIT (already complete)
    GEOCORR   OMIT (already complete)
    WALKCORR  OMIT (already complete)
    DEADCORR  OMIT (already complete)
    PHACORR   OMIT (already complete)
    DOPPCORR  OMIT (already complete)
    XTRACTAB= lref$x1v17418l_1dx.fits
    DISPTAB = lref$x1v17415l_disp.fits
    BRFTAB  = lref$x1u1459il_brf.fits
    FLATCORR  OMIT
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$x1u1459hl_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Airglow region for LY_ALPHA is X: 2283 to 2468, Y: 509 to 570
    Airglow region for OI_1304 is X: 3370 to 3603, Y: 509 to 570
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 28-May-2014 01:01:17 EDT
    Input     lightcurves/out/lc1va0zgq_2_flt_a.fits
    Incounts  lightcurves/out/lc1va0zgq_2_counts_a.fits
    Output    lightcurves/out/lc1va0zgq_2_x1d_a.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$x1v17418l_1dx.fits
    DISPTAB = lref$x1v17415l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$x1v17416l_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVA spectrum was found at y = 536.49 vs. nominal y = 540.31
        error estimate for y location = 13.70, FWHM = 5.36
    Spectrum will be extracted at y = 540.31
    SPWCSTAB= lref$x1v17417l_spwcs.fits
    End   28-May-2014 01:01:27 EDT
    CALCOS version 2.21 (2014-02-18)
    numpy version 1.8.0
    astropy version 0.3
    Begin 28-May-2014 01:01:27 EDT
    Input file = lightcurves/lc1va0zgq_3_corrtag_a.fits
    
    TIME-TAG calibration -- 28-May-2014 01:01:29 EDT
    Input     lightcurves/lc1va0zgq_3_corrtag_a.fits
    OutTag    lightcurves/out/lc1va0zgq_3_corrtag_a.fits
    OutFlt    lightcurves/out/lc1va0zgq_3_flt_a.fits
    OutCounts lightcurves/out/lc1va0zgq_3_counts_a.fits
    OutFlash  lightcurves/out/lc1va0zgq_3_lampflash_a.fits
    DETECTOR  FUV, segment A
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G140L, CENWAVE 1105, FPOFFSET 1
    APERTURE  PSA
    
    HVTAB   = lref$wa41551al_hv.fits
    BADTCORR  OMIT (already complete)
    RANDCORR  OMIT (already complete)
    TEMPCORR  OMIT (already complete)
    GEOCORR   OMIT (already complete)
    WALKCORR  OMIT (already complete)
    DEADCORR  OMIT (already complete)
    PHACORR   OMIT (already complete)
    DOPPCORR  OMIT (already complete)
    XTRACTAB= lref$x1v17418l_1dx.fits
    DISPTAB = lref$x1v17415l_disp.fits
    BRFTAB  = lref$x1u1459il_brf.fits
    FLATCORR  OMIT
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$x1u1459hl_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Airglow region for LY_ALPHA is X: 2283 to 2468, Y: 509 to 570
    Airglow region for OI_1304 is X: 3370 to 3603, Y: 509 to 570
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 28-May-2014 01:02:04 EDT
    Input     lightcurves/out/lc1va0zgq_3_flt_a.fits
    Incounts  lightcurves/out/lc1va0zgq_3_counts_a.fits
    Output    lightcurves/out/lc1va0zgq_3_x1d_a.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$x1v17418l_1dx.fits
    DISPTAB = lref$x1v17415l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$x1v17416l_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVA spectrum was found at y = 536.52 vs. nominal y = 540.31
        error estimate for y location = 7.53, FWHM = 5.26
    Spectrum will be extracted at y = 540.31
    SPWCSTAB= lref$x1v17417l_spwcs.fits
    End   28-May-2014 01:02:15 EDT



    boxsize=15
    
    fig, (ax1, ax2) = plt.subplots(2, figsize=(15,10), sharex=True)
    
    hdu = fits.open('lightcurves/lc1va0zgq_x1d.fits')
    ax1.plot(hdu[1].data['wavelength'].ravel(), convolve(hdu[1].data['flux'].ravel(), np.ones(boxsize)/boxsize), color='k', label='Archive Dataset')
    ax1.set_ylim(0, 1e-14)
    ax1.set_xlim(1100, 1650)
    ax1.set_title('lc1va0zgq_x1d.fits')
    ax1.set_ylabel('Flux')
    
    leg = ax1.legend(shadow=True, numpoints=1, markerscale=3)
    # set the linewidth of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)
    
    
    before = fits.open('lightcurves/out/lc1va0zgq_1_x1d.fits')
    during = fits.open('lightcurves/out/lc1va0zgq_2_x1d.fits')
    after = fits.open('lightcurves/out/lc1va0zgq_3_x1d.fits')
    ax2.plot(before[1].data['wavelength'].ravel(), convolve(before[1].data['flux'].ravel(), np.ones(boxsize)/boxsize), label='Before')
    ax2.plot(during[1].data['wavelength'].ravel(), convolve(during[1].data['flux'].ravel(), np.ones(boxsize)/boxsize), label='During')
    ax2.plot(after[1].data['wavelength'].ravel(), convolve(after[1].data['flux'].ravel(), np.ones(boxsize)/boxsize), label='After')
    
    ax2.set_ylim(0, 1e-14)
    ax2.set_xlim(1100, 1650)
    ax2.set_xlabel('Wavelength $(\AA)$')
    ax2.set_ylabel('Flux')
    leg = ax2.legend(shadow=True, numpoints=1, markerscale=3)
    # set the linewidth of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)
    
    fig.subplots_adjust(hspace=0)
    fig.savefig('eclipse_split.pdf', bbox_inches='tight')


![png](splittag_tutorial_files/splittag_tutorial_7_0.png)



    
