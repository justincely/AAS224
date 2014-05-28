

    import os
    import glob
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from astropy.io import fits
    
    from calcos import calcos
    from costools import timefilter


    #-- iPython notebook magic command to make plots appear within the 
    #-- browser instead of opening a separate plotting window.
    %matplotlib inline


    mpl.rcParams['font.size'] = 16

### Understanding Airglow


    fig = plt.figure(figsize=(15, 15))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    
    ax1.set_ylabel('FUVA', fontsize=16)
    ax2.set_ylabel('FUVB', fontsize=16)
    
    for item in glob.glob('airglow/lbh*x1dsum.fits'):
        hdu = fits.open(item)
        ax1.plot(hdu[1].data['WAVELENGTH'].ravel(), hdu[1].data['FLUX'].ravel(), label=os.path.split(item)[1][:9])
        ax2.plot(hdu[1].data['WAVELENGTH'].ravel(), hdu[1].data['FLUX'].ravel(), label=os.path.split(item)[1][:9])
    
    for axis in [ax1, ax2]:
        axis.grid(True)
        axis.set_xlabel("Wavelength", fontsize=16)
        axis.legend(shadow=True, numpoints=1, loc='upper left', fontsize=16)
        
    ax1.set_ylim(0, 5e-14)
    ax1.set_xlim(1295, 1315)
    
    ax2.set_ylim(0, 2e-14)
    ax2.set_xlim(1190, 1230)




    (1190, 1230)




![png](timefilter_tutorial_files/timefilter_tutorial_4_1.png)



    fig = plt.figure(figsize=(15, 15))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    
    ax1.set_ylabel('Lyman Alpha (cnts/s)', fontsize=16)
    ax2.set_ylabel('Oxygen I (cnts/s)', fontsize=16)
    
    
    for item in glob.glob('airglow/lbh*corrtag_b.fits'):
        hdu = fits.open(item)
        ax1.plot(hdu['timeline'].data['SUN_ALT'], hdu['timeline'].data['LY_ALPHA'], 'o', label=os.path.split(item)[1][:9])
    
    for item in glob.glob('airglow/lbh*corrtag_a.fits'):
        hdu = fits.open(item)
        ax2.plot(hdu['timeline'].data['SUN_ALT'], hdu['timeline'].data['OI_1304'],'o', label=os.path.split(item)[1][:9])
        
    for axis in [ax1, ax2]:
        axis.grid(True)
        axis.axvspan(xmin=-45, xmax=0, ymin=0, ymax=1, label='Orbital Night', alpha=.5, color='grey')
        axis.set_xlim(-45, 90)
        axis.set_xlabel("Sun altitude above horizon (degrees)", fontsize=16)
        axis.legend(shadow=True, numpoints=1, loc='upper left', fontsize=16)
        
    fig.savefig('airglow_strength.pdf', bbox_inches='tight')


![png](timefilter_tutorial_files/timefilter_tutorial_5_0.png)


### Dealing with airglow in data


    fig = plt.figure(figsize=(17, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    for dataset in glob.glob('airglow/lbry*x1d.fits'):
        hdu = fits.open(dataset)
        ax.plot(hdu[1].data['wavelength'].ravel(), hdu[1].data['flux'].ravel(), label=hdu[0].header['rootname'])
        
    ax.annotate('Geo-coronal Oxygen I', 
                xycoords='data',
                xy=(1302, 1e-13),
                xytext=(1310, 1.5e-13),
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=16)
    
    ax.annotate('', 
                xycoords='data',
                xy=(1305, 1e-13),
                xytext=(1310, 1.5e-13),
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=16)
        
    ax.set_ylim(0, 2e-13)
    ax.set_xlim(1280, 1340)
    ax.set_ylabel('Flux')
    ax.set_xlabel('Wavelength $(\AA)$')
    ax.legend(shadow=True, numpoints=1, fontsize=16)
    
    fig.savefig('spectrum_before.pdf', bbox_inches='tight')


![png](timefilter_tutorial_files/timefilter_tutorial_7_0.png)



    os.environ['lref'] = '/grp/hst/cdbs/lref/'


    if not os.path.exists('out_timefilter'):
        os.mkdir('out_timefilter')
    
    for dataset in glob.glob('airglow/lbry*corrtag*.fits'):
        filepath, filename = os.path.split(dataset)
        print "Filtering ", filename
        timefilter.TimelineFilter(input=dataset, 
                                  output='out_timefilter/' + filename,
                                  filter='SUN_ALT > 0')

    Filtering  lbry01i6q_corrtag_a.fits
    Filtering  lbry01i6q_corrtag_b.fits
    Filtering  lbry01icq_corrtag_a.fits
    Filtering  lbry01icq_corrtag_b.fits
    Filtering  lbry01ifq_corrtag_a.fits
    Filtering  lbry01ifq_corrtag_b.fits
    Filtering  lbry01ilq_corrtag_a.fits
    Filtering  lbry01ilq_corrtag_b.fits



    # This step will take awhile as we re-un CalCOS on the filtered datasets.
    
    for dataset in glob.glob('out_timefilter/lbry*corrtag_a*'):
        calcos(dataset, outdir='out_calcos/')

    CALCOS version 2.21 (2014-02-18)
    numpy version 1.8.0
    astropy version 0.3
    Begin 27-May-2014 21:55:27 EDT
    Input file = out_timefilter/lbry01i6q_corrtag_a.fits
    
    TIME-TAG calibration -- 27-May-2014 21:55:31 EDT
    Input     out_timefilter/lbry01i6q_corrtag_a.fits
    OutTag    out_calcos/lbry01i6q_corrtag_a.fits
    OutFlt    out_calcos/lbry01i6q_flt_a.fits
    OutCounts out_calcos/lbry01i6q_counts_a.fits
    OutFlash  out_calcos/lbry01i6q_lampflash_a.fits
    DETECTOR  FUV, segment A
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G130M, CENWAVE 1291, FPOFFSET 0
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
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    BRFTAB  = lref$s7g1700el_brf.fits
    FLATCORR  OMIT (already complete)
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$w4i1707ol_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Airglow region for OI_1304 is X: 2328 to 2894, Y: 457 to 518
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 27-May-2014 21:56:34 EDT
    Input     out_calcos/lbry01i6q_flt_a.fits
    Incounts  out_calcos/lbry01i6q_counts_a.fits
    Output    out_calcos/lbry01i6q_x1d_a.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$u8k1433ql_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVA spectrum was found at y = 490.32 vs. nominal y = 486.38
        error estimate for y location = 4.22, FWHM = 13.86
    Spectrum will be extracted at y = 486.38
    
    TIME-TAG calibration -- 27-May-2014 21:56:43 EDT
    Input     out_timefilter/lbry01i6q_corrtag_b.fits
    OutTag    out_calcos/lbry01i6q_corrtag_b.fits
    OutFlt    out_calcos/lbry01i6q_flt_b.fits
    OutCounts out_calcos/lbry01i6q_counts_b.fits
    OutFlash  out_calcos/lbry01i6q_lampflash_b.fits
    DETECTOR  FUV, segment B
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G130M, CENWAVE 1291, FPOFFSET 0
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
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    BRFTAB  = lref$s7g1700el_brf.fits
    FLATCORR  OMIT (already complete)
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$w4i1707ol_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Airglow region for LY_ALPHA is X: 9029 to 9214, Y: 517 to 578
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 27-May-2014 21:57:43 EDT
    Input     out_calcos/lbry01i6q_flt_b.fits
    Incounts  out_calcos/lbry01i6q_counts_b.fits
    Output    out_calcos/lbry01i6q_x1d_b.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$u8k1433ql_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVB spectrum was found at y = 548.81 vs. nominal y = 547.27
        error estimate for y location = 6.08, FWHM = 15.17
    Spectrum will be extracted at y = 547.27
    SPWCSTAB= lref$uai1737rl_spwcs.fits
    End   27-May-2014 21:57:57 EDT
    CALCOS version 2.21 (2014-02-18)
    numpy version 1.8.0
    astropy version 0.3
    Begin 27-May-2014 21:57:57 EDT
    Input file = out_timefilter/lbry01icq_corrtag_a.fits
    
    TIME-TAG calibration -- 27-May-2014 21:57:59 EDT
    Input     out_timefilter/lbry01icq_corrtag_a.fits
    OutTag    out_calcos/lbry01icq_corrtag_a.fits
    OutFlt    out_calcos/lbry01icq_flt_a.fits
    OutCounts out_calcos/lbry01icq_counts_a.fits
    OutFlash  out_calcos/lbry01icq_lampflash_a.fits
    DETECTOR  FUV, segment A
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G130M, CENWAVE 1291, FPOFFSET 1
    APERTURE  PSA
    
    HVTAB   = lref$wa41551al_hv.fits
    BADTCORR  OMIT (already complete)
    RANDCORR  OMIT (already complete)
    TEMPCORR  OMIT (already complete)
    GEOCORR   OMIT (already complete)
    WALKCORR  OMIT (already complete)
    DEADCORR  OMIT (already complete)
    Warning:  No GTI table found in raw file.
    Warning:  exposure time in header was 0.000
        exptime has been corrected to 1176.224
    PHACORR   OMIT (already complete)
    DOPPCORR  OMIT (already complete)
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    BRFTAB  = lref$s7g1700el_brf.fits
    FLATCORR  OMIT (already complete)
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    Warning:  No GTI table found in raw file.
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$w4i1707ol_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Airglow region for OI_1304 is X: 2328 to 2894, Y: 457 to 518
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 27-May-2014 21:58:58 EDT
    Input     out_calcos/lbry01icq_flt_a.fits
    Incounts  out_calcos/lbry01icq_counts_a.fits
    Output    out_calcos/lbry01icq_x1d_a.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$u8k1433ql_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVA spectrum was found at y = 531.81 vs. nominal y = 486.38
        error estimate for y location = 999.00, FWHM = -1.00
    Spectrum will be extracted at y = 486.38
    
    TIME-TAG calibration -- 27-May-2014 21:59:06 EDT
    Input     out_timefilter/lbry01icq_corrtag_b.fits
    OutTag    out_calcos/lbry01icq_corrtag_b.fits
    OutFlt    out_calcos/lbry01icq_flt_b.fits
    OutCounts out_calcos/lbry01icq_counts_b.fits
    OutFlash  out_calcos/lbry01icq_lampflash_b.fits
    DETECTOR  FUV, segment B
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G130M, CENWAVE 1291, FPOFFSET 1
    APERTURE  PSA
    
    HVTAB   = lref$wa41551al_hv.fits
    BADTCORR  OMIT (already complete)
    RANDCORR  OMIT (already complete)
    TEMPCORR  OMIT (already complete)
    GEOCORR   OMIT (already complete)
    WALKCORR  OMIT (already complete)
    DEADCORR  OMIT (already complete)
    Warning:  No GTI table found in raw file.
    Warning:  exposure time in header was 0.000
        exptime has been corrected to 1176.224
    PHACORR   OMIT (already complete)
    DOPPCORR  OMIT (already complete)
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    BRFTAB  = lref$s7g1700el_brf.fits
    FLATCORR  OMIT (already complete)
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    Warning:  No GTI table found in raw file.
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$w4i1707ol_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Airglow region for LY_ALPHA is X: 9029 to 9214, Y: 517 to 578
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 27-May-2014 22:00:09 EDT
    Input     out_calcos/lbry01icq_flt_b.fits
    Incounts  out_calcos/lbry01icq_counts_b.fits
    Output    out_calcos/lbry01icq_x1d_b.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$u8k1433ql_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVB spectrum was found at y = 592.33 vs. nominal y = 547.27
        error estimate for y location = 999.00, FWHM = -1.00
    Spectrum will be extracted at y = 547.27
    SPWCSTAB= lref$uai1737rl_spwcs.fits
    End   27-May-2014 22:00:31 EDT
    CALCOS version 2.21 (2014-02-18)
    numpy version 1.8.0
    astropy version 0.3
    Begin 27-May-2014 22:00:31 EDT
    Input file = out_timefilter/lbry01ifq_corrtag_a.fits
    
    TIME-TAG calibration -- 27-May-2014 22:00:34 EDT
    Input     out_timefilter/lbry01ifq_corrtag_a.fits
    OutTag    out_calcos/lbry01ifq_corrtag_a.fits
    OutFlt    out_calcos/lbry01ifq_flt_a.fits
    OutCounts out_calcos/lbry01ifq_counts_a.fits
    OutFlash  out_calcos/lbry01ifq_lampflash_a.fits
    DETECTOR  FUV, segment A
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G160M, CENWAVE 1589, FPOFFSET 0
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
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    BRFTAB  = lref$s7g1700el_brf.fits
    FLATCORR  OMIT (already complete)
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$w4i1707ol_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 27-May-2014 22:01:32 EDT
    Input     out_calcos/lbry01ifq_flt_a.fits
    Incounts  out_calcos/lbry01ifq_counts_a.fits
    Output    out_calcos/lbry01ifq_x1d_a.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$u8k1433ql_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVA spectrum was found at y = 484.94 vs. nominal y = 479.95
        error estimate for y location = 7.06, FWHM = 4.87
    Spectrum will be extracted at y = 479.95
    
    TIME-TAG calibration -- 27-May-2014 22:01:41 EDT
    Input     out_timefilter/lbry01ifq_corrtag_b.fits
    OutTag    out_calcos/lbry01ifq_corrtag_b.fits
    OutFlt    out_calcos/lbry01ifq_flt_b.fits
    OutCounts out_calcos/lbry01ifq_counts_b.fits
    OutFlash  out_calcos/lbry01ifq_lampflash_b.fits
    DETECTOR  FUV, segment B
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G160M, CENWAVE 1589, FPOFFSET 0
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
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    BRFTAB  = lref$s7g1700el_brf.fits
    FLATCORR  OMIT (already complete)
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$w4i1707ol_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 27-May-2014 22:02:46 EDT
    Input     out_calcos/lbry01ifq_flt_b.fits
    Incounts  out_calcos/lbry01ifq_counts_b.fits
    Output    out_calcos/lbry01ifq_x1d_b.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$u8k1433ql_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVB spectrum was found at y = 541.42 vs. nominal y = 539.42
        error estimate for y location = 3.83, FWHM = 5.86
    Spectrum will be extracted at y = 539.42
    SPWCSTAB= lref$uai1737rl_spwcs.fits
    End   27-May-2014 22:03:03 EDT
    CALCOS version 2.21 (2014-02-18)
    numpy version 1.8.0
    astropy version 0.3
    Begin 27-May-2014 22:03:03 EDT
    Input file = out_timefilter/lbry01ilq_corrtag_a.fits
    
    TIME-TAG calibration -- 27-May-2014 22:03:06 EDT
    Input     out_timefilter/lbry01ilq_corrtag_a.fits
    OutTag    out_calcos/lbry01ilq_corrtag_a.fits
    OutFlt    out_calcos/lbry01ilq_flt_a.fits
    OutCounts out_calcos/lbry01ilq_counts_a.fits
    OutFlash  out_calcos/lbry01ilq_lampflash_a.fits
    DETECTOR  FUV, segment A
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G160M, CENWAVE 1589, FPOFFSET 1
    APERTURE  PSA
    
    HVTAB   = lref$wa41551al_hv.fits
    BADTCORR  OMIT (already complete)
    RANDCORR  OMIT (already complete)
    TEMPCORR  OMIT (already complete)
    GEOCORR   OMIT (already complete)
    WALKCORR  OMIT (already complete)
    DEADCORR  OMIT (already complete)
    Warning:  No GTI table found in raw file.
    Warning:  exposure time in header was 0.000
        exptime has been corrected to 1486.208
    PHACORR   OMIT (already complete)
    DOPPCORR  OMIT (already complete)
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    BRFTAB  = lref$s7g1700el_brf.fits
    FLATCORR  OMIT (already complete)
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    Warning:  No GTI table found in raw file.
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$w4i1707ol_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 27-May-2014 22:04:13 EDT
    Input     out_calcos/lbry01ilq_flt_a.fits
    Incounts  out_calcos/lbry01ilq_counts_a.fits
    Output    out_calcos/lbry01ilq_x1d_a.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$u8k1433ql_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVA spectrum was found at y = 525.32 vs. nominal y = 479.95
        error estimate for y location = 999.00, FWHM = -1.00
    Spectrum will be extracted at y = 479.95
    
    TIME-TAG calibration -- 27-May-2014 22:04:22 EDT
    Input     out_timefilter/lbry01ilq_corrtag_b.fits
    OutTag    out_calcos/lbry01ilq_corrtag_b.fits
    OutFlt    out_calcos/lbry01ilq_flt_b.fits
    OutCounts out_calcos/lbry01ilq_counts_b.fits
    OutFlash  out_calcos/lbry01ilq_lampflash_b.fits
    DETECTOR  FUV, segment B
    EXPTYPE   EXTERNAL/SCI
    OPT_ELEM  G160M, CENWAVE 1589, FPOFFSET 1
    APERTURE  PSA
    
    HVTAB   = lref$wa41551al_hv.fits
    BADTCORR  OMIT (already complete)
    RANDCORR  OMIT (already complete)
    TEMPCORR  OMIT (already complete)
    GEOCORR   OMIT (already complete)
    WALKCORR  OMIT (already complete)
    DEADCORR  OMIT (already complete)
    Warning:  No GTI table found in raw file.
    Warning:  exposure time in header was 0.000
        exptime has been corrected to 1486.208
    PHACORR   OMIT (already complete)
    DOPPCORR  OMIT (already complete)
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    BRFTAB  = lref$s7g1700el_brf.fits
    FLATCORR  OMIT (already complete)
    WAVECORR  OMIT (already complete)
    BRSTCORR  OMIT
    Warning:  No GTI table found in raw file.
    DQICORR   PERFORM (complete, but repeat)
    BPIXTAB = lref$w4i1707ol_bpix.fits
    GSAGTAB = lref$w9o1334ml_gsag.fits
    STATFLAG  T
    Update the TIMELINE extension.
    Warning:  spt file not found, so TIMELINE extension is incomplete
    
    Spectral Extraction -- 27-May-2014 22:05:46 EDT
    Input     out_calcos/lbry01ilq_flt_b.fits
    Incounts  out_calcos/lbry01ilq_counts_b.fits
    Output    out_calcos/lbry01ilq_x1d_b.fits
    
    Info:  find-target option = no
    X1DCORR   PERFORM
    XTRACTAB= lref$w5g1439ql_1dx.fits
    DISPTAB = lref$v3g18194l_disp.fits
    HELCORR   PERFORM
    BACKCORR  PERFORM
    STATFLAG  T
    FLUXCORR  PERFORM
    FLUXTAB = lref$u8k1433ql_phot.fits
    TDSCORR   PERFORM
    TDSTAB  = lref$w7h1935dl_tds.fits
    FUVB spectrum was found at y = 584.89 vs. nominal y = 539.42
        error estimate for y location = 999.00, FWHM = -1.00
    Spectrum will be extracted at y = 539.42
    SPWCSTAB= lref$uai1737rl_spwcs.fits
    End   27-May-2014 22:06:09 EDT



    fig = plt.figure(figsize=(17, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    for dataset in glob.glob('out_calcos/lbry*x1d.fits'):
        hdu = fits.open(dataset)
        ax.plot(hdu[1].data['wavelength'].ravel(), hdu[1].data['flux'].ravel(), label=hdu[0].header['rootname'])
        
    ax.set_ylim(0, 2e-13)
    ax.set_xlim(1280, 1340)
    ax.set_ylabel('Flux')
    ax.set_xlabel('Wavelength $(\AA)$')
    ax.legend(shadow=True, numpoints=1, fontsize=16)
    
    fig.savefig('spectrum_after.pdf', bbox_inches='tight')


![png](timefilter_tutorial_files/timefilter_tutorial_11_0.png)



    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    for dataset in glob.glob('out_calcos/lbry*corrtag_a.fits'):
        hdu = fits.open(dataset)
        file_path, file_name = os.path.split(dataset)
        
        ax.plot(hdu['timeline'].data['time'], hdu['timeline'].data['sun_alt'], lw=3)
        ax.annotate(file_name[:9], 
                    xycoords='data',
                    xy=(hdu['timeline'].data['time'][-1], hdu['timeline'].data['sun_alt'][-1]),
                    xytext=(hdu['timeline'].data['time'][-1], hdu['timeline'].data['sun_alt'][-1]),
                    fontsize=16)
        
    ax.axhspan(ymin=-20, ymax=0, xmin=0, xmax=1, color='grey', alpha=.5, label='Orbital Day')
    ax.set_ylim(-20, 120)
    ax.set_xlabel('Time into observation (seconds)')
    ax.set_ylabel('Sun Altitude (degrees)')




    <matplotlib.text.Text at 0x10e254a50>




![png](timefilter_tutorial_files/timefilter_tutorial_12_1.png)



    
