================================================================================
Title: Updated Photon Scattering Coefficients (g_values) for Mercury's 
    Exospheric Species
Authors: Killen R.M, Vervack Jr R.J., Burger M.H.
================================================================================
Description of contents: This tar.gz archive contains the data behind Figure
    14 of the accepted manuscript listed above. Specifically, this archive
    contains the updated g_values for the three FUV OI lines at 1302.17,
    1304.86, and 1306.03 A, which is provided in these data file(s):
    
    OI_1302_gvalues.txt
    OI_1304_gvalues.txt
    OI_1306_gvalues.txt
    
System requirements: The files are plain text ASCII files with a four line 
    header followed by four columns of data. The four columns of data are
        
    1. Heliocentric relative velocity, ranging from -50 to +50 km/s 
        in 0.1 km/s increments
    2. g_value in photons/atom/s for each line of the triplet.
        
    The typical header looks like:
    
                      		 g_value (photons/atom/s)
    Velocity (km/s)         1302.17 A 1304.86 A  1306.03 A      <= rest frame wavelength of line


Additional comments: Figure 14b is the spectral irradiance file: wavelength
    (A), Irrdiance (Watts/m2/Angstrom). For the oxygen triplet lines between
    1302-1306 A, we generated line shape profiles according to Gladstone
    (1992) and then scaled them to the appropriate irradiance levels in the
    same SORCE SOLSTICE spectrum used for the sulfur triplet lines.

================================================================================
