## Data location
### Primary lightcone catalogue
The primary galaxy lightcone catalogue is located on COSMA server at Durham, at:
[/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky](/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky)
This is just over 1 Tbyte in size and hence is too large to use other than by applications that are also running on the COSMA server.

### Compact lightcone catalogue
A reduced size lightcone has been generated from this primary dataset, which can be used for testing purposes and at 2.8 Gbytes, is small enough for download to personal off-site computers.  This can be found at:
[/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky_Reduced/galaxy_lightcone_M_limited.h5](/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky_Reduced/galaxy_lightcone_M_limited.h5)
This lightcone only extends out to a redshift of $z<=0.4$ and a luminosity cut-off has been applied such that only galaxies with a luminosity $L<= 21$ are included in the catalogue.

### Matter power spectra
Although not strictly relevant to the lightcone itself, this document also describes some ancillary data that was generated as part of the process of building the lightcone: specifically, a set of matter power spectrum measurements for each of the snapshots that were include din the galaxy lightcone.  These are described in more detail in the section: [matter power spectrum](#power_spectrum)

### File structure
The lightcone data is stored in files conforming to the hdf5 file format.
#### Primary dataset
Data for each shell derived from a specific Gadget snapshot is stored in a separate file, with a name of the form:
```
galaxy_lightcone.snapxx
```
where xx is the snapshot number ranging from 42 - 63.
Each of these files contains 8 datasets, corresponding to each of the 8 octants that comprise the complete spherical volume of the lightcone, with a dataset name of the form:
```
octant_x
```
where x is the octant number in the range 0 - 7
#### Reduced dataset
All the data for the entire reduced dataset is stored in a single dataset ('galaxies') , contained within a single file.

### Data fields
* r - co-moving radial distance $(Mpc)$
* Dec - declination in degrees $(\pm 90^\circ)$
* RA - right ascension  $(\pm 180^\circ)$
* z - redshift
* L - luminosity ratio $(L/L^*)$
* RSD - redshift distortion $(Mpc)$

Data stored as a *gal* datatype, defined as:
```python
#define galaxy datatype
gal = np.dtype([('r', np.float32),('RA', np.float32),('Dec', np.float32),('z', np.float32),('RSD', np.float32),('L', np.float32)])
```
### Co-ordinate system
Positional data is stored in celestial co-ordinates, i.e. radial co-moving distance to observer $r\;(Mpc/h)$, declination from equatorial plane $Dec\; (\pm 90^\circ)$ , and right ascension $RA\;(\pm 180^\circ)$.
## <a name="reading_lightcone"></a>Reading lightcone data
### Prerequisites
Since the lightcone is stored in hdf5 format, it is necessary to load a compatible hdf5 library in order to access the lightcone dataset.  Additionally, other utilities may prove useful in analysing this data.  For Python users, the following imports are recommended:
```python
import h5py
import healpy as hp
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.interpolate import UnivariateSpline as spl
from astropy.cosmology import Planck15 as cosmo, z_at_value
```
### Primary dataset
The primary lightcone dataset consists of 21 files, one for each Gadget snapshot, ranging from 63 ($z=0$)  to 42 ($z\simeq1$).  Each file is divided into 8 hdf5 datasets, one for each full-sky octant.  The datasets are named 'octant_0' though to 'octant_7'.  The largest of these files for snapshot 42 is over 180 Gbytes in size and cannot therefore be read in one chunk, whereas an individual octant at about 22.5 Gbytes is just about manageable. 
The following python example illustrates how the full dataset can be read and analysed, in this case for the purposes of constructing raytracing lens planes.
```python
# Build surface mass density shells from lightcone
import numpy as np
import h5py
import math
import healpy as hp
from datetime import datetime
import os.path
from fast_histogram import histogram1d, histogram2d

# read particle data from lightcone 
def read_octant(snap, oct):
    fpath = '/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky/'
    lens_path = '/cosma6/data/dp004/dc-boot5/lensplanes/OnionSkin/'
    fname = 'galaxy_lightcone.snap{0:02d}'.format(snap)
    key = 'octant_{0:0d}'.format(oct)
    print('Started, snapshot =', snap,  datetime.now(), flush = True)
    with h5py.File(fpath+fname,'r') as fi:
        gals = fi[key]
        print('Opened dataset, octant =', oct, datetime.now(), flush = True)
        r = gals['r']/h
        phi = gals['RA'] * np.pi / 180
        theta = (gals['Dec'] + 90) * np.pi / 180 
        print('Read variables, number of particles = ',len(r), datetime.now(), flush = True)
           
    lp = r//delta # calculate lens plane number from radial comoving distance
    lp = lp.astype(np.int, copy = False)
    lps = np.unique(lp)   
    global first_lp 
    first_lp = lps.min()
    last_lp = lps.max()
    nlp = len(lps)
    print('Calculated lensplanes:', lps, datetime.now(), flush = True)
    
    pix = hp.ang2pix(nside,theta,phi)
    print('Determined pixel', datetime.now(), flush = True)
    
    sigma = np.zeros((nlp, npix))  # seem to need this to reserve memory and prevent segmentation fault
    sigma = histogram2d(lp, pix, range = [[first_lp, last_lp+1], [0, npix+1]],  bins=[nlp,npix])
    print('Sigma histogram completed', datetime.now(), flush = True)
    
    return sigma

## set resolution
O_DE = 0.6914
o_b  = 0.022161
o_c = 0.11889
h = np.sqrt((o_b + o_c)/(1-O_DE))
O_m = (o_b + o_c)/(h*h)
delta = 200 # lens plane thickness in Mpc
nside = 8192 
npix = hp.nside2npix(nside)
dens_fac =  2.69e11 * const.M_sun /h * npix  / (4 * np.pi)
print('Nside =', nside, ' npix =',npix, 'h = ',h,'Omega_m = ', O_m)

lens_path = '/cosma6/data/dp004/dc-boot5/lensplanes/OnionSkin/'
for snap in range(43, 42, -1):
    for oct in range(8):
        if oct==0:
            sigma = read_octant(snap,oct)
        else:
            sigma += read_octant(snap,oct)
        
    nlp = sigma.shape[0]
    for l in range(nlp):     
        # write sigma to file
        plane = l + first_lp
        sigma_file = lens_path + 'Lens_plane_{0:0d}_nside_{1:0d}'.format(plane, nside)
        # check whether file already exists
        if os.path.isfile(sigma_file):
            sigma0 = hp.read_map(sigma_file, dtype = np.int)
            # add new sigma to old sigma
            sigma[l] += sigma0
            
        hp.write_map(sigma_file, sigma[l], dtype = np.float32, overwrite = True)
        print ('sigma written to lens plane {0:0d}'.format(plane), flush = True)
```
### Reduced dataset

The following code snippet can be used to read lightcone data from the reduced dataset file into a series of numpy arrays for subsequent processing and analysis.
``` python
# read data from test file
outpath = '/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky_Reduced/'
fname = outpath + 'galaxy_lightcone_M_limited.h5'
with h5py.File(fname,'r') as fi:
    # open the galaxies dataset
    gals = fi['galaxies']
    # extract each data field into separate arrays
    L = gals['L']
    zz = gals['z']
    r = gals['r']
    ra = gals['RA']
    dec = gals['Dec']
    RSD = gals['RSD']
print('Finished reading {0:01d} galaxies'.format(len(r)))
```


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE4MDkyMzc2NTYsMTcyNTA3NzA3OCwtMj
A4NjQzNDE5OCwtMTczMzE3MTc4MCw4MTY3NDM1MTIsLTExNTY5
NjA4NDIsMTAwNDgwNjM0Miw0MDUwMzc3ODIsLTEwNDMzNDgwOD
AsLTIxMzQ0NDY4NTQsMTMwNzAzNTg1LDE1MDg3MzIxNjIsOTAw
NjIzMTgzLC02NDYyMTk2MDcsLTE1NDIyODU5MDIsLTE0ODY4Nz
kzOTIsMTY3MTYwNDg4LC0xNDM1NzU2NjYxLDE5MDAyNTUwODAs
MTM4NzQzMzI4MV19
-->