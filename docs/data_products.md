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
The primary lightcone dataset consists of 21 files, one for each Gadget snapshot, ranging from 63 ($z=0$)  to 42 ($z\simeq1$).  Each file is divided into 8 hdf5 datasets, one for each full-sky octant.  The datasets are named 'octant_0' though to 'octant_7'.  The largest of these files for snapshot 42 is over 180 Gbytes in size and cannot therefore be 

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
eyJoaXN0b3J5IjpbNjA4NzgwMzgwLDE3MjUwNzcwNzgsLTIwOD
Y0MzQxOTgsLTE3MzMxNzE3ODAsODE2NzQzNTEyLC0xMTU2OTYw
ODQyLDEwMDQ4MDYzNDIsNDA1MDM3NzgyLC0xMDQzMzQ4MDgwLC
0yMTM0NDQ2ODU0LDEzMDcwMzU4NSwxNTA4NzMyMTYyLDkwMDYy
MzE4MywtNjQ2MjE5NjA3LC0xNTQyMjg1OTAyLC0xNDg2ODc5Mz
kyLDE2NzE2MDQ4OCwtMTQzNTc1NjY2MSwxOTAwMjU1MDgwLDEz
ODc0MzMyODFdfQ==
-->