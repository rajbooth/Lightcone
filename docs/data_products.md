## Data location
The primary galaxy lightcone catalogue is located on COSMA server at Durham, at:
[/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky](/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky)
This is just over 1 Tbyte in size and hence is too large to use other than by applications that are also running on the COSMA server.

A reduced size lightcone has been generated from this primary dataset, which can be used for testing purposes and at 2.8 Gbytes, is small enough for download to personal off-site computers.  This can be found at:
[/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky_Reduced/galaxy_lightcone_M_limited.h5](/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky_Reduced/galaxy_lightcone_M_limited.h5)
This lightcone only extends out to a redshift of $z<=0.4$ and a luminosity cut-off has been applied such that only galaxies with a luminosity $L<= 21$ are included in the catalogue.

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
If it intended to visualise the lightcone data then the use of the *yt* visualisation package is recommended for Python users.
```python
import yt
from yt.units import parsec, Msun
from yt.visualization.volume_rendering.api import Scene, VolumeSource, Camera
```
### Primary dataset
*[to follow]*
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
## Visualising the lightcone
![Luminosity limited galaxy lightcone](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-19.png)
The fact that the lightcone dataset is stored as a hdf5 file makes the process of extracting a sub-set of the data for visualisation purposes relatively simple.  Using the h5py Python library makes it possible to use the standard numpy array indexing syntax to define a specific data slice to extract from the entire dataset.  Extending the python code snippet in [reading lightcone data](#reading_lightcone), we add a slice filter to include only those galaxies that lie within a specified absolute declination:
```python
# read data from test file
depth = 5 # depth of slice in degrees
fname = outpath + 'galaxy_lightcone_M_limited.h5'
with h5py.File(fname,'r') as fi:
    # open the galaxies dataset
    g = fi['galaxies']
    # construct a filter to inlude only those galaxies with a luminosity greater than the redshift dependent cut-off
    f0 = tuple([(abs(g['Dec'])<depth)]) #note: in Python 3, an array filter must be a tuple rather than another array
    gals = g[f0]
    # extract each data field into separate arrays
    L = gals['L']
    zz = gals['z']
    r = gals['r']
    ra = gals['RA']
    dec = gals['Dec']
    RSD = gals['RSD']
print('Finished reading {0:01d} galaxies'.format(len(r)))
```
In order to create a 2D image of our lightcone slice, we first need to convert from celestial to Cartesian co-ordinates.  
```python
# convert from polar to cartesian coords
theta = np.radians(dec)
phi = np.radians(ra)
x = r * np.cos(theta) * np.sin(phi)
y = r * np.cos(theta) * np.cos(phi)
z = r * np.sin(theta)
```
We now need to calculate the intensity value to display in our lightcone image.  There are various ways this could be achieved.  One approach would simply be to do a count of the total number of galaxies that fall within an image pixel.  A more sophisticated approach that we employ here is to use the assigned galaxy luminosities as a proxy for galaxy mass.
The code snippet below offers two alternative approaches to assigning luminosity to an image pixel:
- Nearest Grid Point (NGP)
- Cloud in Cell (CIC)
```python
# set density deposition mode
mode = 'NGP' # alternative is 'CIC'
size = 2000 # size if image array/2
X = x * size/x_max + size-1 # rescale x axis
Y = y * size/y_max + size-1 # rescale y axis

# Calculate cell in which each galaxy is located on 2D grid
i = np.floor(X)
j = np.floor(Y)
i = i.astype(int) # x index of cell that contains particle
j = j.astype(int) # y index of cell that contains particle

# define image mesh
# fill it with low non-zero value as we will want to use a log scale
img = np.full((2*size,2*size),0.0001)

# calcuate total luminosity in each image grid cell by summing over all galaxies
for k,l in enumerate(L):
    if mode=='CIC':
        d_x = X[k] - i[k]  # x offset of particle in cell
        d_y = Y[k] - j[k]  # y offset of particle in cell
        t_x = 1 - d_x
        t_y = 1 - d_y

        img[i[k],j[k]]    += l * t_x * t_y
        if i[k]<2*size-1: img[i[k]+1,j[k]]  += l * d_x * t_y
        if j[k]<2*size-1: img[i[k],j[k]+1]  += l * t_x * d_y
        if i[k]<2*size-1 and j[k]<2*size-1: img[i[k]+1,j[k]+1]+= l * d_x * d_y
    else:       # if NGP mode
        img[i[k],j[k]]+= l
        
L_max = np.amax(img)
L_min = np.amin(img)
print ('L_max = ', L_max, 'L_min = ', L_min)
```
Finally, we can display the image array created in the previous step by invoking the 

> MatPlotLib imshow

 method, as illustrated here:
```python
Dpi = 600  					# set image resolution 
Cmap = cm.get_cmap('hot')	# colour map to use
Interp = 'None'				# anti-aliasing type, if required
Norm = cm.colors.LogNorm(vmin=0.01, vmax=L_max) # normalise luminosity values on log scale
ext = (-r_max, r_max, -r_max, r_max)
fig = plt.subplots(figsize = (12,12))
plt.imshow(img,  cmap = Cmap, norm = Norm, interpolation = Interp, extent = ext)
cbar = plt.colorbar(shrink = 0.8)
cbar.set_label('Luminosity')
plt.xlabel('x (Mpc)')
plt.ylabel('y (Mpc)')
plt.savefig('FullSky_Galaxy_Slice_M0-{5:0.0f}_mode-{0:s}_interp-{1:s}_res-{2:0d}_dpi-{3:0d}_cmap-{4:s}'.format(mode,Interp, size*2, Dpi, cmap_name, M0), dpi = Dpi)
plt.show()
```
![Lightcone](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-19_mode-CIC_interp-kaiser_res-4000_dpi-600_blue.png)
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTUwODczMjE2Miw5MDA2MjMxODMsLTY0Nj
IxOTYwNywtMTU0MjI4NTkwMiwtMTQ4Njg3OTM5MiwxNjcxNjA0
ODgsLTE0MzU3NTY2NjEsMTkwMDI1NTA4MCwxMzg3NDMzMjgxLC
0zNjM5NjUyODEsNTIzMzAwMTYsOTA4MTE5NDIwLC0xODU2Njc2
MDMsLTE4NTY2NzYwMywtNTQ4ODA2NDk2XX0=
-->