## Overview
The galaxy lightcone catalogue provides a dataset of over 2 billion galaxies extending out to a redshift of $z \simeq 0.85$, contained in a co-moving cubic volume of $6000 Mpc^3/h^3$.  In addition to the positional co-ordinates of each galaxy, the lightcone dataset also assigns a luminosity value to each galaxy such that the overall luminosity distribution fits a [Schechter function](#luminosity).  This enables the user to apply a magnitude filter to the galaxy dataset to extract only those galaxies that will be visible in an observational galaxy survey for a given instrument sensitivity.
The dataset also holds the redshift space distortion factor associated with each galaxy, derived from its peculiar velocity 
## Building the lightcone

### Source data
Gadget snapshots.
```
```
### <a name="luminosity"></a>The luminosity function
![Minimum luminosity](https://github.com/rajbooth/Lightcone/raw/master/images/Min_Lum_Redshift.png)
$$M = m - 5 (\log_{10}D_L - 1) $$

where $M$ is absolute magnitude and $m$ is apparent magnitude

 $$m = M + 5 (\log_{10}D_L - 1)$$ 

$$\rightarrow D_L = 10^{0.2(m-M) - 1}$$ 

where \(D_L = (1+z) D_C \)

Schechter:
\( M^*_J = -21.4 \)

$$M_{bol} = M_\odot -2.5 \log_{10} \frac{L_*}{L_\odot}$$

$$ \rightarrow \frac{L_*}{L_\odot} = 10^{0.4(M_{bol,\odot} - M_{bol,*})} $$
![enter image description here](https://github.com/rajbooth/Lightcone/raw/master/images/Cumulative_Probability_Distribution.png)
Efstathieu (Les Houches Lectures):

"In a survey limited by apparent magnitude (i.e. a flux limited samplegalaxy of luminosity $L$ can be seen out to a distance $d_{max}$ (redshift $z_{max}$) given by:"

$$5 \log d_{max}(L) = m_{lim} - M_\odot - 25 + 2.5 \log(L/L_\odot) + kz_{max}  $$

$$ \rightarrow \frac{L}{L_\odot} = d_{max}^2 10^{0.4(M_\odot - m_{lim} + 25)} $$
![Luminosity distribution](https://github.com/rajbooth/Lightcone/raw/master/images/Luminosity_Distribution.png)

### Redshift space distortion

## Data products
The primary galaxy lightcone catalogue is located on COSMA server at Durham, at:
[Galaxy Full Sky Lightcone](/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky)
This is just over 1 Tbyte in size and hence is too large to use other than by applications that are also running on the COSMA server.

A reduced size lightcone has been generated from this primary dataset, which can be used for testing purposes and at 2.8 Gbytes, is small enough for download to personal off-site computers.  This can be found at:
[Luminosity and redshift limited lightcone](/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky_Reduced/galaxy_lightcone_M_limited.h5)
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
* Dec - declination in degrees $(\pm 90 \degree)$
* RA - right ascension  $(\pm 180 \degree)$
* z - redshift
* L - luminosity ratio $(L/L^*)$
* RSD - redshift distortion $(Mpc)$

Data stored as a *gal* datatype, defined as:
```python
#define galaxy datatype
gal = np.dtype([('r', np.float32),('RA', np.float32),('Dec', np.float32),('z', np.float32),('RSD', np.float32),('L', np.float32)])
```
### Co-ordinate system

## Reading lightcone data
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
*[to follow - transforming to Cartesian co-ordinates and visualisation in yt]*

<!--stackedit_data:
eyJoaXN0b3J5IjpbMTI3MjIzOTMwMSwtMTQxNTExODU5OSwxOD
ExNDI5Nzk1LC0xNDczNTM5Mzg3LC0xMDUzMjcyMDI4LDE0MTg5
NzY0MDEsNDY1NDU3NzcyLDE3MzA5NjQwNiwtNzYxMzA3Mjc2LD
c1MzM3NTY3NywyMDgzMDU5NjEyLC0xNjM1NjQxMDc1LDE4MDYz
MTc1Myw4NzE5ODU1NjMsOTE4MTk5NDUxLDE0ODA4MzM0LC05Mz
c5ODg2MTgsNjAwNTQxODc4LC0xODYxODk0MDg2XX0=
-->