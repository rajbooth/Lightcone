


# The halo lightcone

## Overview
In addition to the galaxy lightcone catalogue described elsewhere in this documentation, a supplementary halo lightcone has been built using the group data generated from the GADGET simulation  dark matter snapshots.

## Lightcone build
The process of building the halo lightcone is broadly similar to that used for the galaxy lightcone, in that the full lightcone is constructed by joining together eight octants from the simulation box to construct a full sky volume. 

## Lightcone data
### Data location
The entire halo lightcone is stored in a single $7 Gbyte$ hdf5 data file, which can be found at:
[/cosma6/data/dp004/dc-boot5/Lightcone/Halo_FullSky](/cosma6/data/dp004/dc-boot5/Lightcone/Halo_FullSky)

### File structure
The lightcone file comprises 22 datasets, one for each of the 22 GADGET snapshots that were used in its construction, ranging from $63 (z=0)$ to $41 (z\simeq 1.4)$. The dataset key names are of the form *snapshot_xx*  where *xx* is the snapshot number.

### Data fields
The lightcone file provides positional data for each halo in both celestial and Cartesian coordinates.  The following data fields are included in the dataset:

| Parameter |         Description     | Data type| Units |
|----------|-----------------------------|-----------|-------|
| cm | position of halo centre of mass (x,y,z) | vect | Mpc/h |
| vel | velocity of halo centre of mass (x,y,z) | vect |km/s |
| pos | position of central halo particle (x,y,z) | vect | Mpc/h |
| m | halo mass  |float32| $10^{10} M_{\odot}$ |
| r | co-moving distance of central halo particle  |float32 |Mpc/h |
| ra | right ascention of central halo particle  |float32 | $(\pm 90^\circ)$ |
| dec | declination of central halo particle  |float32 |$(\pm 90^\circ)$ |
| zz | redshift of halo  |float32 |Mpc/h |
| vr | radial line-of-sight velocity of halo  |float32)|km/s |
| vel_disp | velocity dispersion  |float32 |km/s |

## Visualisation

## Reading the lightcone file
The following code snippet illustrates how to read from the halo lightcone file and view a slice though the lightcone.
```python
# define data structures
vect = np.dtype([('x', np.float32), ('y', np.float32), ('z', np.float32)])
halo2 = np.dtype([('cm', vect), ('vel', vect), ('pos', vect),
                  ('mass', np.float32), ('r', np.float32), ('ra', np.float32), ('dec', np.float32), 
                  ('zz', np.float32), ('vr', np.float32), ('vel_disp', np.float32)])

# specify file location
file = '/cosma6/data/dp004/dc-boot5/Lightcone/Halo_FullSky/halo_lightcone'
# initialise halo data array
haloes = np.empty(0,dtype=halo2)
with h5py.File(file, 'r') as f:
    # read each snapshot dataset into halo array
    for snap in f.keys():
        haloes = np.append(haloes,f[snap])
        print(snap, haloes.shape)

depth = 0.2 # set thickness of data slice
r_max = 3000 # cut off radial distance
# filter slice haloes from data set based on depth and r_max
h = haloes[abs(haloes['pos']['z'])<depth ]  
h = h[h['r']<r_max]

# populate halo parameter arrays
x = h['pos']['x']
y = h['pos']['y']
m = h['mass']
vd = h['vel_disp']

# plot halo slice as a scatter graph
area = m/100  # make area of blob proportional to halo mass
col = np.log(vd) # make colour proportional to velocity dispersion
fig = plt.subplots(figsize = (12,12))
plt.scatter(x, y, s=area, c=col, ec = None, alpha=0.8)
#plt.savefig('Halo lightcone.png')
plt.show()
```
![Halo Lightcone](https://raw.githubusercontent.com/rajbooth/Lightcone/master/images/Halo_lightcone.png)
The fact that the halo lightcone also contains positional data in celestial coordinates means that it is easy to plot whole-sky data for a radial shell through the lightcone, as in this code snippet.

```python
import healpy as hp
import matplotlib.pyplot as plt
# specify file location
file = '/cosma6/data/dp004/dc-boot5/Lightcone/Halo_FullSky/halo_lightcone'
snap = 44 # specify sdnapshot to be viewed
nside = 256 # set pixel resolution
npix = hp.nside2npix(nside)
ds_name = 'snapshot_{0:02d}'.format(snap)
with h5py.File(file, 'r') as f:
    # read snapshot shell into pixel map
    shell = f[ds_name]
    pix = hp.ang2pix(nside, shell['ra'], shell['dec'], lonlat = True)
# convert pixel map into a density map
dens_map = np.histogram(pix, npix)[0]

# use Healpy to plot halo shell
fig, ax = plt.subplots(1,1, figsize=(12,12))
hp.mollview(dens_map, norm = 'hist', title = 'Snapshot {0:02d}'.format(snap), hold=True, xsize = 3600)
```
![Halo Lightcone Shell](https://raw.githubusercontent.com/rajbooth/Lightcone/master/images/Halo_lightcone_shell.png)
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTg4Mzg4OTcyOSwtMTM3NTA5Mzc4LC0xNz
E0MTAwNzAyLC0xODY0NDUxNjgsLTE0NTA5ODczODYsMjEyOTA2
NTUxNSwyMDMzMjU0NjA2LDE1MTQyOTM2NywtNzc2MDk1NDIyLD
E5Mjg2NTMyMzJdfQ==
-->