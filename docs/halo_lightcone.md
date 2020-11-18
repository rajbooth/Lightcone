


# The halo lightcone

## Overview
In addition to the galaxy lightcone catalogue described elsewhere in this documentation, a supplementary halo lightcone has been built using the group data generated from the GADGET simulation  dark matter snapshots.

## Lightcone build
The process of building the halo lightcone is broadly similar to that used for the galaxy lightcone, in that the full lightcone is constructed by joining together eight octants from the simulation box to construct a full sky volume. 

## Lightcone data
### Data location
The entire halo lightcone is stored in a single $7 Gbyte$ hdf5 data file, which can be found at:
[Halo Lightcone](Canonical%20Path%09/cosma6/data/dp004/dc-boot5/Lightcone/Halo_FullSky)

### File structure
The lightcone file comprises 22 datasets, one for each of the 22 GADGET snapshots that were used in its construction, ranging from $63 (z=0)$ to $41 (z\simeq 1.4)$. The dataset key names are of the form *snapshot_xx*  where *xx* is the snapshot number.

### Data fields
The lightcone file provides positional data for each halo in both celestial and Cartesian coordinates.  The following data fields are included in the dataset:
| Parameter |         Description     | Data type| Units |
|----------|-----------------------------|-----------|-------|
| cm| position of halo centre of mass (x,y,z) | vect | Mpc/h |
| vel| velocity of halo centre of mass (x,y,z) | vect |km/s |
| pos| position of central halo particle (x,y,z) | vect | Mpc/h
| r| co-moving distance of central halo particle  |(float32) |Mpc/h |
| ra | right ascention of central halo particle  |(float32) | $(\pm 90^\circ)$ 
| dec| co-moving distance of central halo particle  |(float32) |$(\pm 90^\circ)$ 
| zz | redshift of halo  |(float32) |Mpc/h |
| vr | radial line-of-sight velocity of halo  |(float32) |km/s |
| vel_disp | velocity dispersion  |(float32) |km/s |

## Visualisation

## Reading the lightcone file
The following code snippet illustrates how to read from the halo lightcone file and view a slice though the lightcone.
```python
# define data structures
vect = np.dtype([('x', np.float32), ('y', np.float32), ('z', np.float32)])
halo2 = np.dtype([('cm', vect), ('vel', vect), ('pos', vect),
                  ('mass', np.float32), ('r', np.float32), ('ra', np.float32), ('dec', np.float32), 
                  ('zz', np.float32), ('vr', np.float32), ('vel_disp', np.float32)])
```

![enter image description here](https://raw.githubusercontent.com/rajbooth/Lightcone/master/images/Halo_lightcone.png)

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTIwODU5NjExMjQsLTc3NjA5NTQyMiwxOT
I4NjUzMjMyXX0=
-->