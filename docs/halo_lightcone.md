


# The halo lightcone

## Overview
In addition to the galaxy lightcone catalogue described elsewhere in this documentation, a supplementary halo lightcone has been built using the group data generated from the GADGET simulation  dark matter snapshots.

## Lightcone build
The process of building the halo lightcone is broadly similar to that used for the galaxy lightcone, in that the full lightcone is constructed by joining together eight octants from the simulation box to construct a full sky volume. 

## Lightcone data
### Data location
The entire halo lightcone is stored in a single hdf5 data file, which can be found at:
[Halo Lightcone](Canonical%20Path%09/cosma6/data/dp004/dc-boot5/Lightcone/Halo_FullSky)

### Data structure
The lightcone file comprises 22 datasets, one for each of the 42 GADGET snapshots that were used in its construction, ranging from $63 (z=0)$ to $42 (z\simeq 1.4)$.

## Visualisation
![enter image description here](https://raw.githubusercontent.com/rajbooth/Lightcone/master/images/Halo_lightcone.png)

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTEwMDE5MjczNjUsMTkyODY1MzIzMl19
-->