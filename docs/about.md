## Overview

## Building the lightcone

### Source data
Gadget snapshots.

### The luminosity function

$$M = m - 5 (\log_{10}D_L - 1) $$

where $M$ is absolute magnitude and $m$ is apparent magnitude

 $$m = M + 5 (\log_{10}D_L - 1)$$ 

$$\rightarrow D_L = 10^{0.2(m-M) - 1}$$ 

where \(D_L = (1+z) D_C \)

Schechter:
\( M^*_J = -21.4 \)

$$M_{bol} = M_\odot -2.5 \log_{10} \frac{L_*}{L_\odot}$$

$$ \rightarrow \frac{L_*}{L_\odot} = 10^{0.4(M_{bol,\odot} - M_{bol,*})} $$

Efstathieu (Les Houches Lectures):

"In a survey limited by apparent magnitude (i.e. a flux limited samplegalaxy of luminosity $L$ can be seen out to a distance $d_{max}$ (redshift $z_{max}$) given by:"

$$5 \log d_{max}(L) = m_{lim} - M_\odot - 25 + 2.5 \log(L/L_\odot) + kz_{max}  $$

$$ \rightarrow \frac{L}{L_\odot} = d_{max}^2 10^{0.4(M_\odot - m_{lim} + 25)} $$


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
where x is 
### Data fields

### Co-ordinate system

## Reading lightcone data
``` python
def add_snapshot(snap):
    ngals = 0
    fname = 'galaxy_lightcone.snap{0:02d}'.format(snap)
    with h5py.File(fpath+fname,'r') as fi:
        for k in fi.keys():
            print (k, fi[k].shape) 
            if True: # k=='octant_0':
                g = fi[k]
                gals = g[(g['Dec']>85) & (g['Dec']<95)]
                L.extend(gals['L'])
                r.extend(gals['r'])
                ra.extend(gals['RA'])
                dec.extend(gals['Dec'])
                RSD.extend(gals['RSD'])
                zz.extend(gals['z'])
                ngals += len(gals)
    return ngals
```
## Visualising the lightcone
![Slice through galaxy lightcone](https://github.com/rajbooth/Lightcone/raw/master/images/particle_lightcone_Particle_z_particle_mass.png)
The file explorer is accessible using the button in left corner of the navigation bar. You can create a new file by clicking the **New file** button in the file explorer. You can also create folders by clicking the **New folder** button.


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTIxNDUwOTQ5NTIsMTczMDk2NDA2LC03Nj
EzMDcyNzYsNzUzMzc1Njc3LDIwODMwNTk2MTIsLTE2MzU2NDEw
NzUsMTgwNjMxNzUzLDg3MTk4NTU2Myw5MTgxOTk0NTEsMTQ4MD
gzMzQsLTkzNzk4ODYxOCw2MDA1NDE4NzgsLTE4NjE4OTQwODZd
fQ==
-->