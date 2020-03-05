
## Visualising the lightcone
![Luminosity limited galaxy lightcone](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-19.png)
### Generating a slice plot
The fact that the lightcone dataset is stored as a hdf5 file makes the process of extracting a sub-set of the data for visualisation purposes relatively simple.  Using the h5py Python library makes it possible to use the standard numpy array indexing syntax to define a specific data slice to extract from the entire dataset.  Extending the python code snippet in [reading lightcone data](#reading_lightcone), we add a slice filter to include only those galaxies that lie within a specified absolute declination:
```python
# read data from test file
depth = 510 # depth of slice in degreesMpc
fname = outpath + 'galaxy_lightcone_M_limited.h5'
with h5py.File(fname,'r') as fi:
    # open the galaxies dataset
    g = fi['galaxies']
    # construct a filter to inlude only those galaxies with a luminosity greater than the redshift dependent cut-off
    fz0 = tuple([(abs(g['Dec'])g['r'] * np.sin(np.radians(abs(g['Dec']))) # calculate z co-ordinate of particle
    f0 = tuple([(z0<depth)]) #note: in Python 3, an array filter must be a tuple rather than another array
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
      if mode=='CIC' and i[k]<2*size-1 and j[k]<2*size-1:
        d_x = X[k] - i[k]  # x offset of particle in cell
        d_y = Y[k] - j[k]  # y offset of particle in cell
        t_x = 1 - d_x
        t_y = 1 - d_y

        img[i[k],j[k]]    += l * t_x * t_y
        if i[k]<2*size-1: img[i[k]+1,j[k]]  += l * d_x * t_y
        if j[k]<2*size-1: img[i[k],j[k]+1]  += l * t_x * d_y
        if i[k]<2*size-1 and j[k]<2*size-1: img[i[k]+1,j[k]+1]+= l * d_x * d_y
    else:       # if NGP mode or on grid boundary
        img[i[k],j[k]]+= l
        
L_max = np.amax(img)
L_min = np.amin(img)
print ('L_max = ', L_max, 'L_min = ', L_min)
```
Finally, we can display the image array created in the previous step by invoking the *MatPlotLib imshow* method, as illustrated here:
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

![Lightcone](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-18_mode-CIC_interp-None_res-4000_dpi-600_cmap-hot.png)

### Applying a luminosity filter
We can now make use of the luminosity data that is associated with each galaxy in order to filter only only those galaxies that will be visible in a survey using an instrument of a given magnitude sensitivity, $M_0$ .
This can be achieved simply by adding a luminosity filter to the filter expression used when reading in data from the lightcone dataset, as shown here:
```python
    f0 = tuple([(z0<depth) & (g['L'] > z2L(g['z']))])
```
![Luminosity filtered](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-18_mode-CIC_interp-kaiser_res-4000_dpi-600_cmap-blue2_new.png)

### Plotting galaxy density
Use galaxy counts rather than luminosity value.
![NGP plot, no interpolation](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-18_mode-NGP_interp-None_res-4000_dpi-600_cmap-bone.png)
NGP plot, no interpolation

![CIC plot, no interpolation](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-18_mode-CIC_interp-None_res-4000_dpi-600_cmap-bone.png)
CIC plot, no interpolation

![NGP plot, gaussian interpolation](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-18_mode-NGP_interp-gaussian_res-4000_dpi-600_cmap-bone.png)
NGP plot, gaussian interpolation

![NGP plot, kaiser interpolation](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-18_mode-NGP_interp-kaiser_res-4000_dpi-600_cmap-bone.png)
NGP plot, kaiser interpolation

### Redshift Space Distortion

![Redshift Space Distortion](https://github.com/rajbooth/Lightcone/raw/master/images/FullSky_Galaxy_Slice_M0-19_mode-CIC_interp-None_res-4000_dpi-600_cmap-blue2.png)

## Sky maps
The full lightcone dataset readily lends itself to the generation of full sky density or luminosity maps.  The simplest approach is to read an entire snapshot file, which is inherently comprised of all galaxies within a given redshift range. Using the  *ang2pix* function in *healpy* Python implementation of Healpix, we can assign each galaxy to a pixel on a sky map, then sum the galaxy luminosities in each pixel to generate the final map.  The following code snippet illustrates this process.
```python
## set resolution
nside = 1024
npix = hp.nside2npix(nside);
print('Nside =', nside, ' npix =',npix)
# initialise map arrays
m = np.full((20,npix),hp.UNSEEN)

snaps = [62,59,56,53]

for snap in snaps:
    # Initialise data arrays
    L = []
    zz = []
    r = []
    ra = []
    dec = []
    RSD = []
    
    add_snapshot(snap)

    ra = np.asarray(ra)
    dec = np.asarray(dec)
    
    phi = ra * np.pi / 180
    theta = (dec + 90) * np.pi / 180 

    pix = hp.ang2pix(nside,theta,phi)
    for i, px in enumerate(pix):
        sn = snap-44
        if (m[sn][px] == hp.UNSEEN):
            m[sn][px] = L[i]
        else:
            m[sn][px] += L[i]
    print('Done making map for snapshot ', snap)
```

![enter image description here](https://github.com/rajbooth/Lightcone/raw/master/images/Galaxy_Shell_snap=60.png)
Having generated a healpix sky map, we can use the *anafast* function within healpy to calculate the angular power spectrum $(C_l)$ associated with this matter distribution.
```python
cl = hp.anafast(m[snap-44])
ell = np.arange(len(cl))
fig, ax = plt.subplots(1, 1,figsize=(12,8))
ax.loglog(ell, cl * ell * (ell+1) /(2 * np.pi), label = 'snap {0:02d}'.format(snap))
ax.set_xlabel('l')
ax.set_ylabel('$C_l.l(l+1)/2\pi$')
ax.legend()
ax.set_ylim(1e-5,1e-1)
ax.set_xlim(10, 4e3)
plt.title('Matter Angular Power Spectrum')
plt.savefig('Angular_Power_Spectrum_snap-{0:02d}.png'.format(snap))
plt.show()
```
![Angular power spectrum](https://github.com/rajbooth/Lightcone/raw/master/images/Angular_Power_Spectrum_snap-62.png)
## <a name="power_spectrum"></a>Matter power spectrum
The power spectrum for each Gadget snapshot can be read from the hdf5 data file and plotted using a code snippet such as:
```python
def Pk_s(snap):
    fname = '/cosma6/data/dp004/dc-boot5/Lightcone/Power_Spectrum/powerspec_{0:03d}.npy'.format(snap)
    with h5py.File(fname,'r') as f:
        Pk = f['Pk0'][()]	# get 3D matter power spectrum P(k)
        k = f['k'][()]  	# get k numbers from powerspec file
    return k, Pk

for snap in range(63,42,-5):
    k, Pk = Pk_s(snap)
    plt.loglog(k, Pk, label = 'Snap {0:d}'.format(snap))
plt.xlabel('k')
plt.ylabel('$P(k) [Mpc^3 h^{-3}]$')
plt.title('Gadget snapshot matter power spectrum')
plt.show()
```
![Matter power spectrum](https://github.com/rajbooth/Lightcone/raw/master/images/gadget_snapshot_Pk_1.png)
This can be compared to the power spectrum calculated by evolving back the initial $z=0$ power spectrum according to the methodology described in  [Smith, R.E>](https://doi.org/10.1093/mnras/stz890)
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTg5MjE5NTkyMCwyMjA1MjkzMjldfQ==
-->