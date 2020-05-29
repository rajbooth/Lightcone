## Methodology

### Source data
The data used in this lightcone catalogue originates in the N-body simulations that were carried out as part of the [Daemmerung simulations](https://doi.org/10.1093/mnras/stz890), using the Gadget 2 N-body code.  Specifically, the lightcone uses a subset of the snapshots generated for the fiducial run for 'large-box' component of these simulations. The large-box runs tracked $N=2048^3$ dark matter particles in a comoving box of size $L = 3000 Gpc/h$ , with a mass per particle of $m_p = 2.69 \times 10^{11} M_\odot/h$.
The principal cosmological parameters used for the fiducial run are summarised in the table below.

| Parameter |         Description     | Value |
|----------|-----------------------------|-----------|
| $\Omega_{DE}$ | Dark energy fraction of $\Omega$ | 0.6914 |
|$\omega_c \equiv \Omega_c h^2$	 | Cold dark matter fraction of $\Omega$  | 0.11889 |
|$\omega_b\equiv \Omega_b h^2$  | Baryonic matter fraction of $\Omega$ |0.022161 |
| $w_0$ | Dark energy EoS parameter | -1.0 |
| $w_a$ | Dark energy EoS parameter | 0.0 |
| $n_s$ | Spectral index | 0.9611 |
| $A_s$ | Curvature fluctuation amplitude | 2.14818 $(\times 10^{-9})$ |

The lightcone includes snapshots 42 - 63 from the fiducial Gadget run, spanning  the redshift range $0 < z < 0.85$

![Gadget snapshot redshifts](https://github.com/rajbooth/Lightcone/raw/master/images/Gadget%20snaphots%20vs%20redshift.png)

(Note: the Gadget snapshots in this simulation are not linearly spaced in redshift, as shown in the figure above.)

### Calculation of data fields in lightcone
The lightcone build process takes as its input the particle position and velocity data read from the Gadget snapshot files. This is first converted from Cartesian to polar coordinates:
```python
# calculate comoving radial distance, RA and Dec
r = np.sqrt(x*x + y*y + z*z)
dec = np.rad2deg(np.arcsin(z/r))
ra = np.rad2deg(np.arctan2(y,x))
```
Since the original Cartesian coordinates are in units of $Mpc/h$, the derived co-moving radial distance $r$ will also be in these units.

The redshift $z$ used in the lightcone dataset is calculated from  the comoving  radial distance $r$ by reference to a spline lookup function
$z = d2z(r)$, where $d2z()$ is defined as
```python
for d in range(1,3000):
	d_c[d] = d
	z[d] = z_at_value(cosmo.comoving_distance, d * uni.Mpc /h)    
# create spine for quick lookup
d2z = spl(d_c,z)
```
This makes use of the astropy python library function *z_at_value* to determine the value of z corresponding to the comoving distance $d_C$ for the given cosmology defined in *cosmo*. This has been initialised with the cosmological parameters used in the Gadget run (see above).  
The astropy *comoving_distance* method uses the standard formula for calculating this quantity, i.e.
$\displaystyle \chi =\dfrac{c}{H_0} \int_0^z {\dfrac{dz^\prime}{E(z^\prime)}}$$
where $E(z)$ is the Hubble parameter evolution function for $\Lambda CDM$ cosmology, defined as
$E(z) = \sqrt{\Omega_M(1+z)^3 + \Omega_k(1+z)^2 + \Omega_\Lambda}$

Note that as 
###  <a name="luminosity"></a> The luminosity function

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
The redshift space distortion for a galaxy at redshift $z$ is given, to a good approximation, by

$$r = (1+z) \frac{v_\parallel }{H(z)}$$
where $v_\parallel$ is the peculiar velocity of the galaxy in the line-of-sight direction.
This redshift distortion factor is calculated for every galaxy in the lightcone, based on the radial peculiar velocity of the corresponding particle in the simulation snapshot.
<!--stackedit_data:
eyJoaXN0b3J5IjpbMzU5NjcyNjE0LDY4NDkxMzg3MiwtNzIwMz
Y5MywtOTk0MjUxNzYsMTM4NzEyOTE1LDE5NTg3MzU1MTEsMTM2
MjcyMzYwNywtMjEwNTY0OTAzMiwxMjg5OTEzODc0LDQ1MjQ2NT
Y4NywtMTg2OTMxOTA4OSwtNTIxOTIzNDkxLDg4NTMxNTEyOCwx
NzMwMDU0NTA5LDgxMDkyMzA1Ml19
-->