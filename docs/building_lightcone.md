## Methodology

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
The redshift space distortion for a galaxy at redshift $z$ is given, to a good approximation, by

$$r = (1+z) \frac{v_\parallel }{H(z)}$$
where $v_\parallel$ is the peculiar velocity of the galaxy in the line-of-sight direction.
This redshift distortion factor is calculated for every galaxy in the lightcone, based on the radial peculiar velocity of the corresponding particle in the simulation snapshot.
<!--stackedit_data:
eyJoaXN0b3J5IjpbODg1MzE1MTI4LDE3MzAwNTQ1MDksODEwOT
IzMDUyXX0=
-->