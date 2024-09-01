# equipotential_surface

Code Information
================

***equipotential_surface*** generates an equipotential surface using a global gravity model and measurements in spherical coordinates. The code is written in $\mathrm{C}$ with the source file name **equipotential_surface.c**.

Instructions
---

To run the code, set the global gravity model parameters:

- The number of fully normalized dimensionless spherical harmonic coefficients $\left(\bar{C}_{nm}, \bar{S}_{nm}\right)$, which is dependent on the maximum degree $N$

- Geocentric gravitational  constant $GM \left(\mathrm{s^{-2}}\right)$

- Gravitational constant $G \left(\mathrm{kg^{-1}\;s^{-2}}\right)$

- Reference radius of the global gravity model $R \left(\mathrm{m}\right)$

- Angular velocity $\omega \left(\mathrm{1 / s}\right)$

- Density of topgraphy $\sigma \left(\mathrm{kg/m^{3}}\right)$

- Convergence tolerance, which is used to end the iterative process

These parameters are defined as macros in lines **25-31**, and you can set their values based on your model. In our example, we use the EGM2008 gravity model parameters:

```c
#define SCOEFF 2401333 /* Number of the the fully normalized dimensionless spherical harmonic coefficients (Cnm, Snm) */
#define GM 0.3986004415e15 /* Geocentric gravitational constant (s^-2) */
#define G 6.673e-11 /* Gravitational constant (kg^-1 s^-2) */
#define R 0.63781363e7 /* Reference radius of the GGM (m) */
#define OMEGA 7.29211500000e-5 /* Angular velocity (1 / s) */
#define DENSITY 2670.0 /* Density (kg / m^3) */
#define CONVTOL 1e-5 /* Convergence tolerance */
```

Furthermore, for the maximum degree $N$, the names of the input files and their headers (header lines) must be specified in lines **64-68**.

For example:

```c
N = 2190; /* The maximum degree of the ALFs */
char file1[BUFSIZ] = "EGM2008.gfc"; /* The name of the global gravity model file */
header1 = 22; /* The header of the global gravity model file */
char file2[BUFSIZ] = "spherical_coords_geoid_step_10.txt"; /* The name of the measurements file */
header2 = 1; /* The header of the measurements file */
```

Input file format
---------

Two files are required to run the code: one for the global gravity model and one for the measurements. The global gravity model file has the coefficients $\bar{C}_{nm}$ and $\bar{S}_{nm}$ whereas the measurements file must have the spherical coordinates (longitude $\lambda$, co-latitude $\theta$ and radius $r$) for each point.

For example, we utilize the [ICGEM's](https://icgem.gfz-potsdam.de/home) [EGM2008](https://icgem.gfz-potsdam.de/getmodel/gfc/c50128797a9cb62e936337c890e4425f03f0461d7329b09a8cc8561504465340/EGM2008.gfc) model, and the input file snippets are shown below:

<center><b>Global gravity model file:</b></center>

```
Pavlis, N.K., S.A. Holmes, S.C. Kenyon, and J.K. Factor:
An Earth Gravitational Model to Degree 2160: EGM2008,
presented at the 2008 General Assembly of the European Geosciences Union,
Vienna, Austria, April 13-18, 2008.


product_type                gravity_field
modelname                   EGM2008
earth_gravity_constant      0.3986004415E+15
radius                      0.63781363E+07
max_degree                  2190
errors                      calibrated
norm                        fully_normalized
tide_system                 tide_free

url                   http://earth-info.nima.mil/GandG/



key     L    M             C                       S                    sigma C             sigma S
end_of_head ============================================================================================
gfc     0    0    1.0d0                    0.0d0                    0.0d0               0.0d0
gfc     2    0   -0.484165143790815e-03    0.000000000000000e+00    0.7481239490e-11    0.0000000000e+00
gfc     2    1   -0.206615509074176e-09    0.138441389137979e-08    0.7063781502e-11    0.7348347201e-11
gfc     2    2    0.243938357328313e-05   -0.140027370385934e-05    0.7230231722e-11    0.7425816951e-11
gfc     3    0    0.957161207093473e-06    0.000000000000000e+00    0.5731430751e-11    0.0000000000e+00
gfc     3    1    0.203046201047864e-05    0.248200415856872e-06    0.5726633183e-11    0.5976692146e-11
gfc     3    2    0.904787894809528e-06   -0.619005475177618e-06    0.6374776928e-11    0.6401837794e-11
gfc     3    3    0.721321757121568e-06    0.141434926192941e-05    0.6029131793e-11    0.6028311182e-11
gfc     4    0    0.539965866638991e-06    0.000000000000000e+00    0.4431111968e-11    0.0000000000e+00
gfc     4    1   -0.536157389388867e-06   -0.473567346518086e-06    0.4568074333e-11    0.4684043490e-11
gfc     4    2    0.350501623962649e-06    0.662480026275829e-06    0.5307840320e-11    0.5186098530e-11
gfc     4    3    0.990856766672321e-06   -0.200956723567452e-06    0.5631952953e-11    0.5620296098e-11
gfc     4    4   -0.188519633023033e-06    0.308803882149194e-06    0.5372877167e-11    0.5383247677e-11
```

<center><b>Measurements file:</b></center>

```
Longitude      Co-latitude    radius
-180.000000    10.066021      6357406.238    
-170.000000    10.066021      6357405.262    
-160.000000    10.066021      6357403.199    
-150.000000    10.066021      6357402.968    
-140.000000    10.066021      6357401.725    
-130.000000    10.066021      6357400.421    
-120.000000    10.066021      6357402.978    
-110.000000    10.066021      6357408.320    
-100.000000    10.066021      6357406.319    
-90.000000     10.066021      6357411.476    
-80.000000     10.066021      6357414.724    
-70.000000     10.066021      6357409.845
```

Run the code
---------------

To compile the code, type the following command at the Linux command line:

```bash
gcc -Wall -O3 equipotential_surface.c -o equipotential_surface -lm
```

This program creates a single executable with the name *equipotential_surface*. Execute the application using the following command on the Linux command line:

```bash
./equipotential_surface
```

Below is an excerpt of the outcomes:

```
N = 2190
S = 2401336 ALFs
Points = 612
Degrees of freedom = 611

Computation of the gravity potential and its derivatives for each point: Start (2s)
Computation of the gravity potential and its derivatives for each point: End (11s)
Computation of the elements of the matrix B
Inversion of the matrix M
Computation of the vector w
Computation of the radial distances di
iteration 0 completed
sigma0 = +/- 0.29496922    

Computation of the gravity potential and its derivatives for each point: Start (11s)
Computation of the gravity potential and its derivatives for each point: End (19s)
Computation of the elements of the matrix B
Inversion of the matrix M
Computation of the vector w
Computation of the radial distances di
iteration 1 completed
sigma0 = +/- 0.00000002    

Computation of the gravity potential and its derivatives for each point: Start (19s)
Computation of the gravity potential and its derivatives for each point: End (27s)
Computation of the elements of the matrix B
Inversion of the matrix M
Computation of the vector w
Computation of the radial distances di
iteration 2 completed
sigma0 = +/- 0.00000010    
n = 612 measurements
r = 611 degrees of freedom
iterations = 3
sigma0 = +/- 0.0000    

Longitude      Co-latitude    Radius              di                  W0

-180.0000      10.0660        6357406.161         0.000               62636852.475        
-170.0000      10.0660        6357405.184         0.000               62636852.475        
-160.0000      10.0660        6357403.122         0.000               62636852.475        
-150.0000      10.0660        6357402.891         0.000               62636852.475        
-140.0000      10.0660        6357401.648         -0.000              62636852.475        
-130.0000      10.0660        6357400.344         -0.000              62636852.475        
-120.0000      10.0660        6357402.901         0.000               62636852.475        
-110.0000      10.0660        6357408.242         0.000               62636852.475        
-100.0000      10.0660        6357406.242         0.000               62636852.475        
-90.0000       10.0660        6357411.427         0.000               62636852.475        
-80.0000       10.0660        6357414.798         0.000               62636852.475        
-70.0000       10.0660        6357409.768         0.000               62636852.475        
-60.0000       10.0660        6357422.143         0.000               62636852.475
```

The final spherical coordinates determine the equipotential surface, and the potential value $W_{0}$ is constant for all points, as shown in the last column of the findings.

*Individuals who are unfamiliar with the Linux command line may find it easier to compile and execute code using an Integrated Development Environment (IDE).
