# jarvis

Software to carry out Principal Component Analysis (PCA) of experimental (or simulated) data.

This software was originally developed by Jorge Quintanilla as part of two collaborations involving Stuart Gibson, Robert Twyman, Dylan Barker, Gunnar Moller and Tymoteusz Tula. For licensing information, see COPYING. You can cite this code using the following DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4266743.svg)](https://doi.org/10.5281/zenodo.4266743)

This software is optimised for Octave but can be easily adapted for use with Matlab.

## Requirements

A working installation of Octave (the software has been tested on Octave 3.4.3; the source code has comments explaining how it can be adapted for use on Matlab).

## Installation

Copy the .m files to your working directory.

(Alternatively, copy the .m files to a directory of your choice. In that case you need to use the addpath command within Octave to tell Octave where to look for the .m files.)  

## Usage

cd to your working directory.

Start Octave.

The following Octave functions will be available:

PCA, gen_A, gen_image_matrix, recons, scores

The usage of each of these functions is documented in comments within its source code.

### Example

The following example illustrates the use of the Octave functions in https://github.com/quintquant/jarvis to analyse data. We will use two-dimensional neutron scattering images as an example, but the data could be of any kind, e.g. a time series or a five-dimensional data set (say, temperature and pressure measured at each point in a volume). We will start with 5 observations and show that they can be described efficiently by their two coordinates in a two-dimensional space defined by two principal components.

#### Generation of the simulated observation data

For this example, the neutron scattering data will be simulated. To perform the simulation we will also need the Octave functions in https://github.com/quintquant/magneto.

Unless otherwise specified, all displayed code is Octave code. We will assume that we are issuing commands within an Octave session with both sets of functions available to Octave.

Let us start by setting the model parameters for the neutron scattering simulations. The model is described in detail in

> H. R. Irons *et al.*, *Phys. Rev. B* **96**, 224408 (2017),
> https://doi.org/10.1103/PhysRevB.96.224408.

We will simulate `num_obs=5` observations the neutron scattering function of a ring-shaped magnetic cluster with `N=2` atoms and anisotropic, nearest-neighbour magnetic interactions with parameters `Gamma=0.6` and `Delta=0.0`:
```
num_obs=5;
N=2*ones(1,num_obs);
Gamma=0.6*ones(1,num_obs);
Delta=0.0*ones(1,num_obs);
```
`num_obs` values of the magnetic field `h` and temperature `T` will be chosen pseudo-randomly from the intervals [0,2] and [0.01,1.99], respectively:
```
h=2*rand(1,num_obs);         
T=0.01+1.99*rand(1,num_obs);
```
So in this example when we say "one observation" we mean one diffuse magnetic neutron scattering function *S*(**q**)  correpsonding to one value of `T` and one value of `h`.
Finally, we will compute *S*(**q**) as a function of the two-dimensional **q** vector on an 8 by 8 grid. This is set by the parameter `Nk=4` (the grid is `2*Nk` by `2*Nk`). The values of each of the two vector components *q*<sub>x</sub> and *q*<sub>y</sub> will range from `-4.5` to `4.5`, which is set by `q_max`:
```
Nk=     4*ones(1,num_obs);
q_max=  4.5*ones(1,num_obs);
```

We now generate the `A` matrix by using the function `gen_A` from `quintquant/jarvis`, calling `gen_pic_mat_wrapped` from `quintquant/magneto`. We will generate a matrix `A` with 5 columns, one for each observation (your output will differ from mine due to the use of pseudo-random numbers):
```
octave:> [A,imax,jmax]=gen_A("gen_pic_mat_wrapped",[N;Gamma;Delta;h;T;q_max;Nk]);
octave:> A
A =

   0.93039   0.92494   0.98081   0.97464   0.89701
   0.92126   0.91591   0.97539   0.96902   0.88523
   0.91097   0.90573   0.96928   0.96268   0.87195
   0.90336   0.89820   0.96476   0.95800   0.86213
   0.90336   0.89820   0.96476   0.95800   0.86213
   0.91097   0.90573   0.96928   0.96268   0.87195
   0.92126   0.91591   0.97539   0.96902   0.88523
   0.93039   0.92494   0.98081   0.97464   0.89701
   0.71382   0.68762   0.93486   0.90661   0.56847
   0.67063   0.64488   0.90921   0.88001   0.51271
   0.60795   0.58286   0.87200   0.84141   0.43182
   0.54768   0.52323   0.83621   0.80430   0.35404
   0.54768   0.52323   0.83621   0.80430   0.35404
   ...
```
Each observation has been generated from an individual call to `gen_pic_mat_wrapped`. The latter function outputs a matrix representing the intensity of neutron scattering at different points in reciprocal space i.e. *S*(**q**). To speed up the calculation and keep the output small we have evaluated *S*(**q**) on an 8 x 8 grid in reciprocal space. However, the reader is encouraged to try also with finer grids. `gen_A` has converted that matrix into a single column of `A`, corresponding to one observation (one value of `T` and `h`).

`gen_A` has also given as its output the indices `imax` and `jmax` representing the size of the scattering matrix:
```
octave:> imax
imax =  8
octave:> jmax
jmax =  8
```
We can use these values to reconstruct the matrix corresponding to each observation from the corresponding column in the matrix `A`. For instance, take the third column of `A`:
```
octave:> A(:,3)
ans =

   0.98081
   0.97539
   0.96928
   0.96476
   0.96476
   0.96928
   0.97539
   0.98081
   0.93486
   0.90921
   0.87200
   0.83621
   0.83621
   0.87200
   ...
```
Using the function `gen_image_matrix` we can convert it back to its original matrix form:
```
octave:> Am=gen_image_matrix(A(:,3),imax,jmax)
Am =

   0.98081   0.93486   0.98727   1.01193   1.01193   0.98727   0.93486   0.98081
   0.97539   0.90921   0.98119   1.01427   1.01427   0.98119   0.90921   0.97539
   0.96928   0.87200   0.96812   1.02208   1.02208   0.96812   0.87200   0.96928
   0.96476   0.83621   0.94591   1.07286   1.07286   0.94591   0.83621   0.96476
   0.96476   0.83621   0.94591   1.07286   1.07286   0.94591   0.83621   0.96476
   0.96928   0.87200   0.96812   1.02208   1.02208   0.96812   0.87200   0.96928
   0.97539   0.90921   0.98119   1.01427   1.01427   0.98119   0.90921   0.97539
   0.98081   0.93486   0.98727   1.01193   1.01193   0.98727   0.93486   0.98081
```
This can now be plotted, e.g. by saving the matrix to a temporary file from Octave...
```
octave:> save "tmp.dat" Am;
```
...and then using GNUplot (note the following line has to be executed within a GNUplot session, not an Octave session; obviously you can use your favourite plotting program instead):
```
gnuplot> set pm3d;unset surface;splot "tmp.dat" matrix
```
![Plot of scattering function](/EXAMPLE_PLOT.png)

#### Princial component analysis of the simulated data

Let us now Carry out the Principal Component Analysis (PCA) of the data we have generated:
```
[X,U,S,scree_data]=PCA(A);
```

This has generated the following matrices:
* `X` is the same as `A` except each observation has had its average substracted (we substract the average *along* the column, not *across* columns - which is not the same as other PCA implementations - for the standard PCA analsysis, use functon PCA_std instead of PCA):
```
octave:> X
X =

   0.02220874   0.02135174   0.01542678   0.01544632   0.02733571
   0.01307913   0.01231798   0.01000546   0.00982413   0.01555249
   0.00279520   0.00214203   0.00389868   0.00349110   0.00227945
  -0.00481490  -0.00538818  -0.00062034  -0.00119535  -0.00754260
  -0.00481490  -0.00538818  -0.00062034  -0.00119535  -0.00754260
  ...
```
The substraction of the average is equivalent to substracting a uniform background from each image [in the case of *S*(**q**) considered here the average, taken in this way, should always have the same value according to a well-known sum rule - any variation is due to the finite grid and finite domain of integration].
* `U` contains the principal components. The Mth column of `U` is the Mth principal component. The total number of principal components is the same as the number of observations. The principal components are given in order of increasing variance in the data set that they capture: U(:,1) is the PC capturing the greatest amount of variance in the data set, U(:,2) is the PC that is second best at capturing the variance in the data set, etc.
* `S` gives the scores of the principal components. In our case it is a `5 x 5` (square) matrix. `S(i,j)` is the score of the jth observation on the ith principal component, i.e. the amplitude corresponding to the ith PC in the decomposition of the jth observation in principal components. Mathematically, `X(:,m) = U*S(:,m)` i.e.
```
X(:,m) = U(:,1)*S(1,m) + U(:,2)*S(2,m) + ... + U(:,M)*S(M,m)
```
In our case it takes the following form:
```
octave:209> S
S =

  -1.6124269611185  -1.7426423831502  -0.4294089049551  -0.5755247653701  -2.3943811716927
   0.0172974703210  -0.0160780470129   0.1293995132061   0.1048207822820  -0.0483485595805
   0.0000000154967   0.0000000150678   0.0000000101576   0.0000000102984   0.0000000194342
   0.0000000081819   0.0000000087439   0.0000000025331   0.0000000032003   0.0000000119396
   0.0000000370701   0.0000000404351   0.0000000085397   0.0000000121781   0.0000000558380
```
Each column represents an observation, and the five numbers in that column represent the coordinates of that observation along the five axes defined by the five principal components. Note that only the first two elements of each column have numbers of any significant size. That's because PCA chooses the axes in such way that most of the variance is captured by by just a few principal components.
* `scree_data` contains the variances (or, more precisely, squared deviation from the origin of coordinates) of the data set along the axes defined by the principal components:
```
octave:> scree_data
scree_data =

   1.1885e+01
   3.0627e-02
   1.9231e-15
   4.0773e-16
   1.3733e-17
```
We can also compute these variances "by hand" as the sums of the squares of the coordinates of our observations along the axes defined by the principal components. Let us do this for the first two principal components:
```
octave:> sum(S(1,:).^2)
ans =  11.885
octave:> sum(S(2,:).^2)
ans =  0.030627
```
This indeed coincides with the first two numbers obtained above.

Note that the first variance is large because we did not re-center the data on the average across observations - therefore, the first principal component is essentially that average. Our variances are defined as the squared deviations from the origin. We could of course have carried out the substraction when we built the matrix A, but then the data we would be working with would look quite different from the raw experimental data.

The main reason PCA is useful for our data set is because of the rapid fall in variance. In our example, only the first two principal components matter. This means that we can reconstruct any observation quite accurately with just the first two principal components. The function `recons` will reconstruct the `X` matrix using a specified number of principal components. Let us do this using 1, 2, and 3 principal components and calculate the squared deviation from the original `X` matrix:
```
octave:> X1=recons(U,S,1);
octave:> X2=recons(U,S,2);
octave:> X3=recons(U,S,3);
```
Let us now check how much each of these reconstructed `X` matrices differ from the original one:
```
octave:> sum(sum((X1-X).^2))
ans =  0.030627
octave:> sum(sum((X2-X).^2))
ans =    6.0281e-29
octave:> sum(sum((X3-X).^2))
ans =    7.6459e-29
```
We can see that from only two principal components the reconstruction is essentially perfect. This justifies describing our 5-observation data set by points in a two-dimensional plane, using the scores of each observation against the two principal components. In this representation, our scores are simply given by
```
octave:> S(1:2,:)
ans =

  -1.612427  -1.742642  -0.429409  -0.575525  -2.394381
   0.017297  -0.016078   0.129400   0.104821  -0.048349
```
i.e. they define 5 points on a two-dimensional plane.

## To do

* Add a flag to select the form of averaging (over observations vs within each observation - currently only the latter is implemented)
* Test on Matlab and resolve any incompatibilities (the software should be compatible with both Octave and Matlab)
* Improve documentation (including providing a muon spin relaxation example)
* Make into a package.
