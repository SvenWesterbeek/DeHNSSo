#
# DeHNSSo

Delft Harmonic Navier-Stokes Solver - A nonlinear stability solver for flow problems in complex 2-dimensional domains.

_Thank you for your interest in using DeHNSSo. This page is intended to guide you through using DeHNSSo easily and fast._

© 2023, Sven Westerbeek & Marios Kotsonis

## Table of Contents

1. [How to use DeHNSSo](#how-to-use-dehnsso)
2. [Installation](#installation)
3. [Input files](#input-files)
   1. [BF](#bf)
   2. [Grid](#grid)
   3. [Stab](#stab)
   4. [Opt](#opt)
4. [Output files](#output-files)
   1. [StabGrid](#stabgrid)
   2. [StabRes](#stabres)
   3. [BF](#bf2)
5. [Example cases](#example-cases)
   1. [Blasius boundary layer](#blasius-boundary-layer-tollmien-schlichting-instabilities)
   2. [Swept-wing boundary layer: stationary crossflow instability](#swept-wing-boundary-layer-stationary-crossflow-instability)
   3. [Swept-wing boundary layer: interaction of a stationary CFI with a hump](#swept-wing-boundary-layer-interaction-of-a-stationary-cfi-with-a-hump)
   4. [Swept-wing boundary layer: interaction of a stationary CFI with a step](#swept-wing-boundary-layer-interaction-of-a-stationary-cfi-with-a-step)
6. [License](#license)
7. [Common errors and solutions](#common-errors-and-solutions)
8. [FAQ](#faq)
9. [Data visualization examples](#data-visualization-examples)
10. [DeHNSSo results](#dehnsso-results)
11. [How to cite DeHNSSo](#how-to-cite-dehnsso)
12. [References](#references)

## How to use DeHNSSo <a id="how-to-use-dehnsso"></a>

The first step to using DeHNSSo is ensuring that you have MATLAB installed. It is recommended to use at least version R2022b. For help with installing Matlab, please follow the instructions on the [MathWorks web page](https://nl.mathworks.com/products/matlab.html). Then, you will need to download DeHNSSo from the DeHNSSo GitHub as will be explained below. DeHNSSo comes with 4 standard examples that can be run immediately. To perform custom simulations, please familiarize yourself with the input formats through these examples and adjust them accordingly.

## Installation <a id="installation"></a>

To start using DeHNSSo, please download the files from this GitHub [Repository](https://github.com/SvenWesterbeek/DeHNSSo). The current version of DeHNSSo was made using MATLAB R2022b. Support for earlier versions is not guaranteed. This file contains several folders, namely: Callers, Tools, Data Files, and Documentation. All folders must remain together in a directory of your choosing.

The base flow data for the fourth example case "CFI over a step in a swept-wing BL" is too large for this repository and can be found in the data [Repository](https://PLACEHOLDER) repository of J. Casacuberta instead.

### Input files <a name="input-files"></a>

DeHNSSo requires several inputs that describe the basic state flow field, the numerical domain and grid, and an inflow boundary condition among others. This data should be provided in dimensionless form. The reference values should be provided in BF. The contents of these structs should carry specific names. These structs can be summarized as follows:

1. BF - Base flow domain, grid, velocity field, and reference values
2. Grid - numerical domain and discretization
3. Stab - Mode specifications, spectral truncation, and inflow data
4. Opt - Solver options such as outflow buffer specifications and inflow amplitude growth rate.


These structures will be discussed in detail below to help you make your own input files.

### BF <a id="bf"></a>

BF contains the base flow data, the grid on which it is defined, and the reference values for nondimensionalization. This data can be provided in a structured, unstructured grid, or vector format.

This data does not need to be presented on the same grid as presented in the Grid struct as the data will be interpolated onto the numerical grid within DeHNSSo using the griddata function (method = 'cubic'). The numerical domain cannot exceed the domain of the base flow. In the table below, the required contents of BF are described. In short, BF.X and BF.Y describe the locations where base flow quantities are described. BF.U, BF.V and BF.W describe the velocities in $x$, $y$ and $z$ respectively. Then, the BF.dxU, BF.dxV, BF.dxW, BF.dyU, BF.dyV, and BF.dyW describe the streamwise and wall-normal derivatives of the aforementioned velocities. Lastly, BF.lref, BF.Uref, BF.nu and BF.Re are the reference values and corresponding Reynolds number defined as .

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| _BF.X_ | Streamwise grid locations | [-] | $(nx,ny)$ |
| _BF.Y_ | Wall-normal grid locations | [-] | $(nx,ny)$ |
| _BF.U_ | Streamwise velocity | [-] | $(nx,ny)$|
| _BF.V_ | Wall-normal velocity | [-] | $(nx,ny)$ |
| _BF.W_ | Spanwise velocity | [-] | $(nx,ny)$ |
| _BF.dxU_ | First-order streamwise derivative of BF.U | [-] | $(nx,ny)$ |
| _BF.dxV_ | First-order streamwise derivative of BF.V | [-] | $(nx,ny)$ |
| _BF.dxW_ | First-order streamwise derivative of BF.W | [-] | $(nx,ny)$ |
| _BF.dyU_ | First-order wall-normal derivative of BF.U | [-] | $(nx,ny)$ |
| _BF.dyV_ | First-order wall-normal derivative of BF.V | [-] | $(nx,ny)$ |
| _BF.dyW_ | First-order wall-normal derivative of BF.W | [-] | $(nx,ny)$ |
| _BF.lref_ | Reference length | [m] | $1$ |
| _BF.Uref_ | Reference velocity | [m/s] | $1$ |
| _BF.nu_ | Kinematic viscosity | [m^2/s] | $1$ |
| _BF.Re_ | Reynolds number | [-] | $1$ |

### Grid <a id="grid"></a>

The Grid structure contains the information on the numerical domain and grid. All contents of this structure are summarized in the table below.

_Grid.wall_ presents the CHNS with the bottom wall coordinates via 2 rows of data. The first row contains the $x$-coordinates and the second row supplies the corresponding $y$-coordinates; This matrix can be of any size. However, it is preferred to be highly refined around any geometric wall features to ensure the interpolation is performed well. Sharp features should not be accounted for in wall as they will be incorporated via an embedded boundary method.

_Grid.mode_ allows for several built-in grid generations to be used. The options currently available are:

- "equidistant"
  - Creates a grid following equidistant streamwise discretization and an eta distribution of collocation points following Malik for a user-defined median collocation point $y\_i$.
- "refined"
  - Creates a grid with a streamwise refined grid based on a Gaussian distribution following user-defined inputs _Grid.mug, Grid.sig and Grid.ag._ Wall-normal distribution follows Malik for a user-defined median collocation point $y\_i$.
- "curved"
  - Creates an equidistant distribution of streamwise grid points over a curved surface with straight $\eta$ axes wall-normal to that surface. The $\eta$ axes are thus not parallel to the global y axes. Wall-normal collocation points are clustered near the wall following Malik's mapping for a user-defined median collocation point $y\_i$.
- "wallorthogonal"
  - Creates a grid elliptically with orthogonality at the wall following exactly the eta distribution of the mapping of Malik according to a user-defined $y\_i$. The streamwise distribution is nearly equidistant but can be slightly adjusted for the sake of orthogonality.

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| _Grid.nx_ | Number of streamwise stations | [-] | $1$ |
| _Grid.ny_ | Number of wall-normal stations | [-] | $1$ |
| _Grid.wall_ | Smooth wall description | [-] | ($nx\_{wall}$,2) |
| _Grid.H_ | Domain height | [-] | $1$ |
| _Grid.y\_i_ | Median collocation point height | [-] | $1$ |
| _Grid.S_ | Domain start | [-] | $1$ |
| _Grid.L_ | Domain length | [-] | $1$ |
| _Grid.mode_ | Grid generation mode | [-] | String |
| _Grid.ft_ | Flat top flag [0, 1], enforces the top boundary to be flat if 1 (default = 1) | [-] | $1$ |
| _Grid.mug_ | Refinement peak location [S S+L] | [-] | $1$ |
| _Grid.sig_ | Refinement variance (Gaussian) [0-1] (default = 1) | [-] | $1$ |
| _Grid.ag_ | Refinement strength [0-1] (default = 0) | [-] | $1$ |
| _Grid.StepX_ | Step location (default = 0) | [-] | $1$ |
| _Grid.StepH_ | Step Height (default = 0) | [-] | $1$ |
| _Grid.ystretch_ | Wall-normal distribution stretching factor | [-] | $1$ |
| Grid.StepType | Sharp geometry type "FFS", "BFS", "GAP", "HUMP" | [-] | String |

## <a id="stab"></a> stab

The _Stab_ structure is used to define the mode ensemble of interest and present the solver with inflow conditions.

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| _Stab.N_ | Spectral truncation of beta modes | [-] | $1$ |
| _Stab.M_ | Spectral truncation of omega modes | [-] | $1$ |
| _Stab.A0_ | Initial amplitude of all modes | [-] | $((2N+1) \times (2M+1),1)$ |
| _Stab.omega\_0_ | Fundamental frequency | [-] | $1$ |
| _Stab.beta\_0_ | Fundamental spanwise wavelength | [-] | $1$ |
| _Stab.IC_ | Initialization method "ILST","ZERO'', 'LOAD" | [-] | string |
| _Stab.bcwx_ | Inhomogeneous boundary condition locations | [-] | $(any,1)$ |
| _Stab.bcw_ | Inhomogeneous boundary conditions (default = 0's) | [-] | $($ Stab.bcwx $, 3\times(2N+1)\times(2M+1))$ |
| _Stab.u0_ | Normalized streamwise perturbation velocity at x\_0 | [-] | $(3 \times (2N+1) \times (2M+1)),ny)$ |
| _Stab.v0_ | Normalized wall-normal perturbation velocity at x\_0 | [-] | $(3 \times (2N+1) \times (2M+1)),ny)$ |
| _Stab.w0_ | Normalized spanwise perturbation velocity at x\_0 | [-] | $(3 \times (2N+1) \times (2M+1)),ny)$ |
| _Stab.p0_ | Normalized perturbation pressure at x\_0 | [-] | $(3 \times (2N+1) \times (2M+1)),ny)$ |
| _Stab.y0_ | Wall-normal distribution of inflow perturbation data | [-] | $(1,ny)$ |

_Stab.IC_ sets the mode initialization method. Currently, two methods are implemented. "ILST" calls a routine that finds the solution to the local eigenvalue problem at the inflow for all modes that have a nonzero amplitude (presented in _Stab.A0_). These results are then normalized with the maximum streamwise perturbation velocity and multiplied by the respective initialization amplitude. "ZERO" instead means no inflow condition is supplied. This generally means that the user intends to simulate the receptivity problem by supplying inhomogeneous boundary conditions. "LOAD" instead uses the perturbation profiles presented in _Stab.u0, Stab.v0, Stab.w0, Stab.p0_ to define the inflow boundary condition. The initial perturbation data is interpolated onto the numerical grid within the solver and can thus be supplied on any distribution of points consistent with _Stab.y0_.

Stab.bcw is used to define inhomogeneous wall conditions in the streamwise, wall-normal and spanwise velocity components per mode defined at the streamwise locations presented in _Stab.bcwx._

## Opt <a id="opt"></a>

The final structure put into the solver contains solver-specific options.

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| Opt.xb | Buffer starting location [0 1] | [-] | $1$ |
| Opt.kappa | Buffer strength [1 -\>] | [-] | $1$ |
| Opt.nltbufxb | Nonlinear term buffer starting location | [-] | $1$ |
| Opt.Th | Nonlinear convergence threshold | [-] | $1$ |
| Opt.Sweep | Output intermediate results flag (true = $1$, false = $0$) | [-] | $1$ |
| Opt.AFg | Amplitude factor growth rate (default = 1.1) | [-] | $1$ |
| Opt.Conv | Convergence criterion (default = 1e-4) | [-] | $1$ |
| Opt.ConvF | Convergence criterion relaxation during ramping (default = 100) | [-] | $1$ |
| Opt.AMAX | Maximum amplitude for initializing ramping procedure (default = 0.1) | [-] | $1$ |

## output files <a id="output-files"></a>

The solver returns the following structures with outputs:

1. StabGrid; The numerical grid generated within the solver as well as all grid transformations
2. StabRes; Stability calculation results
3. BF; Base flow values interpolated on the numerical grid

## StabGrid <a id="stabgrid"></a>

The StabGrid structure contains the numerical grid generated in the solver on which the simulation results are defined. The streamwise location increases with the column index while the wall-normal location decreases with the row index. The StabGrid contains both the physical ($x$, $y$) and the computational ($\xi$, $\eta$) grid. The transformation coefficients are also presented in this structure.

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| StabGrid.xw | Global Streamwise wall coordinate | [-] | $(nx)$ |
| StabGrid.yw | Global Wall-normal wall coordinate | [-] | $(nx)$ |
| StabGrid.x | Global streamwise coordinate | [-] | $(nx,ny)$ |
| StabGrid.y | Global wall-normal coordinate | [-] | $(nx,ny)$ |
| StabGrid.xi | Computational streamwise coordinate | [-] | $(nx,ny)$ |
| StabGrid.eta | Computational wall-normal coordinate | [-] | $(nx,ny)$ |
| StabGrid.xix | dxi/dx transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.xiy | dxi/dy transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.xixx | ddxi/dxx transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.xiyy | ddxi/dyy transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.etax | deta/dx transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.etay | deta/dy transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.etaxx | ddeta/dxx transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.etayy | ddeta/dyy transformation coefficient | [-] | $(nx,ny)$ |

## StabRes <a id="stabres"></a>

The StabRes structure contains all the stability results defined on the locations defined by StabGrid.x, StabGrid.y. Additionally, some key factors are calculated that are commonly used in stability analysis for comparison purposes even when they might not be used in HNS (such as the streamwise wavenumber $\alpha$). Results are shown nondimensionally and normalized. In other words, amplitudes are extracted from the perturbation shape functions of $\hat{u}$, $\hat{v}$, $\hat{w}$, and $\hat{p}$ based on the maximum of the absolute streamwise velocity perturbation value $u$.

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| StabRes.A | Perturbation amplitudes based on maximum u | [-] | $(nx,nf)$ |
| StabRes.u | Streamwise perturbation velocities | [-] | $(nx,ny,nf)$ |
| StabRes.v | Wall-normal perturbation velocities | [-] | $(nx,ny,nf)$ |
| StabRes.w | Spanwise perturbation velocities | [-] | $(nx,ny,nf)$ |
| StabRes.p | Perturbation pressures | [-] | $(nx,ny,nf)$ |
| StabRes.beta | Spanwise wavenumber per mode | [-] | $(nf)$ |
| StabRes.omega | Angular frequency per mode | [-] | $(nf)$ |
| StabRes.alpha | Streamwise wavenumber per mode | [-] | $(nx,nf)$ |

If an amplitude sweep is performed, intermediate results are also presented via additional outputs in the StabRes struct. These fields will contain the word "sweep". The "iter" in the size of these fields corresponds to the iteration. Not that this concerns only one fully converged result per inflow amplitude.

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| StabRes.Asweep | Perturbation amplitudes based on maximum u | [-] | $(nx,nf,iter)$ |
| StabRes.usweep | Streamwise perturbation velocities | [-] | $(nx,ny,nf,iter)$ |
| StabRes.vsweep | Wall-normal perturbation velocities | [-] | $(nx,ny,nf,iter)$ |
| StabRes.wsweep | Spanwise perturbation velocities | [-] | $(nx,ny,nf,iter)$ |
| StabRes.psweep | Perturbation pressures | [-] | $(nx,ny,nf,iter)$ |
| StabRes.alphasweep | Streamwise wavenumber per mode | [-] | $(nx,nf,iter)$ |
| StabRes.beta | Spanwise wavenumber per mode | [-] | $(nf)$ |
| StabRes.omega | Angular frequency per mode | [-] | $(nf)$ |

## BF <a id="bf"></a>

The BF structure is both an input and output. In the output, the structure is appended with the result of the base flow interpolation for easier post-processing. The following quantities are added:

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| BF.Ur | Streamwise base flow velocity | [-] | $(nx,ny)$ |
| BF.Vr | Wall-normal base flow velocity | [-] | $(nx,ny)$ |
| BF.Wr | Spanwise base flow velocity | [-] | $(nx,ny)$ |
| BF.dxUr | $x$-derivative of streamwise base flow velocity | [-] | $(nx,ny)$ |
| BF.dxVr | $x$-derivative of wall-normal base flow velocity | [-] | $(nx,ny)$ |
| BF.dxWr | $x$-derivative of spanwise base flow velocity | [-] | $(nx,ny)$ |
| BF.dyUr | $y$-derivative of streamwise base flow velocity | [-] | $(nx,ny)$ |
| BF.dyVr | $y$-derivative of wall-normal base flow velocity | [-] | $(nx,ny)$ |
| BF.dyWr | $y$-derivative of spanwise base flow velocity | [-] | $(nx,ny)$ |

# Example cases <a id="example-cases"></a>

Some example cases are presented here so that users can test DeHNSSo and get familiar with the required inputs through some canonical cases. Additionally, these cases can be adjusted easily to represent cases relevant to the user. The details of the inputs are not disclosed here. They can be found in the respective caller. Instead, the cases are shortly described below.

## Blasius boundary layer: Tollmien-Schlichting instabilities <a id="blasius-boundary-layer-tollmien-schlichting-instabilities"></a>

The development of Tollmien-Schlichting instabilities in a Blasius boundary layer is considered in the first case. This case was previously considered in Bertolotti et al. (1992), Chang et al. (1993), and Herbert (1997).

### Geometry, outflow buffer, and discretization

The flow setup can be described by a constant external velocity of 10 m/s (and zero pressure gradient) over a domain ranging from $x = 0.243$ to $x = 1.7336$ m. The domain height is set to $H = 0.06$ m.

### Reference quantities

All quantities will be Nondimensionalized by (a combination of) the reference velocity $U\_0=10$ m/s and $l\_{ref} = 6.075 \times 10^{-4}$ m. The kinematic viscosity is $1.518$ m^2/s such that the global $Re = U\_0 \times l\_{ref}/\nu = 400$ in accordance with the references

### Nonlinear mode ensemble and amplitude ramping

The inflow conditions are comprised of the solution to the local eigenvalue problem at the inflow for $90$ Hz corresponding to $\omega = 0.0344$ superimposed with an initial streamwise perturbation amplitude of $A = 0.00125 \sqrt{2}$. This inflow amplitude does not require any amplitude ramping.

For this example case, the spectral domain is truncated at $M = 5$ (and $N = 0$). Equal to the number of modes presented in the references. Higher harmonics and the mean flow distortion are not presented at the inflow and rise naturally downstream of the inflow through nonlinear forcing.

### Base Flow

The base flow is the solution to the incompressible boundary layer equation found using an in-house solver. The base flow simulation was performed on a fine equidistant numerical grid of 2000 streamwise stations by 100 collocation points and stored here in the BF structure. The boundary layer solver is not part of the package that comes with DeHNSSo.

Running this case exactly as given will result in the data used to plot Figure 8 from Westerbeek et al. (2023).

## Swept-wing boundary layer: Stationary crossflow instability <a id="swept-wing-boundary-layer-stationary-crossflow-instability"></a>

In the second example, the stability of a swept-wing boundary layer is assessed nonlinearly for stationary crossflow instabilities. The boundary layer is simulated on a flat plate mimicking the experiments of Rius-Vidales et al. (2021). This is done by imposing a fitted external velocity distribution on the top boundary of the base flow simulation as presented in Casacuberta et al. (2022). The external velocity is given by the equation:

$U_e(x) = 0.0023 \ln(x)^4 + 0.0377 \ln(x)^3 + 0.1752 \ln(x)^2 + 0.5303 \ln(x) + 1.874$

and holds for all cases described hereafter.

### Domain description and reference values

The problem is described additionally by the streamwise coordinate ranging from $x = 220$. The domain height is set to $H = 0.02$ m. The reference length, defined as the Blasius length at the inflow is $l\_{ref} = 2.1394 \times 10^{-4}$ m. The external velocity at the inflow is $U\_{ref} = 15.1$ m/s. The kinematic viscosity of the flow is $1.47 \times 10^{-5}$ such that the global Reynolds number is $220$.

### Base Flow

The base flow data is an interpolated version of DNS results using INCA Hickel and Adams (2008), Hickel et al.(2014) for the problem as described in the previous paragraph. This data was kindly provided by J. Casacuberta who performed a full DNS simulation for this problem in J. Casacuberta et al. (2022).

### Numerical domain and buffer

The domain is discretized in 1272 equidistant stations and 50 wall-normal collocation points. The outflow buffer is initialized from 85% of the domain and an additional damping of the nonlinear terms starts from 80% of the domain.

### Nonlinear mode ensemble, reference data, and ramping

The spectral domain is truncated at five harmonics for this case ($N=5$, $M=0$) and the mean flow distortion. Only the fundamental mode, characterized by a spanwise wavelength of $\lambda_z = 7.5$ mm or $\beta = 0.18$, is introduced at the inflow as the solution to the local eigenvalue problem. The inflow amplitude of $A = 3.5 \times 10^{-2}$ is imposed on the result. For this case, amplitude ramping is required. The inflow amplitude is increased by 10% each iteration after an amplitude reduction ensures that the linear simulation result is capped at AMAX ($= 0.1$, the default value). Reference data for the stability solution is presented by J. Casacuberta (2021) as well as NPSE solutions for the current problem provided by the authors of the current code and also previously shown in the aforementioned work.

## Swept-wing boundary layer: Interaction of Stationary crossflow instability with a hump <a id="swept-wing-boundary-layer-interaction-of-a-stationary-cfi-with-a-hump"></a>

This third simulation considers the interaction of a stationary CFI with a smooth hump. This case was previously examined in both Westerbeek 2023a and the article on the current solver Westerbeek 2023b.

### Domain description, base flow, and reference values

The same external flow velocity as for the flat plate case is prescribed at the top boundary of the base flow calculation. For this problem, the base flow calculation is performed in COMSOL (COMSOL 2022). The domain is described by $x = 220-2165$ and a domain height of $H=89$. Normalized by the reference length $l\_{ref} = 2.14 \times 10^{-4}$ m with is the Blasius length at the inflow. This leads to the same reference values as for the flat plate case. The external velocity at the inflow is $U\_{ref} = 15.1$ m/s. The kinematic viscosity of the flow is $1.47 \times 10^{-5}$ such that the global Reynolds number is $220$.

The hump is symmetrical and centered around $x\_m = 859$ and relatively shallow. Consequently, the base flow features no flow separation. The problem is thus not too challenging and does not demand much refinement around the hump. Still, a slight streamwise refinement is introduced (using Grid.mode "xrefined") around the hump for users to get familiar and play with.

### Nonlinear mode ensemble, reference data, and ramping

The stability of this flow problem is assessed linearly ($N=1$, $M=0$) given that few solvers are able to perform such simulations. Only the fundamental CFI, characterized by $\beta = 0.18$ is introduced at the inflow. The reference data is provided by J. Franco of DLR using AHLNS (see Franco 2018) on the same base flow. The AHLNS is physically equivalent to the current HNS for linear simulations. No amplitude ramping is required as this concerns a linear simulation.

The default value of $x_b = 0.85$ is used to define a buffer region covering the last 15% of the domain.

## Swept-wing boundary layer: Interaction of Stationary crossflow instability with a Forward-Facing Step <a id="swept-wing-boundary-layer-interaction-of-a-stationary-cfi-with-a-step"></a>

This last simulation

# License<a id="license"></a>

The contents in the Data files directory are licensed under a  **CC-BY 4.0**  (see CC-BY-4.0 file and [here](/https://joinup.ec.europa.eu/licence/creative-commons-attribution-40-international-cc-40). The codes and any other file in this repository are licensed under a  **GPLv3 license**  (see GPLv3 file and [here](/https://www.gnu.org/licenses/gpl-3.0.en.html).

Copyright notice:

Technische Universiteit Delft hereby disclaims all copyright interest in the program "DeHNSSo". DeHNSSo is a MATLAB tool to solve Nonlinear Harmonic Navier-Stokes problems written by the Author(s).

Henri Werij, Dean of Faculty of Aerospace Engineering, Technische Universiteit Delft.

© 2023, Sven Westerbeek & Marios Kotsonis

# Common errors and solutions <a id="common-errors-and-solutions"></a>
Some errors might be the result of improper use rather than a problem with the solver. In this section, we will post some hints as to what might have gone wrong when certain errors appear. Additionally, DeHNSSo might be adjusted to provide error messages with solution hints within MATLAB directly.

# FAQ<a id="faq"></a>
Questions that we receive concerning DeHNSSo will be used to improve DeHNSSo. This concerns either the code itself, the commenting, the manual and/or posting the questions with answers in this section.

# DeHNSSo results<a id="dehnsso-results"></a>
Here we will post results obtained using DeHNSSo. If you have worked with DeHNSSo and have figures that you would like to post. Please contact us at:
S.H.J.Westerbeek@tudelft.nl
M.Kotsonis@tudelft.nl

# How to cite DeHNSSo<a id="how-to-cite-dehnsso"></a>
If you have used DeHNSSo or any subroutine within this GitHub. Please cite the GitHub and the journal article:


# References<a id="references"></a>
