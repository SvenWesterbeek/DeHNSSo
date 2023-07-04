<pre>
 ██████████            █████   █████ ██████   █████  █████████   █████████          
░░███░░░░███          ░░███   ░░███ ░░██████ ░░███  ███░░░░░███ ███░░░░░███         
 ░███   ░░███  ██████  ░███    ░███  ░███░███ ░███ ░███    ░░░ ░███    ░░░   ██████ 
 ░███    ░███ ███░░███ ░███████████  ░███░░███░███ ░░█████████ ░░█████████  ███░░███
 ░███    ░███░███████  ░███░░░░░███  ░███ ░░██████  ░░░░░░░░███ ░░░░░░░░███░███ ░███
 ░███    ███ ░███░░░   ░███    ░███  ░███  ░░█████  ███    ░███ ███    ░███░███ ░███
 ██████████  ░░██████  █████   █████ █████  ░░█████░░█████████ ░░█████████ ░░██████ 
░░░░░░░░░░    ░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░  ░░░░░░░░░   ░░░░░░░░░   ░░░░░░  
                                                                                          
 
</pre>                                                                  
                                                                          
                                                                          

Delft Harmonic Navier-Stokes Solver - A nonlinear stability solver for flow problems in complex two-dimensional domains.

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
12. [Acknowledgments](#acknowledgments)
13. [References](#references)

## How to use DeHNSSo <a id="how-to-use-dehnsso"></a>

The first step towards using DeHNSSo is ensuring that you have MATLAB installed **[#CMM:is this sufficient? what toolboxes are necessary? i found a symbolic toolbox necessity in grid. what others?]**. It is recommended to use at least version R2022b. For help with installing Matlab, please follow the instructions on the [MathWorks web page](https://nl.mathworks.com/products/matlab.html). Then, you will need to download DeHNSSo from the DeHNSSo repository **[#CMM: update with correct name, possibly GitLab]** as will be explained below. DeHNSSo comes with 4 representative test cases that can be executed immediately. To perform custom or modified simulations, please familiarize yourself with the input formats through these examples and adjust them accordingly.

## Installation <a id="installation"></a>

To start using DeHNSSo, please download the files from this GitHub [Repository](https://github.com/SvenWesterbeek/DeHNSSo) **[#CMM: update this accordingly. possibly GitLab]**.  The current version of DeHNSSo was made using MATLAB R2022b. Support for earlier versions is not guaranteed. This file contains several folders, namely: Callers, Tools, Data Files, and Documentation. All folders must remain together in a directory of your choosing.

The base flow data for the fourth example case "CFI over a step in a swept-wing BL" is too large for this repository and can be found in the release snapshot 4TU [Repository](https://PLACEHOLDER) adjusted for use in DeHNSSo. The raw files can be found in the [Repository](https://PLACEHOLDER) of J. Casacuberta instead. 

### Input files <a name="input-files"></a>

DeHNSSo requires several inputs in the form of MATLAB structure arrays (_struct_) that describe the basic state flow field, the numerical domain and grid, and an inflow boundary condition among others. This data should be provided in dimensionless form. The reference values should be provided in BF. The contents of these structure arrays should carry specific names. These can be summarized as follows:

1. _BF_ - Base flow domain, grid, velocity field, and reference values
2. _Grid_ - numerical domain and discretization
3. _Stab_ - Mode specifications, spectral truncation, and inflow data
4. _Opt_ - Solver options such as outflow buffer specifications and inflow amplitude rate of increase.

These structure arrays will be discussed in detail below to help you make your own input files.

### BF <a id="bf"></a>

_BF_ contains the base flow data, the grid on which it is defined, and the reference values for nondimensionalization. This data can be provided in a structured, unstructured grid, or vector format.

This data does not need to be presented on the same grid as presented in the _Grid_ struct as the data will be interpolated onto the numerical grid within DeHNSSo using the griddata function (method = 'cubic'). The numerical domain cannot exceed the domain of the base flow. In the table below, the required contents of _BF_ are described. In short, _BF.X_ and _BF.Y_ describe the locations where base flow quantities are described. _BF.U_, _BF.V_ and _BF.W_ describe the velocities in $x$, $y$ and $z$ respectively. Then, the _BF.dxU_, _BF.dxV_, _BF.dxW_, _BF.dyU_, _BF.dyV_, and _BF.dyW_ describe the streamwise and wall-normal derivatives of the aforementioned velocities. Lastly, _BF.lref_, _BF.Uref_, _BF.nu_ and _BF.Re_ are the reference values and corresponding Reynolds number defined as $Re = l_{ref}\times U_{ref}/\nu$.

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

The Grid structure contains the information on the numerical domain and grid used to generate both the Cartesian ($x$, $y$) and numerical body-fitted grid ($\xi$, $\eta$).

_Grid.wall_ presents DeHNSSo with the bottom wall coordinates via 2 rows of data. The first row contains the $x$-coordinates and the second row supplies the corresponding $y$-coordinates; This matrix can be of any size. However, it is preferred to be highly refined around any geometric wall features to ensure the interpolation is performed well. Note that sharp surfaces are not explicitly featured in this list. To include a step, define the step via the inputs _Grid.StepX_ and _Grid.StepH_. The step is then accounted for using this data via an embedded boundary method.

_Grid.mode_ allows for several built-in grid generations to be used. The options currently available are:

- "equidistant"
  - Creates a grid following equidistant streamwise discretization and an $\eta$ distribution of collocation points following Malik (1990) for a user-defined median collocation point $y\_i$.
- "refined"
  - Creates a grid with a streamwise refined grid based on a Gaussian distribution following user-defined inputs _Grid.mug, Grid.sig and Grid.ag._ Wall-normal distribution follows Malik (1990) for a user-defined median collocation point $y\_i$.
- "curved"
  - Creates an equidistant distribution of streamwise grid points over a curved surface with straight $\eta$ axes wall-normal to that surface. The $\eta$ axes are thus not parallel to the global $y$ axis. Wall-normal collocation points are clustered near the wall following Malik's (1990) mapping for a user-defined median collocation point $y\_i$.
- "wallorthogonal"
  - Creates a grid using elliptic generation, with orthogonality at the wall following exactly the $\eta$ distribution of the mapping of Malik (1990) according to a user-defined $y\_i$. The streamwise distribution is nearly equidistant but can be slightly adjusted for the sake of orthogonality.

The rest of the inputs to _Grid_ are shown in the table below:


| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| _Grid.nx_ | Number of streamwise stations | [-] | $1$ |
| _Grid.ny_ | Number of wall-normal stations | [-] | $1$ |
| _Grid.wall_ | Smooth wall description | [-] | ($nx\_{wall}$, $2$) |
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
| Grid.StepType | Sharp geometry type ("FFS") | [-] | String |

## <a id="stab"></a> Stab

The _Stab_ structure is used to define the mode ensemble of interest and present the solver with inflow conditions. The spectral truncation can be defined by _Stab.N_ and _Stab.M_ which are the maximum multiples of the fundamental _Stab.omega\_0_ and _Stab.beta\_0_ respectively. 

_Stab.IC_ sets the mode initialization method. Currently, two methods are implemented. "ILST" calls a routine that finds the solution to the local eigenvalue problem at the inflow for all modes that have a nonzero amplitude (presented in _Stab.A0_). Note that not all modes need to be supplied with an amplitude at the inflow. These results are then normalized with the maximum streamwise perturbation velocity and multiplied by the respective initialization amplitude. "ZERO" instead means no inflow condition is supplied. This generally means that the user intends to simulate the receptivity problem by supplying inhomogeneous boundary conditions. "LOAD" instead uses the perturbation profiles presented in _Stab.u0, Stab.v0, Stab.w0, Stab.p0_ to define the inflow boundary condition. The initial perturbation data is interpolated onto the numerical grid within the solver and can thus be supplied on any distribution of points consistent with _Stab.y0_.

_Stab.bcw_ is used to define inhomogeneous wall conditions in the streamwise, wall-normal, and spanwise velocity components per mode defined at the streamwise locations presented in _Stab.bcwx._

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| _Stab.N_ | Spectral truncation of beta modes | [-] | $1$ |
| _Stab.M_ | Spectral truncation of omega modes | [-] | $1$ |
| _Stab.A0_ | Initial amplitude of all modes | [-] | $((2N+1) \times (2M+1),1)$ |
| _Stab.omega\_0_ | Fundamental frequency | [-] | $1$ |
| _Stab.beta\_0_ | Fundamental spanwise wavelength | [-] | $1$ |
| _Stab.IC_ | Initialization method "ILST","ZERO'', 'LOAD" | [-] | string |
| _Stab.bcwx_ | Inhomogeneous boundary condition locations | [-] | $(any,1)$ |
| _Stab.bcw_ | Inhomogeneous boundary conditions (default = 0's) | [-] | $( any, 3\times(2N+1)\times(2M+1))$ |    **[#CMM: in the code you also have top boundary condition. add descriptions here]**
| _Stab.u0_ | Normalized streamwise perturbation velocity at x\_0 | [-] | $(3 \times (2N+1) \times (2M+1)),ny)$ | **[#CMM: make it more clear what normalised means]**
| _Stab.v0_ | Normalized wall-normal perturbation velocity at x\_0 | [-] | $(3 \times (2N+1) \times (2M+1)),ny)$ |
| _Stab.w0_ | Normalized spanwise perturbation velocity at x\_0 | [-] | $(3 \times (2N+1) \times (2M+1)),ny)$ |
| _Stab.p0_ | Normalized perturbation pressure at x\_0 | [-] | $(3 \times (2N+1) \times (2M+1)),ny)$ |
| _Stab.y0_ | Wall-normal distribution of inflow perturbation data | [-] | $(1,ny)$ |

## Opt <a id="opt"></a>

The input structure to DeHNSSo contains solver-specific options. All of these have default options which will work for most cases. The user can choose to overwrite these to improve solver convergence, speed and numerical behaviour for specific cases if necessary.

The buffer can be adjusted using the inputs for the starting location (_Opt.xb_) and strength (_Opt.kappa_) as well as the start of the buffer on nonlinear terms (_Opt.nltbufxb_). The expected amplitude that a mode needs to have can be adjusted via the option _Opt.Th_.  High inflow amplitudes will likely be damped to improve the odds of converging the nonlinear terms. The maximum amplitude that a mode can reach linearly in the first iteration can be adjusted via _Opt.AMAX_. The applied damping results in a lower inflow amplitude. This amplitude is increased every iteration by rate of  _Opt.AFg_. The results of intermediate steps can be saved by supplying the _Opt.Sweep_ parameter with a 1. Intermediate results need not be converged to the criterion supplied by _opt.Conv_. However, if _Opt.Sweep_ $=1$, the input _Opt.ConvF_ allows the user to specify the required convergence of intermediate spteps as a factor of the final convergence criterion (_Opt.Conv_)

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| Opt.xb | Buffer starting location [0 1] (default = 0.85)| [-] | $1$ |
| Opt.kappa | Buffer strength [1 -\>] (default = 6) | [-] | $1$ |
| Opt.nltbufxb | Nonlinear term buffer starting location (default = _Opt.xb_)| [-] | $1$ |
| Opt.Th | Nonlinear mode introduction threshold (default = $10^{-11}$)| [-] | $1$ |
| Opt.Sweep | Output intermediate results flag (true = $1$, false = $0$ (default)) | [-] | $1$ |
| Opt.AFg | Amplitude factor rate of increase (default = 1.1) | [-] | $1$ |
| Opt.Conv | Convergence criterion (default = 1e-4) | [-] | $1$ |
| Opt.ConvF | Convergence criterion relaxation during ramping (default = 100) | [-] | $1$ |
| Opt.AMAX | Maximum amplitude for initializing ramping procedure (default = 0.1) | [-] | $1$ |

## output files <a id="output-files"></a>

The solver returns the following structures with outputs:

1. StabGrid; The numerical grid generated within the solver as well as all grid transformations
2. StabRes; Stability calculation results
3. BF; Base flow values interpolated on the numerical grid

The output structs are discussed in more detail below.

## StabGrid <a id="stabgrid"></a>

The StabGrid structure contains the numerical grid generated in the solver on which the simulation results are defined. The streamwise location increases with the column index while the wall-normal location decreases with the row index. The StabGrid contains both the physical ($x$, $y$) and the computational ($\xi$, $\eta$) grid. The transformation coefficients are also presented in this structure.

| Name | Content | Unit | Size |
| --- | --- | --- | --- |
| StabGrid.x | Global streamwise coordinate | [-] | $(nx,ny)$ |
| StabGrid.y | Global wall-normal coordinate | [-] | $(nx,ny)$ |
| StabGrid.xi | Computational streamwise coordinate | [-] | $(nx,ny)$ |
| StabGrid.eta | Computational wall-normal coordinate | [-] | $(nx,ny)$ |
| StabGrid.xix | $\frac{\partial \xi}{\partial x}$ transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.xiy | $\frac{\partial \xi}{\partial y}$ transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.xixx | $\frac{\partial^2 \xi}{\partial x^2}$ transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.xiyy | $\frac{\partial^2 \xi}{\partial y^2}$ transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.etax | $\frac{\partial \eta}{\partial x}$ transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.etay | $\frac{\partial \eta}{\partial y}$ transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.etaxx | $\frac{\partial^2 \eta}{\partial x^2}$ transformation coefficient | [-] | $(nx,ny)$ |
| StabGrid.etayy | $\frac{\partial^2 \eta}{\partial y^2}$ transformation coefficient | [-] | $(nx,ny)$ |

## StabRes <a id="stabres"></a>

The StabRes structure contains all the stability results defined on the locations defined by StabGrid.x, StabGrid.y. Additionally, some key factors are calculated that are commonly used in stability analysis for comparison purposes even when they might not be used in DeHNSSo (such as the streamwise wavenumber $\alpha$). Results are shown nondimensionally and normalized. In other words, amplitudes are extracted from the perturbation shape functions of $\hat{u}$, $\hat{v}$, $\hat{w}$, and $\hat{p}$ based on the maximum of the absolute streamwise velocity perturbation value $u$.

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

The BF structure is both an input and output. In the output, the structure is appended with the result of the base flow interpolation indicated by the subscript $_r$ for easier post-processing. The following quantities are added:

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

The development of Tollmien-Schlichting instabilities in a Blasius boundary layer is considered in the first case. This case was previously considered by Bertolotti _et al._ (1992), Chang _et al._ (1993), Herbert(1993, Agard report 793, p.4-19, figure 19), and Herbert (1997).

### Reference values
All quantities are nondimensionalized by the reference velocity $U\_0=10$ m/s and $l\_{ref} = 6.075 \times 10^{-4}$ m. The kinematic viscosity is $1.518$ m^2/s such that the global $Re = U\_0 \times l\_{ref}/\nu = 400$ in accordance with Bertolotti _et al._ (1992), Chang _et al._ (1993), Herbert(1993, Agard report 793, p.4-19), and Herbert (1997).

### Domain description

The flow setup can be described by a constant external velocity of 10 m/s (and zero pressure gradient) over a domain (including the buffer) ranging from $x = 400$ to $x = 2854$. The outflow buffer is initiated at 85% of the domain. The domain height is set to $H = 99$. The flow is fully two-dimensional, thus spanwise velocity is zero.

### Base Flow

The base flow is the solution to the incompressible boundary layer equation found using an in-house solver that employs second-order backward streamwise discretization and wall-normal derivatives are calculated using a spectral collocation method with Chebyshev polynomial bases. The base flow simulation was performed on a fine equidistant numerical grid of 5000 streamwise stations by 100 collocation points and stored in the BF structure loaded here.

### Nonlinear mode ensemble and amplitude ramping

The inflow conditions are comprised of the solution to the local eigenvalue problem at the inflow for $90.6$ Hz corresponding to $\omega = 0.0344$ superimposed with an initial streamwise perturbation amplitude of $A = 0.00125 \sqrt{2}$. This inflow amplitude does not require any amplitude ramping.

For this example case, the spectral domain is truncated at $M = 5$ (and $N = 0$) equal to the number of modes presented in Bertolotti _et al._ (1992), Chang _et al._ (1993), Herbert(1993, Agard report 793, p.4-19), and Herbert (1997). Higher harmonics and the mean flow distortion are not presented at the inflow and rise naturally downstream of the inflow through nonlinear forcing.

Running this case exactly as given will result in the data used to plot Figure 8 from Westerbeek _et al._ (2023).

## Swept-wing boundary layer: Stationary Crossflow Instability (CFI) <a id="swept-wing-boundary-layer-stationary-crossflow-instability"></a>

In the second example, the stability of a swept-wing boundary layer is assessed nonlinearly for stationary crossflow instabilities. The boundary layer is simulated on a flat plate mimicking the experiments of Rius-Vidales _et al._ (2021). This is done by imposing a fitted external velocity distribution on the top boundary of the base flow simulation as presented in Casacuberta _et al._ (2022). The external velocity is given by the equation:

$U_e(x) = 0.0023 \ln(x)^4 + 0.0377 \ln(x)^3 + 0.1752 \ln(x)^2 + 0.5303 \ln(x) + 1.874.$

The spanwise velocity $W_e=-1.24$ is constant over the domain. This holds for all cases described hereafter.

### Reference values
The reference length, defined as the Blasius length at the inflow is $l\_{ref} = 2.1394 \times 10^{-4}$ m. The external velocity at the inflow is $U\_{ref} = 15.1$ m/s. The kinematic viscosity of the flow is $1.47 \times 10^{-5}$ such that the global Reynolds number is $220$.

### Domain description

The problem is described additionally by the streamwise coordinate ranging from $x = 220$ to $x = 1542$. The domain height is set to $H = 89$. The outflow buffer is initialized from 85% of the domain and an additional damping of the nonlinear terms starts from 80% of the domain.

### Base Flow

The base flow data is an interpolated version of DNS results for the problem as described in the previous paragraph using the solver INCA (see Hickel and Adams (2008) and Hickel _et al._ (2014)). This data was kindly provided by J. Casacuberta who performed a full DNS simulation for this problem in J. Casacuberta _et al._ (2022).

### Nonlinear mode ensemble and amplitude ramping

The spectral domain is truncated at five harmonics for this case ($N=5$, $M=0$) and the mean flow distortion. Only the fundamental mode, characterized by a spanwise wavelength of $\lambda_z = 7.5$ mm or $\beta = 0.18$, is introduced at the inflow as the solution to the local eigenvalue problem. The inflow amplitude of $A = 3.5 \times 10^{-2}$ is imposed on the result. For this case, amplitude ramping is required. The inflow amplitude is increased by 10% each iteration after an amplitude reduction ensures that the linear simulation result is capped at AMAX ($= 0.1$, the default value). Reference data for the stability solution is presented by J. Casacuberta (2021) as well as NPSE solutions for the current problem provided by the authors of the current code and also previously shown in the aforementioned work.

## Swept-wing boundary layer: Interaction of Stationary Crossflow Instability with a hump <a id="swept-wing-boundary-layer-interaction-of-a-stationary-cfi-with-a-hump"></a>

This third simulation considers the interaction of a stationary CFI with a smooth hump. This case was previously examined in both Westerbeek (2023a) and the article on the current solver Westerbeek (2023b).

### Reference values
Lengths are normalized by the reference length $l\_{ref} = 2.14 \times 10^{-4}$ m with is the Blasius length at the inflow. The external velocity at the inflow is $U\_{ref} = 15.1$ m/s. The kinematic viscosity of the flow is $1.47 \times 10^{-5}$. This leads to the same reference values as for the flat plate case described by a global Reynolds number of $220$.

### Domain description

The domain is described by $x = 220$ to $x = 2150$ and a domain height of $H=89$. 

The shallow hump is symmetrical and centered around $x\_m = 589$. The base flow features no flow separation from the presence of the hump. The flow problem is thus not too challenging and does not demand grid refinement around the hump. Still, a slight streamwise refinement is introduced (using _Grid.mode_ _"xrefined"_) around the hump for users to get familiar with it.

### Base Flow
For this problem, the base flow calculation was performed in COMSOL (COMSOL 2022). The same external flow velocity as for the flat plate case is prescribed at the top boundary of the base flow calculation.

### Nonlinear mode ensemble, reference data, and ramping

The stability of this flow problem is assessed linearly ($N=1$, $M=0$) for computational efficiency. Only the fundamental CFI, characterized by $\beta = 0.18$ is introduced at the inflow as the solution to the local eigenvalue problem. The reference data is provided by J. Franco of DLR using AHLNS (see Franco 2018) on the same base flow. The AHLNS is physically equivalent to the current HNS for linear simulations. No amplitude ramping is required as this concerns a linear simulation.

The default value of $x_b = 0.85$ is used to define a buffer region covering the last 15% of the domain.

## Swept-wing boundary layer: Interaction of Stationary Crossflow Instability with a Forward-Facing Step <a id="swept-wing-boundary-layer-interaction-of-a-stationary-cfi-with-a-step"></a>

This last simulation features a swept-wing boundary layer featuring a step. This flow problem was previously examined in Casacuberta (2022) (largest step) and is an adaptation on the experiments of Rius-Vidales and Kotsonis (2021). 

### Reference values
Lengths are normalized by the same reference length as the previous two cases: $l\_{ref} = 2.14 \times 10^{-4}$ m. The external velocity at the inflow is $U\_{ref} = 15.1$ m/s. The kinematic viscosity of the flow is $1.47 \times 10^{-5}$. This leads to the same reference values as for the flat plate case described by a global Reynolds number of $220$.

### Domain description

The rectangular domain is described by $x = 685$ to $x = 993$ and a height of $H=89$. The step is located at $x = 859$ and has a height of $3.5$. This step wall is not supplied as wall data. Instead, it is accounted for using an embedded boundary method over a flat wall. The default value of $x_b = 0.85$ is used to define a buffer region covering the last 15% of the domain.

### Base Flow
The base flow data for this case is an interpolated version of DNS results calculated using INCA (see Hickel and Adams (2008) and Hickel _et al._(2014)). This data was kindly provided by J. Casacuberta who performed a full DNS simulation for this problem in J. Casacuberta _et al._ (2022). This data was adjusted to fit DeHNSSo's input format. The raw DNS data can be found [here](https://PLACEHOLDER).

### Nonlinear mode ensemble, reference data, and ramping
The stability of this flow problem is assessed linearly ($N=1$, $M=0$). Given that the problem is assessed linearly, no ramping was necessary. Only the fundamental CFI, characterized by $\beta = 0.18$ is introduced at the inflow as the solution to the local eigenvalue problem. The reference data is full transitional DNS data calculated using INCA provided by J. Casacuberta. 


# License<a id="license"></a>

The contents in the Data files directory are licensed under a  **CC-BY 4.0**  (see CC-BY-4.0 file and [here](/https://joinup.ec.europa.eu/licence/creative-commons-attribution-40-international-cc-40)). The codes and any other file in this repository are licensed under a  **GPLv3 license**  (see GPLv3 file and [here](/https://www.gnu.org/licenses/gpl-3.0.en.html)).

Copyright notice:

Technische Universiteit Delft hereby disclaims all copyright interest in the program "DeHNSSo". DeHNSSo is a MATLAB tool to solve Nonlinear Harmonic Navier-Stokes problems written by the Author(s).

Henri Werij, Dean of Faculty of Aerospace Engineering, Technische Universiteit Delft.

© 2023, Sven Westerbeek & Marios Kotsonis

# Common errors and solutions <a id="common-errors-and-solutions"></a>
This section covers some of the most common errors and their solutions. If your specific question is not answered in this section, please send your questions to S.H.J.Westerbeek@tudelft.nl or M.Kotsonis@tudelft.nl and it will be patched if possible. However, some errors might not be a bug and instead, be the result of improper use. In this section, some hints will be posted as to what might have gone wrong when certain errors appear and how to fix them. Additionally, DeHNSSo might be adjusted to provide error messages with solution hints within MATLAB directly.

# FAQ<a id="faq"></a>
Questions that we receive concerning DeHNSSo will be used to improve DeHNSSo. This concerns either the code itself, the commenting and the documentation. We will post the questions with answers in this section.

# DeHNSSo results<a id="dehnsso-results"></a>
Here we will post results obtained using DeHNSSo. If you have worked with DeHNSSo and have figures that you would like to share, please contact us at:
S.H.J.Westerbeek@tudelft.nl
M.Kotsonis@tudelft.nl

# How to cite DeHNSSo<a id="how-to-cite-dehnsso"></a>
If you have used DeHNSSo or any subroutine within this GitHub. Please cite the GitHub and the journal article:

S. Westerbeek, S. Hulshoff, H. Schuttelaars & M. Kotsonis (2023) DeHNSSo: An Efficient Harmonic Navier-Stokes Solver for Non-Linear Stability Problems with Complex 2D Geometric Features. DOI:
DeHNSSo. Delft Harmonic Navier-Stokes Solver. GitHub Repository. DOI: 

# Acknowledgments<a id="acknowledgments"></a>
We would like to appreciate several direct and indirect contributors for providing comments, data or codes. We thank J. Casacuberta and S. Hickel for allowing the use and publication of both the base flow data and reference results. Similarly, we thank J.A. Franco for providing AHLNS results as a reference and allowing us to publish this data in this repository. We appreciate the comments of H. Schuttelaars and S. Hullshof in the making of DeHNSSo. We thank J.A.C. Weidemann and S.C. Reddy for making their MATLAB differentiation matrix suite available and allowing it to be published alongside DeHNSSo.

# References<a id="references"></a>

Malik, M. R. (1990). Numerical methods for hypersonic boundary layer stability. Journal of computational physics, 86(2), 376-413.

Bertolotti, F. P., Herbert, T., & Spalart, P. R. (1992). Linear and nonlinear stability of the Blasius boundary layer. Journal of fluid mechanics, 242, 441-474.

Chang, C. L., Malik, M. R., Erlebacher, G., & Hussaini, M. Y. (1993). Linear and nonlinear PSE for compressible boundary layers. INSTITUTE FOR COMPUTER APPLICATIONS IN SCIENCE AND ENGINEERING HAMPTON VA.

Herbert, T. (1997). Parabolized stability equations. Annual Review of Fluid Mechanics, 29(1), 245-283.

Weideman, J. A., & Reddy, S. C. (2000). A MATLAB differentiation matrix suite. ACM Transactions on Mathematical Software (TOMS), 26(4), 465-519.

Franco, J. A., Hein, S., & Valero, E. (2018). Effect of humps and indentations on boundary-layer transition of compressible flows using the AHLNS methodology. Proc. 6th ECCM–7th ECFD, Paper, 2018.

Rius-Vidales, A. F., & Kotsonis, M. (2021). Impact of a forward-facing step on the development of crossflow instability. Journal of Fluid Mechanics, 924, A34.

Casacuberta, J., Hickel, S., Westerbeek, S., & Kotsonis, M. (2022). Direct numerical simulation of interaction between a stationary crossflow instability and forward-facing steps. Journal of Fluid Mechanics, 943, A46.

Westerbeek, S., Franco Sumariva, J. A., Michelis, T., Hein, S., & Kotsonis, M. (2023). Linear and Nonlinear Stability Analysis of a Three-Dimensional Boundary Layer over a Hump. In AIAA SCITECH 2023 Forum (p. 0678).

Banner created using: https://manytools.org/hacker-tools/ascii-banner/  
