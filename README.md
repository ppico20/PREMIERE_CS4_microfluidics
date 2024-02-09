# AgNPs_PSD: Predicting the Particle Size Distribution of AgNPs using PBM-CFD Simulations

## Overview
This project focuses on providing a numerical framework based on coupled PBM-CFD simulations to predict the spatio-temporal evolution of the particle size distribution (PSD) of silver nanoparticles (AgNPs) in both microchannels and well-mixed batch reactors. The solvers and libraries are based on **[OpenQBMM](https://github.com/OpenQBMM/OpenQBMM)**, an OpenFoam-based implementation of the **[Quadrature Method of Moments](https://www.sciencedirect.com/science/article/pii/S0009250917306590?via%3Dihub)**. More information can be found in our publication in Chemical Engineering Journal:

**<ins>Pico, P.</ins>**, Nathanael, K., Lavino, A.D., Kovalchuk, N.M., Simmons, M.J.H. and Matar, O.K. (2023). "*Silver nanoparticles synthesis in microfluidic and well-mixed reactors: A combined experimental and PBM-CFD study*". Chem. Eng. J., 474, p.145692. **[DOI: 10.1016/j.cej.2023.145692](https://www.sciencedirect.com/science/article/pii/S1385894723044236)**

If you'd like to learn more about AgNPs in general, please have a look at our review paper:

Nathanael, K., **<ins>Pico, P.</ins>**, Kovalchuk, N.M., Lavino, A.D., Simmons, M.J.H. and Matar, O.K. (2022). "*Computational modelling and microfluidics as emerging approaches to synthesis of silver nanoparticles – A review*". Chem. Eng. J., 436, p.135178. **[DOI: 10.1016/j.cej.2022.135178](https://www.sciencedirect.com/science/article/pii/S1385894722006830)**

## Table of Contents
- [Background](#background)
- [Available solvers](#available-solvers)
- [Available libraries](#available-libraries)
- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)
- [Contact](#contact)

## Background

Our framework is based on an extension of the classic **[Finke-Watzky two-step mechanism](https://pubs.acs.org/doi/10.1021/ja9705102)** (FW), which predicts the temporal evolution of reagent consumption and PSD in batch reactors for reduction-based synthesis of metal nanoparticles. This model was originally designed for well-mixed systems in which effects related to transport phenomena in a reactor (i.e., reactants mixing) are not taken into account. Given the critical importance of mixing effects in microreactors, we propose a few modifications to the FW mechanism. A summary of these modifications is as follows:

- Addition of an elementary reduction reaction of the form $SN + \nu_{R}R \xrightarrow{k_{r}} Ag_{(l)}$, where $SN$ denotes the silver precursor, $R$ a generic reducing agent, $Ag_{(l)}$ silver atoms in liquid, and $k_{r}$ is the kinetic constant associated with this reduction reaction. 

- Addition of reactive convection-diffusion equations of the form $\frac{\partial (\rho y_{j})}{\partial t} + \nabla_{\textbf{x}}\cdot(\rho \textbf{u}y_{j}) = \nabla_{\textbf{x}}\cdot(D_{j}\nabla_{\textbf{x}}(\rho y_{j})) + S_{j}$ for species $i$. These equations determine the concentration of each species in the system, which are connected to PBM through models of nucleation and growth.
  
- Inclusion of size-dependent particle diffusion, convection, and agglomeration in the population balance equation.

The following figure shows a schematic of the proposed coupling between PBM and CFD. In essence, on the CFD side of things we solve for the hydrodynamic and reactive aspects of the system (i.e., velocity, pressure, and species concentration fields). On the PBM side, we solve for the PSD (using a univariate number density function) using a few nucleation, growth, and agglomeration models, which themselves depend on velocity, pressure, and species concentration:

![coupling](https://github.com/ppico20/PREMIERE_CS4_microfluidics/blob/master/Coupling_PBM-CFD.png)

## Available solvers

### [newReacting_buoyantPbePimpleFoam:](https://github.com/ppico20/PREMIERE_CS4_microfluidics/tree/master/applications/solvers/newReacting_buoyantPbePimpleFoam)

This is the package's primary solver. Besides solving the continuity and momentum conservation equations, it incorporates five additional reactive convection-diffusion equations for each of the main species present in our system ($SN$: silver precursor, $R$: reducing agent, $Ag_{(l)}$: silver atoms in liquid, $Ag_{s}$: silver nuclei, $Ag_{s2}$: silver nuclei which have grown into AgNPs). It also solves the population balance equation via the quadrature method of moments. In our formulation, we use a univariate number density function with particle size, L, as the internal coordinate. Therefore, we obtain the fiels of the first $k$ moments of the PSD in a specific domain; we can then approximate the PSD with its moments.

### [newReacting_buoyantSimpleFoam:](https://github.com/ppico20/PREMIERE_CS4_microfluidics/tree/master/applications/solvers/newReacting_buoyantSimpleFoam)

This is a testing solver to exclusively solve for velocity, pressure, and species concentration, without activating the PBM. Since applications related to nanoparticles synthesis involve highly diluted systems, the motion of the particles will not have a major influence on the continuous phase. A such, removing the PBM portion will not have an influence on the velocity field.

### [newReacting_pbeFoam:](https://github.com/ppico20/PREMIERE_CS4_microfluidics/tree/master/applications/solvers/reacting_pbeFoam)

This is another testing solver in which all equations are solved in a 'one-cell' domain. What this means is that all spatial terms are removed from the equations, leading to temporal dependencies only. This is equivalent to assuming all reagents are in a state of perfect mixing initially and thus no reaction delays occur due to diffusive or convective mixing.

## Available libraries:

One library, named 'libcompressible', was written to include custom models for nucleation, growth, and a special treatment for particle diffusion in the population balance equation, which we explain below:

### Nucleation model ([nucleation_reaction](https://github.com/ppico20/PREMIERE_CS4_microfluidics/tree/master/libcompressible/nucleationModels/nucleation_reaction)):

This is a nucleation model extended from the F-W mechanism. It calculates the nucleation rate (i.e., rate of silver nuclei produced per unit volume) using the following expression:

$J = N_{av}\left(\frac{d_{m}}{d_{crit}}\right)k_{1}C_{Ag_{(l)}},$

where $N_{av}$ (Avogrado constant), $d_{m}$ (size of element atom), $d_{crit}$ (critical nucleus size), and $k_{1}$ (nucleation kinetic constant) are constants of the model that depend on the nature of the synthesised NPs. These constants are declared in the OpenQBMM dictionary ```constant/populationBalanceProperties```. $C_{Ag_{(l)}} = \rho y_{Ag_{(l)}} / MW_{Ag_{(l)}}$ is the molar concentration of $Ag_{(l)}$, directly dependent on its mass fraction in solution, $y_{Ag_{(l)}}$, obtained from solving the aforementioned convection-diffusion equations. $\rho$ corresponds to the density of the continuous phase, and $MW_{Ag_{(l)}}$ the molecular weight of $Ag_{(l)}$, defined by the user in ```constant/reactionProperties```.

### Growth model ([growth_reaction](https://github.com/ppico20/PREMIERE_CS4_microfluidics/tree/master/libcompressible/growthModels/growth_reaction)):

This is a growth model extended from the F-W mechanism. It calculates the size-dependent growth rate using the following expression: 

$G = \frac{1}{3}k_{2}LC_{Ag_{(l)}}),$

where $k_{2}$ is the growth kinetic constant, also declared by the user in ```constant/populationBalanceProperties```.

### Diffusion model (stokesEinstein):

In solid-liquid systems, 

## Installation


## Usage


## Credits

This project is a collaborative effort between Imperial College London (ICL, numerical and computational part) and University of Birmingham (UB, experimental part). The project's contributors are the following:

- Paula Pico (ICL)
- Konstantina Nathanael (UB)
- Dr. Alessio Lavino (ICL)
- Dr. Nina Kovalchuk (UB)
- Prof. Mark Simmons (UB)
- Prof. Omar Matar (ICL)

## Contact
- p.pico20@imperial.ac.uk - Paula Pico
