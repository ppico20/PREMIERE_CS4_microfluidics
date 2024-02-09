# AgNPs_PSD: Predicting the Particle Size Distribution of AgNPs using PBM-CFD Simulations

## Overview
This project focuses on providing a numerical framework based on coupled PBM-CFD simulations to predict the spatio-temporal evolution of the particle size distribution (PSD) of silver nanoparticles (AgNPs) in both microchannels and well-mixed batch reactors. The solvers and libraries are based on **[OpenQBMM](https://github.com/OpenQBMM/OpenQBMM)**, an OpenFoam-based implementation of the **[Quadrature Method of Moments](https://www.sciencedirect.com/science/article/pii/S0009250917306590?via%3Dihub)**. More information can be found in our publication in Chemical Engineering Journal:

**<ins>Pico, P.</ins>**, Nathanael, K., Lavino, A.D., Kovalchuk, N.M., Simmons, M.J.H. and Matar, O.K. (2023). "*Silver nanoparticles synthesis in microfluidic and well-mixed reactors: A combined experimental and PBM-CFD study*". Chem. Eng. J., 474, p.145692. **[DOI: 10.1016/j.cej.2023.145692](https://www.sciencedirect.com/science/article/pii/S1385894723044236)**

If you'd like to learn more about AgNPs in general, please have a look at our review paper:

Nathanael, K., **<ins>Pico, P.</ins>**, Kovalchuk, N.M., Lavino, A.D., Simmons, M.J.H. and Matar, O.K. (2022). "*Computational modelling and microfluidics as emerging approaches to synthesis of silver nanoparticles – A review*". Chem. Eng. J., 436, p.135178. **[DOI: 10.1016/j.cej.2022.135178](https://www.sciencedirect.com/science/article/pii/S1385894722006830)**

## Table of Contents
- [Background](#background)
- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)
- [Contact](#contact)

## Background

Our framework is based on an extension of the classic **[Finke-Watzky two-step mechanism](https://pubs.acs.org/doi/10.1021/ja9705102)** (FW), which predicts the temporal evolution of reagent consumption and PSD in batch reactors for reduction-based synthesis of metal nanoparticles. This model was originally designed for well-mixed systems in which effects related to transport phenomena in a reactor (i.e., reactants mixing) are not taken into account. Given the critical importance of mixing effects in microreactors, we propose a few modifications to the FW mechanism. A summary of these modifications is as follows:

- Addition of an elementary reduction reaction of the form $SN + \nu_{R}R \xrightarrow{k_{r}} Ag_{(l)}$, where $SN$ denotes the silver precursor, $R$ a generic reducing agent, $Ag_{(l)}$ silver atoms in liquid, and $k_{r}$ is the kinetic constant associated with this reduction reaction. 

- Addition of reactive convection-diffusion equations of the form $\frac{\partial (\rho y_{j})}{\partial t} + \nabla_{\textbf{x}}\cdot(\rho \textbf{u}y_{j}) = \nabla_{\textbf{x}}\cdot(D_{j}\nabla_{\textbf{x}}(\rho y_{j})) + S_{j}$ for species $i$. These equations determine the concentration of each species in the system, which are connected to PBM througb models of nucleation and growth.
  
- Inclusion of size-dependent particle diffusion, convection, and agglomeration in the population balance equation.

![coupling](https://github.com/ppico20/PREMIERE_CS4_microfluidics/blob/master/Coupling_PBM-CFD.png)

## Available solvers

### newReacting_buoyantPbePimpleFoam:

### newReacting_pbeFoam: 

### FW_pbeFoam: 
Validation solver based on the FW two-step mechanism


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
