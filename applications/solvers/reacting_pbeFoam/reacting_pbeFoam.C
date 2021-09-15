/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
-------------------------------------------------------------------------------
2015-06-21 Alberto Passalacqua: Derived solver from chemFoam.
2019-11-29 Alberto Passalacqua: Ported to OpenFOAM+ v1906.
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reacting_pbeFoam

Description
    Solver for population balance problems
    - designed for use on single cell cases to validate kernels
    - single cell mesh created on-the-fly
    - fields created on the fly from the initial conditions

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "OFstream.H"
#include "hexCellFvMesh.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Single-cell solver for the population balance equation."
    );

    argList::noParallel();

    #define CREATE_MESH createSingleCellMesh.H
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createSingleCellMesh.H"
    #include "createFields.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        rho = thermo.rho();

        const volScalarField& y_AgNO3 = mesh.lookupObject<volScalarField>("y_AgNO3");
        const volScalarField& y_tannic = mesh.lookupObject<volScalarField>("y_tannic");
        const volScalarField& y_tri = mesh.lookupObject<volScalarField>("y_tri");
        const volScalarField& y_H2O = mesh.lookupObject<volScalarField>("y_H2O");
        const volScalarField& y_Agl = mesh.lookupObject<volScalarField>("y_Agl");
        const volScalarField& y_Ags = mesh.lookupObject<volScalarField>("y_Ags");
      
        fvScalarMatrix y_AgNO3Eqn
        (
            fvm::ddt(rho,y_AgNO3) == - (mw_AgNO3)*kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_tannic/mw_tannic)
        );

        fvScalarMatrix y_tannicEqn
        (
            fvm::ddt(rho,y_tannic) == - (mw_tannic)*kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_tannic/mw_tannic)
        );

        fvScalarMatrix y_triEqn
        (
            fvm::ddt(rho,y_tri) 
        );

        fvScalarMatrix y_H2OEqn
        (
            fvm::ddt(rho,y_H2O)
        );

        fvScalarMatrix y_AglEqn
        (
            fvm::ddt(rho,y_Agl) == (mw_Agl)*(kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_tannic/mw_tannic)-k1*(rho*y_Agl/mw_Agl)-k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
        );

        fvScalarMatrix y_AgsEqn
        (
            fvm::ddt(rho,y_Ags) == (mw_Ags)*(k1*(rho*y_Agl/mw_Agl)-k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
        );

        fvScalarMatrix y_Ags2Eqn
        (
            fvm::ddt(rho,y_Ags2) == (mw_Ags2)*(k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
        );

        rate_product_Agl_reduction = (mw_Agl)*(kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_tannic/mw_tannic));
        rate_product_Ags_nuc = (mw_Ags)*(k1*(rho*y_Agl/mw_Agl));
        rate_product_Ags2_g = (mw_Ags2)*(k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags));

        y_AgNO3Eqn.solve();
        y_tannicEqn.solve();
        y_triEqn.solve();
        y_H2OEqn.solve();
        y_AglEqn.solve();
        y_AgsEqn.solve();
        y_Ags2Eqn.solve();

        turbulence->validate();
        populationBalance->solve();

        #include "output.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "Number of steps = " << runTime.timeIndex() << endl;
    Info << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
