/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    buoyantSimpleFoam

Group
    grpHeatTransferSolvers

Description
    Steady-state solver for buoyant, turbulent flow of compressible fluids,
    including radiation, for ventilation and heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "EEqn.H"
            #include "pEqn.H"
        }

        turbulence->correct();

        volScalarField& y_AgNO3 = mesh.lookupObjectRef<volScalarField>("y_AgNO3");
        volScalarField& y_reduc = mesh.lookupObjectRef<volScalarField>("y_reduc");
        volScalarField& y_Agl = mesh.lookupObjectRef<volScalarField>("y_Agl");
        volScalarField& y_Ags = mesh.lookupObjectRef<volScalarField>("y_Ags");
        volScalarField& y_H2O = mesh.lookupObjectRef<volScalarField>("y_H2O");

        fvScalarMatrix y_AgNO3Eqn
        (
            fvm::div(phi, y_AgNO3)
            - fvm::laplacian(rho*D_AgNO3, y_AgNO3) == - (mw_AgNO3)*kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_reduc/mw_reduc)
        );

        fvScalarMatrix y_reducEqn
        (
            fvm::div(phi, y_reduc)
            - fvm::laplacian(rho*D_reduc, y_reduc) == - (mw_reduc)*kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_reduc/mw_reduc)
        );

        fvScalarMatrix y_AglEqn
        (
            fvm::div(phi, y_Agl)
            - fvm::laplacian(rho*D_Agl, y_Agl) == (mw_Agl)*(kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_reduc/mw_reduc)-k1*(rho*y_Agl/mw_Agl)-k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
        );

        fvScalarMatrix y_AgsEqn
        (
            fvm::div(phi, y_Ags)
            - fvm::laplacian(rho*D_Ags, y_Ags) == (mw_Ags)*(k1*(rho*y_Agl/mw_Agl)-k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
        );

        fvScalarMatrix y_Ags2Eqn
        (
            fvm::div(phi, y_Ags2)
            - fvm::laplacian(rho*D_Ags2, y_Ags2) == (mw_Ags2)*(k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
        );

        rate_product_Agl_reduction = (mw_Agl)*(kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_reduc/mw_reduc));
        rate_product_Ags_nuc = (mw_Ags)*(k1*(rho*y_Agl/mw_Agl));
        rate_product_Ags2_g = (mw_Ags2)*(k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags));

        y_AgNO3Eqn.solve();
        y_reducEqn.solve();
        y_AglEqn.solve();
        y_AgsEqn.solve();
        y_Ags2Eqn.solve();

        y_H2O = scalar(1) - (y_AgNO3+y_reduc+y_Agl+y_Ags+y_Ags2);


        runTime.write();

        runTime.printExecutionTime(Info);
    
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
