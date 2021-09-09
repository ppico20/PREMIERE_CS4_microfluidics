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
    reacting_buoyantPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of compressible fluids for
    ventilation and heat-transfer.

    Turbulence is modelled using a run-time selectable compressible RAS or
    LES model.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow"
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
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        const volScalarField& y_AgNO3 = mesh.lookupObject<volScalarField>("y_AgNO3");
        const volScalarField& y_tannic = mesh.lookupObject<volScalarField>("y_tannic");
        const volScalarField& y_tri = mesh.lookupObject<volScalarField>("y_tri");
        const volScalarField& y_H2O = mesh.lookupObject<volScalarField>("y_H2O");
        const volScalarField& y_Agl = mesh.lookupObject<volScalarField>("y_Agl");
        const volScalarField& y_Ags = mesh.lookupObject<volScalarField>("y_Ags");
      
        fvScalarMatrix y_AgNO3Eqn
        (
            fvm::ddt(rho,y_AgNO3)
            + fvm::div(phi, y_AgNO3)
            - fvm::laplacian(rho*D_AgNO3, y_AgNO3) == - (mw_AgNO3)*kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_tannic/mw_tannic)
        );

        fvScalarMatrix y_tannicEqn
        (
            fvm::ddt(rho,y_tannic)
            + fvm::div(phi, y_tannic)
            - fvm::laplacian(rho*D_tannic, y_tannic) == - (mw_tannic)*kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_tannic/mw_tannic)
        );

        fvScalarMatrix y_triEqn
        (
            fvm::ddt(rho,y_tri)
            + fvm::div(phi, y_tri)
            - fvm::laplacian(rho*D_tri, y_tri)
        );

        fvScalarMatrix y_H2OEqn
        (
            fvm::ddt(rho,y_H2O)
            + fvm::div(phi, y_H2O)
            - fvm::laplacian(rho*D_H2O, y_H2O)
        );

        fvScalarMatrix y_AglEqn
        (
            fvm::ddt(rho,y_Agl)
            + fvm::div(phi, y_Agl)
            - fvm::laplacian(rho*D_Agl, y_Agl) == (mw_Agl)*(kr*(rho*y_AgNO3/mw_AgNO3)*(rho*y_tannic/mw_tannic)-k1*(rho*y_Agl/mw_Agl)-k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
        );

        fvScalarMatrix y_AgsEqn
        (
            fvm::ddt(rho,y_Ags)
            + fvm::div(phi, y_Ags)
            - fvm::laplacian(rho*D_Ags, y_Ags) == (mw_Ags)*(k1*(rho*y_Agl/mw_Agl)-k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
        );

        fvScalarMatrix y_Ags2Eqn
        (
            fvm::ddt(rho,y_Ags2)
            + fvm::div(phi, y_Ags2)
            - fvm::laplacian(rho*D_Ags2, y_Ags2) == (mw_Ags2)*(k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags))
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

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
