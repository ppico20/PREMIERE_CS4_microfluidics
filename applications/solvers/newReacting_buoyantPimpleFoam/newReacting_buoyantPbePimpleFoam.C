/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
-------------------------------------------------------------------------------
2015-06-21 Alberto Passalacqua: Derived solver from buoyantPimpleFoam and
                                added solution of population balance.
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
    NewReacting_buoyantPbePimpleFoam

Description
    Transient solver for buoyant, turbulent flow of compressible fluids with
    solution of population balance.

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
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for with population balance with buoyancy and "
        "heat transfer."
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

        volScalarField& y_AgNO3 = mesh.lookupObjectRef<volScalarField>("y_AgNO3");
        volScalarField& y_reduc = mesh.lookupObjectRef<volScalarField>("y_reduc");
        volScalarField& y_Agl = mesh.lookupObjectRef<volScalarField>("y_Agl");
        volScalarField& y_Ags = mesh.lookupObjectRef<volScalarField>("y_Ags");
        volScalarField& y_H2O = mesh.lookupObjectRef<volScalarField>("y_H2O");

        #include "ySuSp.H"

        fvScalarMatrix y_AgNO3Eqn
        (
            fvm::ddt(rho,y_AgNO3)
            + fvm::div(phi, y_AgNO3)
            - fvm::laplacian(rho*D_AgNO3, y_AgNO3) == Su_AgNO3 + fvm::Sp(Sp_AgNO3, y_AgNO3)
        );

        fvScalarMatrix y_reducEqn
        (
            fvm::ddt(rho,y_reduc)
            + fvm::div(phi, y_reduc)
            - fvm::laplacian(rho*D_reduc, y_reduc) == Su_reduc + fvm::Sp(Sp_reduc, y_reduc)
        );

        fvScalarMatrix y_AglEqn
        (
            fvm::ddt(rho,y_Agl)
            + fvm::div(phi, y_Agl)
            - fvm::laplacian(rho*D_Agl, y_Agl) == Su_Agl + fvm::Sp(Sp_Agl, y_Agl)
        );

        fvScalarMatrix y_AgsEqn
        (
            fvm::ddt(rho,y_Ags)
            + fvm::div(phi, y_Ags)
            - fvm::laplacian(rho*D_Ags, y_Ags) == Su_Ags + fvm::Sp(Sp_Ags, y_Ags)
        );

        fvScalarMatrix y_Ags2Eqn
        (
            fvm::ddt(rho,y_Ags2)
            + fvm::div(phi, y_Ags2)
            - fvm::laplacian(rho*D_Ags2, y_Ags2) == Su_Ags2 + fvm::Sp(Sp_Ags2, y_Ags2)
        );

        rate_product_Agl_reduction = (mw_Agl)*(kr*(pow(rho*y_AgNO3/mw_AgNO3,nu_AgNO3))*(pow(rho*y_reduc/mw_reduc,nu_reduc)));
        rate_product_Ags_nuc = (mw_Ags)*(k1*(rho*y_Agl/mw_Agl));
        rate_product_Ags2_g = (mw_Ags2)*(k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags));

        y_AgNO3.relax();
        y_reduc.relax();
        y_Agl.relax();
        y_Ags.relax();
        y_Ags2.relax();
        y_H2O.relax();

        y_AgNO3Eqn.solve();
        y_reducEqn.solve();
        y_AglEqn.solve();
        y_AgsEqn.solve();
        y_Ags2Eqn.solve();

        y_H2O = scalar(1) - (y_AgNO3+y_reduc+y_Agl+y_Ags+y_Ags2);

        y_AgNO3.min(1.0);
        y_AgNO3.max(0.0);
        y_reduc.min(1.0);
        y_reduc.max(0.0);
        y_Agl.min(1.0);
        y_Agl.max(0.0);
        y_Ags.min(1.0);
        y_Ags.max(0.0);
        y_Ags2.min(1.0);
        y_Ags2.max(0.0);
        y_H2O.min(1.0);
        y_H2O.max(0.0);

        populationBalance->solve();

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
