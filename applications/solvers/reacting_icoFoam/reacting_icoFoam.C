/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    reacting_icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        const volScalarField& c_AgNO3 = mesh.lookupObject<volScalarField>("c_AgNO3");
        const volScalarField& c_tannic = mesh.lookupObject<volScalarField>("c_tannic");
        const volScalarField& c_tri = mesh.lookupObject<volScalarField>("c_tri");
        const volScalarField& c_H2O = mesh.lookupObject<volScalarField>("c_H2O");
        const volScalarField& c_Agl = mesh.lookupObject<volScalarField>("c_Agl");
        const volScalarField& c_Ags = mesh.lookupObject<volScalarField>("c_Ags");
      
        fvScalarMatrix c_AgNO3Eqn
        (
            fvm::ddt(c_AgNO3)
            + fvm::div(phi, c_AgNO3)
            - fvm::laplacian(D_AgNO3, c_AgNO3) == - (mw_AgNO3)*kr*(c_AgNO3/mw_AgNO3)*(c_tannic/mw_tannic)
        );

        fvScalarMatrix c_tannicEqn
        (
            fvm::ddt(c_tannic)
            + fvm::div(phi, c_tannic)
            - fvm::laplacian(D_tannic, c_tannic) == - (mw_tannic)*kr*(c_AgNO3/mw_AgNO3)*(c_tannic/mw_tannic)
        );

        fvScalarMatrix c_triEqn
        (
            fvm::ddt(c_tri)
            + fvm::div(phi, c_tri)
            - fvm::laplacian(D_tri, c_tri)
        );

        fvScalarMatrix c_H2OEqn
        (
            fvm::ddt(c_H2O)
            + fvm::div(phi, c_H2O)
            - fvm::laplacian(D_H2O, c_H2O)
        );

        fvScalarMatrix c_AglEqn
        (
            fvm::ddt(c_Agl)
            + fvm::div(phi, c_Agl)
            - fvm::laplacian(D_Agl, c_Agl) == (mw_Agl)*(kr*(c_AgNO3/mw_AgNO3)*(c_tannic/mw_tannic)-k1*(c_Agl/mw_Agl)-k2*(c_Agl/mw_Agl)*(c_Ags/mw_Ags))
        );

        fvScalarMatrix c_AgsEqn
        (
            fvm::ddt(c_Ags)
            + fvm::div(phi, c_Ags)
            - fvm::laplacian(D_Ags, c_Ags) == (mw_Ags)*(k1*(c_Agl/mw_Agl)-k2*(c_Agl/mw_Agl)*(c_Ags/mw_Ags))
        );

        fvScalarMatrix c_Ags2Eqn
        (
            fvm::ddt(c_Ags2)
            + fvm::div(phi, c_Ags2)
            - fvm::laplacian(D_Ags2, c_Ags2) == (mw_Ags2)*(k2*(c_Agl/mw_Agl)*(c_Ags/mw_Ags))
        );

       rate_produc_Agl_reduction = (mw_Agl)*(kr*(c_AgNO3/mw_AgNO3)*(c_tannic/mw_tannic));
       rate_product_Ags_nuc = (mw_Ags)*(k1*(c_Agl/mw_Agl));
       rate_product_Ags2_g = (mw_Ags2)*(k2*(c_Agl/mw_Agl)*(c_Ags/mw_Ags));

        c_AgNO3Eqn.solve();
        c_tannicEqn.solve();
        c_triEqn.solve();
        c_H2OEqn.solve();
        c_AglEqn.solve();
        c_AgsEqn.solve();
        c_Ags2Eqn.solve();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
