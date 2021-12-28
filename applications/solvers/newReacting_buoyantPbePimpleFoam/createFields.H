Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

mesh.setFluxRequired(p_rgh.name());

label pRefCell = 0;
scalar pRefValue = 0.0;

if (p_rgh.needReference())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

dimensionedScalar initialMass("initialMass", fvc::domainIntegrate(rho));

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"
#include "createRadiationModel.H"
#include "createFvOptions.H"


const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

Info<< "Reading reactionProperties\n" << endl;

IOdictionary reactionProperties
(
    IOobject
    (
        "reactionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar D_AgNO3
(
     reactionProperties.lookup("D_AgNO3")
);

dimensionedScalar D_reduc
(
     reactionProperties.lookup("D_reduc")
);

dimensionedScalar D_Agl
(
     reactionProperties.lookup("D_Agl")
);

dimensionedScalar D_Ags
(
     reactionProperties.lookup("D_Ags")
);

dimensionedScalar D_Ags2
(
     reactionProperties.lookup("D_Ags2")
);

dimensionedScalar D_H2O
(
     reactionProperties.lookup("D_H2O")
);

dimensionedScalar mw_AgNO3
(
     reactionProperties.lookup("mw_AgNO3")
);

dimensionedScalar mw_reduc
(
     reactionProperties.lookup("mw_reduc")
);

dimensionedScalar mw_Agl
(
     reactionProperties.lookup("mw_Agl")
);

dimensionedScalar mw_Ags
(
     reactionProperties.lookup("mw_Ags")
);

dimensionedScalar mw_Ags2
(
     reactionProperties.lookup("mw_Ags2")
);

dimensionedScalar mw_H2O
(
     reactionProperties.lookup("mw_H2O")
);

dimensionedScalar kr
(
     reactionProperties.lookup("kr")
);

dimensionedScalar k1
(
     reactionProperties.lookup("k1")
);

dimensionedScalar k2
(
     reactionProperties.lookup("k2")
);

Info<< "Reading field y_AgNO3\n" <<endl;
volScalarField y_AgNO3
(
    IOobject
    (
         "y_AgNO3",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field y_reduc\n" <<endl;
volScalarField y_reduc
(
    IOobject
    (
         "y_reduc",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field y_Agl\n" <<endl;
volScalarField y_Agl
(
    IOobject
    (
         "y_Agl",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field y_Ags\n" <<endl;
volScalarField y_Ags
(
    IOobject
    (
         "y_Ags",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field y_Ags2\n" <<endl;
volScalarField y_Ags2
(
    IOobject
    (
         "y_Ags2",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field y_H2O\n" <<endl;
volScalarField y_H2O
(
    IOobject
    (
         "y_H2O",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

    volScalarField rate_product_Agl_reduction
    (
        IOobject
        (
            "rate_product_Agl_reduction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/(dimVolume*dimTime)
     );

    volScalarField rate_product_Ags_nuc
    (
        IOobject
        (
            "rate_product_Ags_nuc",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/(dimVolume*dimTime)
     );

    volScalarField rate_product_Ags2_g
    (
        IOobject
        (
            "rate_product_Ags2_g",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/(dimVolume*dimTime)
     );

Info<< "Reading populationBalanceProperties\n" << endl;

IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

surfaceScalarField phiByRho("phiByRho", phi/fvc::interpolate(rho));

autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        "populationBalance", populationBalanceProperties, phiByRho
    )
);
