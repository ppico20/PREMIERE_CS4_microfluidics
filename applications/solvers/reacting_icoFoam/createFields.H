Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

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

dimensionedScalar D_tannic
(
     reactionProperties.lookup("D_tannic")
);

dimensionedScalar D_tri
(
     reactionProperties.lookup("D_tri")
);

dimensionedScalar D_H2O
(
     reactionProperties.lookup("D_H2O")
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

dimensionedScalar mw_AgNO3
(
     reactionProperties.lookup("mw_AgNO3")
);

dimensionedScalar mw_tannic
(
     reactionProperties.lookup("mw_tannic")
);

dimensionedScalar mw_tri
(
     reactionProperties.lookup("mw_tri")
);

dimensionedScalar mw_H2O
(
     reactionProperties.lookup("mw_H2O")
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


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


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

Info<< "Reading field c_AgNO3\n" <<endl;
volScalarField c_AgNO3
(
    IOobject
    (
         "c_AgNO3",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field c_tannic\n" <<endl;
volScalarField c_tannic
(
    IOobject
    (
         "c_tannic",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field c_tri\n" <<endl;
volScalarField c_tri
(
    IOobject
    (
         "c_tri",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field c_H2O\n" <<endl;
volScalarField c_H2O
(
    IOobject
    (
         "c_H2O",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field c_Agl\n" <<endl;
volScalarField c_Agl
(
    IOobject
    (
         "c_Agl",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field c_Ags\n" <<endl;
volScalarField c_Ags
(
    IOobject
    (
         "c_Ags",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field c_Ags2\n" <<endl;
volScalarField c_Ags2
(
    IOobject
    (
         "c_Ags2",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

    volScalarField rate_produc_Agl_reduction
    (
        IOobject
        (
            "rate_produc_Agl_reduction",
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


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
