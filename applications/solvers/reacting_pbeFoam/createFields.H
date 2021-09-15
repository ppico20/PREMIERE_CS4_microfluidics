    if (mesh.nCells() != 1)
    {
        FatalErrorIn(args.executable())
            << "Solver only applicable to single cell cases"
            << exit(FatalError);
    }

    Info<< nl << "Reading thermophysical properties\n" << endl;

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero),
        p.boundaryField().types()
    );

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

Info<< "Reading field y_tannic\n" <<endl;
volScalarField y_tannic
(
    IOobject
    (
         "y_tannic",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Reading field y_tri\n" <<endl;
volScalarField y_tri
(
    IOobject
    (
         "y_tri",
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


    #include "createPhi.H"

    Info << "Creating turbulence model.\n" << endl;
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

    autoPtr<populationBalanceModel> populationBalance
    (
        populationBalanceModel::New
        (
            "populationBalance",
            populationBalanceProperties,
            phi
        )
    );
