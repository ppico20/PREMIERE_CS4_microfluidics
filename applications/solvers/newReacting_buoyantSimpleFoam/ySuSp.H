volScalarField::Internal Sp_AgNO3
(
    IOobject
    (
        "Sp_AgNO3",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Sp_AgNO3", rho.dimensions()*y_AgNO3.dimensions()/dimTime, 0)
);

volScalarField::Internal Su_AgNO3
(
    IOobject
    (
        "Su_AgNO3",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Su_AgNO3", rho.dimensions()*y_AgNO3.dimensions()/dimTime, 0)
);

volScalarField::Internal Sp_reduc
(
    IOobject
    (
        "Sp_reduc",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Sp_reduc", rho.dimensions()*y_reduc.dimensions()/dimTime, 0)
);

volScalarField::Internal Su_reduc
(
    IOobject
    (
        "Su_reduc",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Su_reduc", rho.dimensions()*y_reduc.dimensions()/dimTime, 0)
);

volScalarField::Internal Sp_Agl
(
    IOobject
    (
        "Sp_Agl",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Sp_Agl", rho.dimensions()*y_Agl.dimensions()/dimTime, 0)
);

volScalarField::Internal Su_Agl
(
    IOobject
    (
        "Su_Agl",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Su_Agl", rho.dimensions()*y_Agl.dimensions()/dimTime, 0)
);

volScalarField::Internal Sp_Ags
(
    IOobject
    (
        "Sp_Ags",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Sp_Ags", rho.dimensions()*y_Ags.dimensions()/dimTime, 0)
);

volScalarField::Internal Su_Ags
(
    IOobject
    (
        "Su_Ags",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Su_Ags", rho.dimensions()*y_Ags.dimensions()/dimTime, 0)
);

volScalarField::Internal Sp_Ags2
(
    IOobject
    (
        "Sp_Ags2",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Sp_Ags2", rho.dimensions()*y_Ags2.dimensions()/dimTime, 0)
);

volScalarField::Internal Su_Ags2
(
    IOobject
    (
        "Su_Ags2",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Su_Ags2", rho.dimensions()*y_Ags2.dimensions()/dimTime, 0)
);

const volScalarField w_AgNO3((mw_AgNO3)*(-nu_AgNO3)*kr*(pow(rho*y_AgNO3/mw_AgNO3,nu_AgNO3))*(pow(rho*y_reduc/mw_reduc,nu_reduc)));

forAll(w_AgNO3, celli)
{
    if (w_AgNO3[celli] > 0.0)
    {
        Sp_AgNO3[celli] -= w_AgNO3[celli]/max(1.0 - y_AgNO3[celli], 1e-9);
        Su_AgNO3[celli] += w_AgNO3[celli]/max(1.0 - y_AgNO3[celli], 1e-9);
    }
    else if (w_AgNO3[celli] < 0.0)
    {
        Sp_AgNO3[celli] += w_AgNO3[celli]/max(y_AgNO3[celli], 1e-9);
    }
}

const volScalarField w_reduc((mw_reduc)*(-nu_reduc)*kr*(pow(rho*y_AgNO3/mw_AgNO3,nu_AgNO3))*(pow(rho*y_reduc/mw_reduc,nu_reduc)));

forAll(w_reduc, celli)
{
    //Sp_reduc[celli] += nu_reduc*w_reduc[celli]/max(y_reduc[celli], 1e-12);
    //Su_reduc[celli] -= w_reduc[celli] -  nu_reduc*w_reduc[celli]/max(y_reduc[celli], 1e-12);

    if (w_reduc[celli] > 0.0)
    {
        Sp_reduc[celli] -= w_reduc[celli]/max(1.0 - y_reduc[celli], 1e-9);
        Su_reduc[celli] += w_reduc[celli]/max(1.0 - y_reduc[celli], 1e-9);
    }
    else if (w_reduc[celli] < 0.0)
    {
        Sp_reduc[celli] += w_reduc[celli]/max(y_reduc[celli], 1e-9);
    }

}

const volScalarField w_Agl((mw_Agl)*(kr*(pow(rho*y_AgNO3/mw_AgNO3,nu_AgNO3))*(pow(rho*y_reduc/mw_reduc,nu_reduc))-k1*(rho*y_Agl/mw_Agl)-k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags)));

forAll(w_Agl, celli)
{
    if (w_Agl[celli] > 0.0)
    {
        Sp_Agl[celli] -= w_Agl[celli]/max(1.0 - y_Agl[celli], 1e-9);
        Su_Agl[celli] += w_Agl[celli]/max(1.0 - y_Agl[celli], 1e-9);
    }
    else if (w_Agl[celli] < 0.0)
    {
        Sp_Agl[celli] += w_Agl[celli]/max(y_Agl[celli], 1e-9);
    }
}

const volScalarField w_Ags((mw_Ags)*(k1*(rho*y_Agl/mw_Agl)-k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags)));

forAll(w_Ags, celli)
{
    if (w_Ags[celli] > 0.0)
    {
        Sp_Ags[celli] -= w_Ags[celli]/max(1.0 - y_Ags[celli], 1e-9);
        Su_Ags[celli] += w_Ags[celli]/max(1.0 - y_Ags[celli], 1e-9);
    }
    else if (w_Ags[celli] < 0.0)
    {
        Sp_Ags[celli] += w_Ags[celli]/max(y_Ags[celli], 1e-9);
    }
}

const volScalarField w_Ags2((mw_Ags2)*(k2*(rho*y_Agl/mw_Agl)*(rho*y_Ags/mw_Ags)));

forAll(w_Ags2, celli)
{
    if (w_Ags2[celli] > 0.0)
    {
        Sp_Ags2[celli] -= w_Ags2[celli]/max(1.0 - y_Ags2[celli], 1e-9);
        Su_Ags2[celli] += w_Ags2[celli]/max(1.0 - y_Ags2[celli], 1e-9);
    }
    else if (w_Ags2[celli] < 0.0)
    {
        Sp_Ags2[celli] += w_Ags2[celli]/max(y_Ags2[celli], 1e-9);
    }
}
