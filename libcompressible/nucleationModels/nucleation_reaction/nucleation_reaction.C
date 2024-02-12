/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
    Modification made by: Paula Pico (Imperial College London)
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

\*---------------------------------------------------------------------------*/

#include "nucleation_reaction.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(nucleation_reaction, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        nucleation_reaction,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::nucleation_reaction::nucleation_reaction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),

    continuousPhase_(dict.lookupOrDefault("continuousPhase", word::null)),

    k1_("k1", pow(dimTime,-1), dict),
    Nav_("Nav", pow(dimMoles,-1.0), dict),
    dcrit_("dcrit", dimLength, dict),
    dm_("dm", dimLength, dict),
    mwAgs_("mwAgs", dimMass*pow(dimMoles,-1.0), dict)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::nucleation_reaction::~nucleation_reaction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::nucleation_reaction::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{


    const volScalarField& y_Agl = mesh_.lookupObject<volScalarField>("y_Agl");
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    scalar abscissaNucleation = dcrit_.value();

    //return ((mwAgs_.value()*k1_.value()*(y_Agl[celli]/mwAgs_.value()))/((3.1415/6.0)*pow3(dcrit_.value())))*pow(abscissaNucleation, momentOrder);

    return ((Nav_.value()*dm_.value())/(dcrit_.value()))*(k1_.value()*(rho[celli]*y_Agl[celli]/mwAgs_.value()))*pow(abscissaNucleation, momentOrder);

}

// ************************************************************************* //
