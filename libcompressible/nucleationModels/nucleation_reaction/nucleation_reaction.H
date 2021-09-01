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

Class
    Foam::populationBalanceSubModels::nucleationModels::nucleation_reaction

Description
    nucleation_reaction nucleation model for soot formation.

SourceFiles
    nucleation_reaction.C

\*---------------------------------------------------------------------------*/

#ifndef nucleation_reaction_H
#define nucleation_reaction_H

#include "nucleationModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                        Class nucleation_reaction Declaration
\*---------------------------------------------------------------------------*/

class nucleation_reaction
:
    public nucleationModel
{

private:

    // Private data

        //- Continuous phase name
        const word continuousPhase_;

        //- Reaction constant
        dimensionedScalar k1_;

        //- dcrit
        dimensionedScalar dcrit_;

        dimensionedScalar mwAgs_;



public:

        //- Runtime type information
        TypeName("nucleation_reaction");


    // Constructors

        //- Construct from components
        nucleation_reaction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~nucleation_reaction();


    // Member Functions

        //- Nucleation model
        virtual scalar nucleationSource
        (
            const label& momentOrder,
            const label celli,
            const label environment = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
