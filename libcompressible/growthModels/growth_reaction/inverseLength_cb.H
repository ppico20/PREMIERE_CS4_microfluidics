/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::growthModels::growth_reaction

Description
    Linear evaporation model.

    \f[
        C_g = -E_v \frac{\pi d_i^3}{6}
    \f]

    where

    \vartable
        E_v         |  Evaporation coefficient [m^2/s]
        d_i         |  Diameter of droplet i [m^3]
    \endvartable

SourceFiles
    growth_reaction.C

\*---------------------------------------------------------------------------*/

#ifndef growth_reaction_H
#define growth_reaction_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class linear Declaration
\*---------------------------------------------------------------------------*/

class growth_reaction
:
    public growthModel
{

private:

    // Private data

        //- Continuous phase name
        const word continuousPhase_;

        //- [mol/m3]: Equilibrium concentration in the liquid-gas-hydrates interface.
        dimensionedScalar c_eq_;

public:

        //- Runtime type information
        TypeName("growth_reaction");


    // Constructors

        //- Construct from components
        growth_reaction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~growth_reaction();


    // Member Functions

        //- Aggregation kernel
        virtual scalar Kg
        (
            const scalar& d,
            const bool lengthBased,
            const label environment = 0
        ) const;

        //- Return phase space convection
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature
        );

        //- Return phase space convection
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
