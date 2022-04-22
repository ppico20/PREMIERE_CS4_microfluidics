/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
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

#include "stokesEinstein.H"
#include "addToRunTimeSelectionTable.H"

#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace diffusionModels
{
    defineTypeNameAndDebug(stokesEinstein, 0);

    addToRunTimeSelectionTable
    (
        diffusionModel,
        stokesEinstein,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::diffusionModels::stokesEinstein
::stokesEinstein
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    diffusionModel(dict, mesh),
    gammaLam_("gammaLam", dimViscosity, dict),
    D0_("D0", dimViscosity, dict),
    L0_("L0", dimLength, dict)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::diffusionModels::stokesEinstein
::~stokesEinstein()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::diffusionModels::stokesEinstein::diffusionSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)

{

    scalar dSource(0);

    const PtrList<volScalarNode>& nodes = quadrature.nodes();
    bool lengthBased = nodes[0].lengthBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return dSource;
    }

    label sizeOrder = momentOrder[sizeIndex];

    bool volumeFraction = nodes[0].useVolumeFraction();
    if (volumeFraction)
    {
        if (lengthBased)
        {
            sizeOrder += 3;
        }
        else
        {
            sizeOrder += 1;
        }
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volScalarNode& node = nodes[pNodeI];

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));
            scalar d = node.d(celli, bAbscissa);
            scalar n =
                node.n(celli, node.primaryWeight()[celli], bAbscissa);

            scalar dSourcei =

                (D0_.value()*L0_.value()*n*pow(bAbscissa, sizeOrder))/(bAbscissa + L0_.value());

               // n
               //*Kg(d, lengthBased)
               //*sizeOrder
               //*pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, nodei)
            {
                if (scalarIndexes[nodei] != sizeIndex)
                {
                    dSourcei *=
                        pow
                        (
                            node.primaryAbscissae()[nodei][celli],
                            momentOrder[scalarIndexes[nodei]]
                        );
                }
            }
            dSource += dSourcei;
        }

        return dSource;
    }

    forAll(nodes, pNodeI)
    {
        const volScalarNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights()[sizeIndex], sNodei)
        {
            scalar bAbscissa =
                max
                (
                    node.secondaryAbscissae()[sizeIndex][sNodei][celli],
                    scalar(0)
                );
            scalar d = node.d(celli, bAbscissa);
            scalar n =
                node.n(celli, node.primaryWeight()[celli], bAbscissa)
               *node.secondaryWeights()[sizeIndex][sNodei][celli];

            scalar dSourcei =

                (D0_.value()*L0_.value()*n*pow(bAbscissa, sizeOrder))/(bAbscissa + L0_.value());

                //n
               //*Kg(d, lengthBased)
               //*sizeOrder
               //*pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    dSourcei *=
                        node.secondaryWeights()[cmpt][sNodei][celli]
                       *pow
                        (
                            node.secondaryAbscissae()[cmpt][sNodei][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            dSource += dSourcei;
        }
    }

    return dSource;

}

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::diffusionModels::stokesEinstein
::momentDiff
(
    const volScalarField& moment
) const
{

    return fvm::laplacian(gammaLam_, moment);
}

// ************************************************************************* //
