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
\*---------------------------------------------------------------------------*/

#include "diffusion_driven.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(diffusion_driven, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        diffusion_driven,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::diffusion_driven
::diffusion_driven
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    growthModel(dict, mesh),

    continuousPhase_(dict.lookupOrDefault("continuousPhase", word::null)),
    D_Agl_("D_Agl", pow(dimLength,2)*pow(dimTime,-1.0), dict),
    vm_("vm", pow(dimLength,3), dict),
    Nav_("Nav", pow(dimMoles,-1.0), dict),
    mwAgl_("mwAgl", dimMass*pow(dimMoles,-1.0), dict)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::diffusion_driven::
~diffusion_driven()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar
Foam::populationBalanceSubModels::growthModels::diffusion_driven::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    scalar gSource(0);

    const volScalarField& y_Agl = mesh_.lookupObject<volScalarField>("y_Agl");
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    const PtrList<volScalarNode>& nodes = quadrature.nodes();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return gSource;
    }
    label sizeOrder = momentOrder[sizeIndex];


    bool lengthBased = nodes[0].lengthBased();
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

    if (sizeOrder < 1)
    {
        return gSource;
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

            scalar gSourcei =
                n
               *(Kg(d, lengthBased)*(D_Agl_.value()*vm_.value()*Nav_.value()*(rho[celli]*y_Agl[celli]/mwAgl_.value())))
               *sizeOrder
               *pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, nodei)
            {
                if (scalarIndexes[nodei] != sizeIndex)
                {
                    gSourcei *=
                        pow
                        (
                            node.primaryAbscissae()[nodei][celli],
                            momentOrder[scalarIndexes[nodei]]
                        );
                }
            }
            gSource += gSourcei;
        }

        return gSource;
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

            scalar gSourcei =
                n
               *(Kg(d, lengthBased)*(D_Agl_.value()*vm_.value()*Nav_.value()*(rho[celli]*y_Agl[celli]/mwAgl_.value())))
               *sizeOrder
               *pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    gSourcei *=
                        node.secondaryWeights()[cmpt][sNodei][celli]
                       *pow
                        (
                            node.secondaryAbscissae()[cmpt][sNodei][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            gSource += gSourcei;
        }
    }

    return gSource;
}



Foam::scalar
Foam::populationBalanceSubModels::growthModels::diffusion_driven::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature
)
{
    scalar gSource(0);

    const PtrList<volVelocityNode>& nodes = quadrature.nodes();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return gSource;
    }
    label sizeOrder = momentOrder[sizeIndex];

    bool lengthBased = nodes[0].lengthBased();
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

    if (sizeOrder < 1)
    {
        return gSource;
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    forAll(nodes, pNodeI)
    {
        const volVelocityNode& node = nodes[pNodeI];

        scalar bAbscissa =
            max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));
        scalar d = node.d(celli, bAbscissa);
        scalar n =
            node.n(celli, node.primaryWeight()[celli], bAbscissa);

        scalar gSourcei =
            n*Kg(d, lengthBased)*sizeOrder*pow(bAbscissa, sizeOrder - 1);

        forAll(scalarIndexes, nodei)
        {
            if (scalarIndexes[nodei] != sizeIndex)
            {
                gSourcei *=
                    pow
                    (
                        node.primaryAbscissae()[nodei][celli],
                        momentOrder[scalarIndexes[nodei]]
                    );
            }
        }
        forAll(velocityIndexes, cmpt)
        {
            gSourcei *=
                pow
                (
                    node.velocityAbscissae()[celli][cmpt],
                    momentOrder[velocityIndexes[cmpt]]
                );
        }
        gSource += gSourcei;
    }

    return gSource;
}


Foam::scalar
Foam::populationBalanceSubModels::growthModels::diffusion_driven::Kg
(
    const scalar& d,
    const bool lengthBased,
    const label environment
) const
{

  if(d < SMALL){

    return 0.0;

  } else {

    return (pow(d,-1.0));
  }
    //return (Cg_.value()/max(d,dmin_.value()));
}

// ************************************************************************* //
