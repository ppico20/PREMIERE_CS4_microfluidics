/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/

#include "univariatePopulationBalance_SE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(univariatePopulationBalance_SE, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        univariatePopulationBalance_SE,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::univariatePopulationBalance_SE
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    univariatePDFTransportModel(name, dict, phi.mesh(), phi, "RPlus"),
    populationBalanceModel(name, dict, phi),
    odeType(phi.mesh(), dict),
    diffusiveTerm_(5),
    aggregation_(dict.lookupOrDefault("aggregation", false)),
    breakup_(dict.lookupOrDefault("breakup", false)),
    growth_(dict.lookupOrDefault("growth", false)),
    nucleation_(dict.lookupOrDefault("nucleation", false)),
    aggregationKernel_(),
    breakupKernel_(),
    growthModel_(),
    diffusionModel_
    (
        Foam::populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    nucleationModel_()
{
    if (aggregation_)
    {
        aggregationKernel_ =
            Foam::populationBalanceSubModels::aggregationKernel::New
            (
                dict.subDict("aggregationKernel"),
                phi_.mesh()
            );
    }

    if (breakup_)
    {
        breakupKernel_ =
            Foam::populationBalanceSubModels::breakupKernel::New
            (
                dict.subDict("breakupKernel"),
                phi_.mesh()
            );
    }

    if (growth_)
    {
        growthModel_ =
            Foam::populationBalanceSubModels::growthModel::New
            (
                dict.subDict("growthModel"),
                phi_.mesh()
            );
    }

    if (nucleation_)
    {
        nucleationModel_ =
            Foam::populationBalanceSubModels::nucleationModel::New
            (
                dict.subDict("nucleationModel"),
                phi_.mesh()
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::~univariatePopulationBalance_SE()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::implicitMomentSource
(
    const volScalarMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}


void
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::updateCellMomentSource(const label)
{}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
)
{
    scalar source(0);

    if (aggregation_)
    {
        source +=
            aggregationKernel_->aggregationSource
            (
                momentOrder,
                celli,
                quadrature,
                environment
            );
    }

    if (breakup_)
    {
        source +=
            breakupKernel_->breakupSource
            (
                momentOrder,
                celli,
                quadrature
            );
    }

    if (growth_)
    {
        source +=
            growthModel_->phaseSpaceConvection
            (
                momentOrder,
                celli,
                quadrature
            );
    }

    if (nucleation_)
    {
        source += nucleationModel_->nucleationSource(momentOrder[0], celli);
    }

    source += diffusiveTerm_[momentOrder][celli];

    return source;
}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::realizableCo() const
{
    return univariatePDFTransportModel::realizableCo();
}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::CoNum() const
{
    return 0.0;
}


bool
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::solveMomentSources() const
{
    if (aggregation_ || breakup_ || growth_ || nucleation_)
    {
        return odeType::solveSources_;
    }

    scalar dSource(0);
    scalar L0(1e-12);
    scalar D0(1e-12);

    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment

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

                (D0*L0*n*pow(bAbscissa, sizeOrder))/(bAbscissa + L0);

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

                (D0*L0*n*pow(bAbscissa, sizeOrder))/(bAbscissa + L0);

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

    diffusiveTerm_[momentOrder] = fvc::laplacian(D0*L0, dSource);

    return false;
}


bool
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::solveMomentOde() const
{
    return odeType::solveOde_;
}


void
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::explicitMomentSource()
{
    odeType::solve(quadrature_, 0);
}


void
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::solve()
{
    univariatePDFTransportModel::solve();
}


bool 
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance_SE
::readIfModified()
{
    odeType::read
    (
        populationBalanceProperties_.subDict(type() + "Coeffs")
    );
    
    return true;
}


// ************************************************************************* //
