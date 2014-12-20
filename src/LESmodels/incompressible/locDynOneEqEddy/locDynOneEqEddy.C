/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "locDynOneEqEddy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESmodels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(locDynOneEqEddy, 0);
addToRunTimeSelectionTable(LESmodel, locDynOneEqEddy, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

volScalarField locDynOneEqEddy::ck
(
    const volTensorField& D,
    const volScalarField& KK
) const
{
    volTensorField LL =
        simpleFilter_(dev(filter_(U() * U()) - (filter_(U()) * filter_(U()))));

    volTensorField MM = simpleFilter_(-2.0*delta()*pow(KK, 0.5)*filter_(D));

    volScalarField ck = 
        simpleFilter_(0.5*(LL && MM))
       /(
            simpleFilter_(MM && MM)
          + dimensionedScalar("small", sqr(MM.dimensions()), VSMALL)
        );

    return 0.5*(mag(ck) + ck);
}


volScalarField locDynOneEqEddy::ce
(
    const volTensorField& D,
    const volScalarField& KK
) const
{
    volScalarField ce =
        simpleFilter_(nuEff()*(filter_(D && D) - (filter_(D) && filter_(D))))
       /simpleFilter_(pow(KK, 1.5)/(2.0*delta()));

    return 0.5*(mag(ce) + ce);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

locDynOneEqEddy::locDynOneEqEddy
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESmodel(typeName, U, phi, transport),
    GenEddyVisc(U, phi, transport),

    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), LESmodelProperties())),
    filter_(filterPtr_())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

locDynOneEqEddy::~locDynOneEqEddy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void locDynOneEqEddy::correct(const tmp<volTensorField>& gradU)
{
    LESmodel::correct(gradU);

    volTensorField D = symm(gradU);

    volScalarField KK = 
        0.5*(filter_(U() & U()) - (filter_(U()) & filter_(U())));
    KK.max(dimensionedScalar("small", KK.dimensions(), SMALL));

    volScalarField P = 2.0*nuSgs_*(D && D);

    solve
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(DkEff(), k_)
    ==
       P
     - fvm::Sp(ce(D, KK)*sqrt(k_)/delta(), k_)
    );

    bound(k_, k0());

    nuSgs_ = ck(D, KK)*sqrt(k_)*delta();
    nuSgs_.correctBoundaryConditions();
}


bool locDynOneEqEddy::read()
{
    if (GenEddyVisc::read())
    {
        filter_.read(LESmodelProperties());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESmodels
} // End namespace Foam

// ************************************************************************* //
