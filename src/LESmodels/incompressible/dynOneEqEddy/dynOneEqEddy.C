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

#include "dynOneEqEddy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESmodels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynOneEqEddy, 0);
addToRunTimeSelectionTable(LESmodel, dynOneEqEddy, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

dimensionedScalar dynOneEqEddy::ck(const volTensorField& D) const
{
    volScalarField KK =
        0.5*(filter_(U() & U()) - (filter_(U()) & filter_(U())));

    volTensorField LL =
        dev(filter_(U() * U()) - (filter_(U()) * filter_(U())));

    volTensorField MM =
        delta()*(filter_(sqrt(k_)*D) - 2*sqrt(KK + filter_(k_))*filter_(D));

    dimensionedScalar MMMM = average(MM && MM);

    if (MMMM.value() > VSMALL)
    {
        return average(LL && MM)/MMMM;
    }
    else
    {
        return 0.0;
    }
}


dimensionedScalar dynOneEqEddy::ce(const volTensorField& D) const
{
    volScalarField KK =
        0.5*(filter_(U() & U()) - (filter_(U()) & filter_(U())));

    volScalarField mm =
        pow(KK + filter_(k_), 1.5)/(2*delta()) - filter_(pow(k_, 1.5))/delta();

    volScalarField ee =
        2*delta()*ck(D)*(filter_(sqrt(k_)*(D && D))
      - 2*sqrt(KK + filter_(k_))*(filter_(D) && filter_(D)));

    dimensionedScalar mmmm = average(mm && mm);

    if (mmmm.value() > VSMALL)
    {
        return average(ee*mm)/mmmm;
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynOneEqEddy::dynOneEqEddy
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESmodel(typeName, U, phi, transport),
    GenEddyVisc(U, phi, transport),

    filterPtr_(LESfilter::New(U.mesh(), LESmodelProperties())),
    filter_(filterPtr_())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynOneEqEddy::~dynOneEqEddy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynOneEqEddy::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);

    volTensorField D = symm(gradU);

    volScalarField P = 2.0*nuSgs_*magSqr(D);

    solve
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(DkEff(), k_)
    ==
       P
     - fvm::Sp(ce(D)*sqrt(k_)/delta(), k_)
    );

    bound(k_, k0());

    nuSgs_ = ck(D)*sqrt(k_)*delta();
    nuSgs_.correctBoundaryConditions();
}


bool dynOneEqEddy::read()
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
