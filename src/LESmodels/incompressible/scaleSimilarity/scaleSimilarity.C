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

#include "scaleSimilarity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESmodels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scaleSimilarity, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

scaleSimilarity::scaleSimilarity
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESmodel(typeName, U, phi, transport),
    filterPtr_(LESfilter::New(U.mesh(), LESmodelProperties())),
    filter_(filterPtr_())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scaleSimilarity::~scaleSimilarity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volTensorField> scaleSimilarity::B() const
{
    return
    (
        (filter_(U() * U()) - (filter_(U()) * filter_(U())))
    );
}


tmp<volScalarField> scaleSimilarity::k() const
{
    return
    (
        0.5*(filter_(U() & U()) - (filter_(U()) & filter_(U())))
    );
}


tmp<volScalarField> scaleSimilarity::epsilon() const
{
    volTensorField D = symm(fvc::grad(U()));

    return
    (
        ((filter_(U() * U()) - (filter_(U())*filter_(U()))) && D)
    );
}


tmp<fvVectorMatrix> scaleSimilarity::divB(volVectorField& U) const
{
    return fvm::Su(divB(), U);
}


tmp<volVectorField> scaleSimilarity::divB() const
{
    return
    (
        fvc::div
        (
            filter_(U())*filter_(U())
          + 0.333*I*tr(filter_(U() * U()) - filter_(U())*filter_(U()))
          - filter_(U() * U()),
            "div(B)"
        )
    );
}


void scaleSimilarity::correct(const tmp<volTensorField>&)
{}


bool scaleSimilarity::read()
{
    if (LESmodel::read())
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
