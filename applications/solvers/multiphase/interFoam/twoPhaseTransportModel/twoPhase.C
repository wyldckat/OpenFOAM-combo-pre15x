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

Class
    twoPhase

\*---------------------------------------------------------------------------*/

#include "twoPhase.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoPhase, 0);

addToRunTimeSelectionTable
(
    transportModel, 
    twoPhase, 
    dictionary
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void twoPhase::calcNu()
{
    // Average kinematic viscosity calculated from averaging dynamic viscosity
    nu_ = (gamma_*rho1_*phase1_->nu() + (1.0 - gamma_)*rho2_*phase2_->nu())
        /(gamma_*rho1_ + (1.0 - gamma_)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhase::twoPhase
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& phaseName
)
:
    transportModel(U, phi, phaseName),

    phase1Name_(phaseTransportProperties_.lookup("phase1")),
    phase2Name_(phaseTransportProperties_.lookup("phase2")),

    phase1_(transportModel::New(U, phi, phase1Name_)),
    phase2_(transportModel::New(U, phi, phase2Name_)),

    rho1_(phase1_().phaseTransportProperties().lookup("rho")),
    rho2_(phase2_().phaseTransportProperties().lookup("rho")),

    gamma_(U_.db().lookupObject<const volScalarField> ("gamma")),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool twoPhase::read()
{
    return transportModel::read()
        && phase1_().read()
        && phase2_().read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
