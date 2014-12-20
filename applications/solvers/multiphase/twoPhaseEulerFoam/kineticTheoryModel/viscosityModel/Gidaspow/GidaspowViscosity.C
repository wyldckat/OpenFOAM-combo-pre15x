/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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

#include "GidaspowViscosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GidaspowViscosity, 0);

    addToRunTimeSelectionTable(viscosityModel, GidaspowViscosity, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::GidaspowViscosity::GidaspowViscosity(const Foam::dictionary& dict)
:
    viscosityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GidaspowViscosity::~GidaspowViscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField Foam::GidaspowViscosity::mua
(
    const Foam::volScalarField& alpha,
    const Foam::volScalarField& Theta,
    const Foam::volScalarField& g0,
    const Foam::dimensionedScalar& rhoa,
    const Foam::dimensionedScalar& da,
    const Foam::dimensionedScalar& e
) const
{
    const scalar piSqrt = pow(M_PI, 0.5);
    volScalarField ThetaSqrt = pow(Theta, 0.5);

    return
        (4.0/5.0)*pow(alpha, 2.0)*rhoa*da*g0*(1.0+e)*ThetaSqrt/piSqrt
      + (1.0/15.0)*ThetaSqrt*piSqrt*rhoa*da*g0*(1.0+e)*pow(alpha, 2.0)
      + (1.0/6.0)*ThetaSqrt*piSqrt*rhoa*da*alpha
      + (10.0/96.0)*ThetaSqrt*piSqrt*rhoa*da/((1.0+e)*g0);
}


// ************************************************************************* //
