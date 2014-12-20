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

#include "HrenyaSinclairConductivity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HrenyaSinclairConductivity, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        HrenyaSinclairConductivity,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::HrenyaSinclairConductivity::HrenyaSinclairConductivity(const Foam::dictionary& dict)
:
    conductivityModel(dict),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    L_(coeffsDict_.lookup("L"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HrenyaSinclairConductivity::~HrenyaSinclairConductivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField Foam::HrenyaSinclairConductivity::kappa
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
    volScalarField lamda = 
        scalar(1) + da/(6.0*sqrt(2.0)*(alpha+ scalar(1.0e-5)))/L_;

    return
        2.0*pow(alpha, 2.0)*rhoa*da*g0*(1.0+e)*ThetaSqrt/piSqrt
      + (9.0/8.0)*ThetaSqrt*piSqrt*rhoa*da*0.25*pow(1.0+e, 2.0)*(2.0*e-1.0)*pow(alpha, 2.0)/(49.0/16.0-33.0*e/16.0)
      + (15.0/16.0)*ThetaSqrt*piSqrt*alpha*rhoa*da*(0.5*e*e+0.25*e-0.75+lamda)/((49.0/16.0-33.0*e/16.0)*lamda)
      + (25.0/64.0)*ThetaSqrt*piSqrt*rhoa*da/((1.0+e)*(49.0/16.0-33.0*e/16.0)*lamda*g0);
}


// ************************************************************************* //
