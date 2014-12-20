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
    BirdCarreau

Description


\*---------------------------------------------------------------------------*/

#include "BirdCarreau.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace transportModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(BirdCarreau, 0);

addToRunTimeSelectionTable
(
    transportModel,
    BirdCarreau,
    dictionary
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
tmp<volScalarField> BirdCarreau::calcNu() const
{
    return 
        nuInf_
      + (nu0_ - nuInf_)
       *pow(1.0 + sqr(k_*strainRate()), (n_ - 1.0)/2.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
BirdCarreau::BirdCarreau
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& phaseName
)
:
    transportModel(U, phi, phaseName),
    BirdCarreauCoeffs_
    (
        phaseTransportProperties_.subDict(typeName + "Coeffs")
    ),
    nu0_(BirdCarreauCoeffs_.lookup("nu0")),
    nuInf_(BirdCarreauCoeffs_.lookup("nuInf")),
    k_(BirdCarreauCoeffs_.lookup("k")),
    n_(BirdCarreauCoeffs_.lookup("n")),
    nu_
    (
        IOobject
        (
            "nu" + phaseName,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool BirdCarreau::read()
{
    if (transportModel::read())
    {
        BirdCarreauCoeffs_ = 
            phaseTransportProperties_.subDict(typeName + "Coeffs");

        BirdCarreauCoeffs_.lookup("nu0") >> nu0_;
        BirdCarreauCoeffs_.lookup("nuInf") >> nuInf_;
        BirdCarreauCoeffs_.lookup("k") >> k_;
        BirdCarreauCoeffs_.lookup("n") >> n_;

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace transportModels
} // End namespace Foam

// ************************************************************************* //
