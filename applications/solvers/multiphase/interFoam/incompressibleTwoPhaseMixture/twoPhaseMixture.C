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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    twoPhaseMixture

\*---------------------------------------------------------------------------*/

#include "twoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void twoPhaseMixture::calcNu()
{
    volScalarField limitedGamma = min(max(gamma_, 0.0), 1.0);

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedGamma*rho1_ + (1.0 - limitedGamma)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhaseMixture::twoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    transportModel(U, phi),

    phase1Name_("phase1"),
    phase2Name_("phase2"),

    nuModel1_
    (
        viscosityModel::New
        (
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rho1_(nuModel1_->viscosityProperties().lookup("rho")),
    rho2_(nuModel2_->viscosityProperties().lookup("rho")),

    U_(U),
    phi_(phi),

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

tmp<volScalarField> twoPhaseMixture::mu() const
{
    volScalarField limitedGamma = min(max(gamma_, 0.0), 1.0);

    return tmp<volScalarField> 
    (
        new volScalarField
        (
            "mu",
            limitedGamma*rho1_*nuModel1_->nu()
          + (1.0 - limitedGamma)*rho2_*nuModel2_->nu()
        )
    );
}


tmp<surfaceScalarField> twoPhaseMixture::muf() const
{
    surfaceScalarField gammaf = min(max(fvc::interpolate(gamma_), 0.0), 1.0);

    return tmp<surfaceScalarField> 
    (
        new surfaceScalarField
        (
            "mu",
            gammaf*rho1_*fvc::interpolate(nuModel1_->nu())
          + (1.0 - gammaf)*rho2_*fvc::interpolate(nuModel2_->nu())
        )
    );
}


tmp<surfaceScalarField> twoPhaseMixture::nuf() const
{
    surfaceScalarField gammaf = min(max(fvc::interpolate(gamma_), 0.0), 1.0);

    return tmp<surfaceScalarField> 
    (
        new surfaceScalarField
        (
            "nu",
            (
                gammaf*rho1_*fvc::interpolate(nuModel1_->nu())
              + (1.0 - gammaf)*rho2_*fvc::interpolate(nuModel2_->nu())
            )/(gammaf*rho1_ + (1.0 - gammaf)*rho2_)
        )
    );
}


bool twoPhaseMixture::read()
{
    if (transportModel::read())
    {
        if (nuModel1_().read(*this) && nuModel2_().read(*this))
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
