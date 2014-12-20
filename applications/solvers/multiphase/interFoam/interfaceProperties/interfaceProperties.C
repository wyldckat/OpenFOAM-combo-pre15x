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

Application
    interfaceProperties

Description
    Properties to aid interFoam :
    1. Correct the gamma boundary condition for dynamic contact angle.
    2. Calculate interface curvature.

\*---------------------------------------------------------------------------*/

#include "interfaceProperties.H"
#include "gammaContactAngleFvPatchScalarField.H"
#include "physicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const scalar interfaceProperties::convertToRad = physicalConstant::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void interfaceProperties::correctContactAngle
(
    fvPatchVectorFieldField& nHatb
) const
{
    const fvMesh& mesh = gamma_.mesh();
    const fvPatchScalarFieldField& gbf = gamma_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (typeid(gbf[patchi]) == typeid(gammaContactAngleFvPatchScalarField))
        {
            const gammaContactAngleFvPatchScalarField& gcap = 
                refCast<const gammaContactAngleFvPatchScalarField>(gbf[patchi]);

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch =
                mesh.Sf().boundaryField()[patchi]
               /mesh.magSf().boundaryField()[patchi];

            scalar theta0 = convertToRad*gcap.theta0();
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = gcap.uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > SMALL)
            {
                scalar thetaA = convertToRad*gcap.thetaA();
                scalar thetaR = convertToRad*gcap.thetaR();

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall =
                    U_.boundaryField()[patchi].patchInternalField()
                  - U_.boundaryField()[patchi];
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall =
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch;

                // Normalise nWall
                nWall /= (mag(nWall) + SMALL);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall = nWall & Uwall;

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            scalarField a12 = nHatPatch & AfHatPatch;

            scalarField b1 = cos(theta);

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det = 1.0 - a12*a12;

            scalarField a = (b1 - a12*b2)/det;
            scalarField b = (b2 - a12*b1)/det;

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}


void interfaceProperties::calculateK()
{
    const fvMesh& mesh = gamma_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of gamma
    volVectorField gradGamma = fvc::grad(gamma_);

    // Interpolated face-gradient of gamma
    surfaceVectorField gradGammaf = fvc::interpolate(gradGamma);
    //gradGammaf -=
    //    2*(mesh.Sf()/mesh.magSf())
    //    *(fvc::snGrad(gamma_) - (mesh.Sf() & gradGammaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv = gradGammaf/(mag(gradGammaf) + deltaN_);
    correctContactAngle(nHatfv.boundaryField());

    /*
    dimensionedVector cen
    (
        "cen",
        dimLength,
        sum(mesh.V()*gamma_*mesh.C())/sum(mesh.V()*gamma_)
    );
    surfaceVectorField nHatfv = -(mesh.Cf() - cen)/mag(mesh.Cf() - cen);
    */

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    //volVectorField nHat = gradGamma/(mag(gradGamma) + deltaN_);
    //nHat.boundaryField() = nHatfv.boundaryField();
    //K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interfaceProperties::interfaceProperties
(
    const volScalarField& gamma,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cGamma_
    (
        readScalar
        (
            gamma.mesh().solutionDict().subDict("PISO").lookup("cGamma")
        )
    ),
    sigma_(dict.lookup("sigma")),

    deltaN_
    (
        "deltaN",
        dimensionSet(0, -1, 0, 0, 0),
        1e-8/pow(average(gamma.mesh().V()), 1.0/3.0)
    ),

    gamma_(gamma),
    U_(U),

    nHatf_
    (
        (
            fvc::interpolate(fvc::grad(gamma_))
           /(mag(fvc::interpolate(fvc::grad(gamma_))) + deltaN_)
        ) & gamma_.mesh().Sf()
    ),

    K_
    (
        IOobject
        (
            "K",
            gamma_.time().timeName(),
            gamma_.mesh()
        ),
        -fvc::div(nHatf_)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
