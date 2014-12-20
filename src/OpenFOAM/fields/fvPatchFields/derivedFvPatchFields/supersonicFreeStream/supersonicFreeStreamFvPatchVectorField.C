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

Description
    Supersonic free stream condition.

    Supersonic outflow is vented according to ???

    Supersonic inflow is assumed to occur according to the Prandtl-Meyer
    expansion process.

    Subsonic outflow is zero-gradiented from inside the domain.

    N.B. This boundary condition is ill-posed if the free-stream flow is
         normal to the boundary.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "supersonicFreeStreamFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

supersonicFreeStreamFvPatchVectorField::supersonicFreeStreamFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF
)
:
    mixedFvPatchVectorField(p, iF),
    UInf_(vector::zero),
    pInf_(0),
    TInf_(0)
{
    refValue() = patchInternalField();
    refGrad() = vector::zero;
    valueFraction() = 1;
}


supersonicFreeStreamFvPatchVectorField::supersonicFreeStreamFvPatchVectorField
(
    const supersonicFreeStreamFvPatchVectorField& ptf,
    const fvPatch& p,
    const vectorField& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    UInf_(ptf.UInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_)
{}


supersonicFreeStreamFvPatchVectorField::supersonicFreeStreamFvPatchVectorField
(
    const fvPatch& p,
    const vectorField& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    UInf_(dict.lookup("UInf")),
    pInf_(readScalar(dict.lookup("pInf"))),
    TInf_(readScalar(dict.lookup("TInf")))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 1;

    if (pInf_ < SMALL)
    {
        FatalIOErrorIn
        (
            "supersonicFreeStreamFvPatchVectorField::"
            "supersonicFreeStreamFvPatchVectorField"
            "(const fvPatch&, const vectorField&, const dictionary&)",
            dict
        )   << "    unphysical pInf specified (pInf <= 0.0)"
            << exit(FatalIOError);
    }
}


supersonicFreeStreamFvPatchVectorField::supersonicFreeStreamFvPatchVectorField
(
    const supersonicFreeStreamFvPatchVectorField& sfspvf,
    const vectorField& iF
)
:
    mixedFvPatchVectorField(sfspvf, iF),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate the patch field
void supersonicFreeStreamFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& pT =
        lookupPatchField<volScalarField, scalar>("T");

    const fvPatchField<scalar>& pp =
        lookupPatchField<volScalarField, scalar>("p");

    const fvPatchField<scalar>& ppsi =
        lookupPatchField<volScalarField, scalar>("psi");

    // Need R of the free-stream flow.  Assume R is independent of location
    // along patch so use face 0
    scalar R = 1.0/(ppsi[0]*pT[0]);


    scalar gamma = 1.296;


    scalar MachInf = mag(UInf_)/sqrt(gamma*R*TInf_);

    if (MachInf < 1.0)
    {
        FatalErrorIn
        (
            "supersonicFreeStreamFvPatchVectorField::updateCoeffs()"
        )   << "    MachInf < 1.0, free stream must be supersonic"
            << abort(FatalError);
    }

    vectorField& Up = refValue();
    valueFraction() = 1;

    // get the near patch internal cell values
    vectorField U = patchInternalField();


    // Find the component of U normal to the free-stream flow and in the
    // plane of the free-stream flow and the patch normal

    // Direction of the free-stream flow
    vector UInfHat = UInf_/mag(UInf_);

    // Normal to the plane defined by the free-stream and the patch normal
    vectorField nnInfHat = UInfHat ^ patchMesh().nf();

    //vectorField UnInf = U - nnInfHat*(nnInfHat & U);
    //vectorField Un = UnInf - UInfHat*(UInfHat & UnInf);
    //vectorField nHatInf =
    //    (Un/(mag(Un) + SMALL)) * sign(patchMesh().nf() & Un);

    // Normal to the free-stream in the plane defined by the free-stream
    // and the patch normal
    vectorField nHatInf = nnInfHat ^ UInfHat;

    // Component of U normal to the free-stream in the plane defined by the
    // free-stream and the patch normal
    vectorField Un = nHatInf*(nHatInf & U);

    // The tangential component is
    vectorField Ut = U - Un;

    // Calculate the Prandtl-Meyer function of the free-stream
    scalar nuMachInf =
        sqrt((gamma + 1)/(gamma - 1))
       *atan(sqrt((gamma - 1)/(gamma + 1)*(sqr(MachInf) - 1)))
      - atan(sqr(MachInf) - 1);


    // Set the patch boundary condition based on the Mach number and direction
    // of the flow dictated by the boundary/free-stream pressure difference

    forAll(Up, facei)
    {
        if (pp[facei] >= pInf_) // If outflow
        {
            // Assume supersonic outflow and calculate the boundary velocity
            // according to ???

            scalar fpp =
                sqrt(sqr(MachInf) - 1)
               /(gamma*sqr(MachInf))*mag(Ut[facei])*log(pp[facei]/pInf_);

            Up[facei] = Ut[facei] + fpp*nHatInf[facei];

            // Calculate the Mach number of the boundary velocity
            scalar Mach = mag(Up[facei])/sqrt(gamma/ppsi[facei]);

            if (Mach <= 1) // If subsonic
            {
                // Zero-gradient subsonic outflow

                Up[facei] = U[facei];
                valueFraction()[facei] = 0;
            }
        }
        else // if inflow
        {
            // Calculate the Mach number of the boundary velocity
            // from the boundary pressure assuming constant total pressure
            // expansion from the free-stream
            scalar Mach =
                sqrt
                (
                    (2/(gamma - 1))*(1 + ((gamma - 1)/2)*sqr(MachInf))
                   *pow(pp[facei]/pInf_, (1 - gamma)/gamma)
                  - 2/(gamma - 1)
                );

            if (Mach > 1) // If supersonic
            {
                // Supersonic inflow is assumed to occur according to the
                // Prandtl-Meyer expansion process

                scalar nuMachf =
                    sqrt((gamma + 1)/(gamma - 1))
                   *atan(sqrt((gamma - 1)/(gamma + 1)*(sqr(Mach) - 1)))
                  - atan(sqr(Mach) - 1);

                scalar fpp = (nuMachInf - nuMachf)*mag(Ut[facei]);

                Up[facei] = Ut[facei] + fpp*nHatInf[facei];
            }
            else // If subsonic
            {
                FatalErrorIn
                (
                    "supersonicFreeStreamFvPatchVectorField::updateCoeffs()"
                )   << "unphysical subsonic inflow has been generated"
                    << abort(FatalError);
            }
        }
    }

    mixedFvPatchVectorField::updateCoeffs();
}


// Write
void supersonicFreeStreamFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("UInf") << UInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("pInf") << pInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("TInf") << TInf_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, supersonicFreeStreamFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
