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

#include "adjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorFvPatchFields.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::adjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& p
)
{
    if (p.needReference())
    {
        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;

        forAll (phi.boundaryField(), patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvPatchScalarField& phip = phi.boundaryField()[patchi];

            if (typeid(phip) != typeid(processorFvPatchScalarField))
            {
                if
                (
                    Up.fixesValue() && typeid(Up)
                 != typeid(inletOutletFvPatchVectorField)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }
                }
                else
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            adjustableMassOut += phip[i];
                        }
                    }
                }
            }
        }

        reduce(massIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());

        scalar massCorr = 1.0;

        if (mag(adjustableMassOut) > SMALL)
        {
            massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }
        else if(mag(fixedMassOut - massIn) > SMALL)
        {
            FatalErrorIn
            (
                "adjustPhi(surfaceScalarField& phi, const volVectorField& U,"
                "const volScalarField& p"
            )   << "Continuity error cannot be removed by adjusting the"
                   " outflow.\nPlease check the velocity boundary conditions"
                   " and/or run potentialFoam to initialise the outflow."
                << exit(FatalError);
        }

        forAll (phi.boundaryField(), patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvPatchScalarField& phip = phi.boundaryField()[patchi];

            if (typeid(phip) != typeid(processorFvPatchScalarField))
            {
                if (!Up.fixesValue())
                {
                    forAll(phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                    }
                }
            }
        }

        return mag(massIn) < SMALL
            && mag(fixedMassOut) < SMALL
            && mag(adjustableMassOut) < SMALL;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
