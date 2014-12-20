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
    Tetrahedral Finite Element matrix check for diagonal dominance.
    Especially useful for debugging complex interfaces

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

template<class Type>
void tetFemMatrix<Type>::check()
{
#   ifdef DEBUGtetFemMatrix
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::check : checking tetFemMatrix<Type>"
            << endl;
    }
#   endif

    Info << "first check diagonal dominance" << endl;
    // Calculate local matrix off-diag sum
    scalarField dummySource(diag().size(), 0.0);

    const scalarField oldLower = ((const lduMatrix&)(*this)).lower();
    const scalarField oldUpper = ((const lduMatrix&)(*this)).upper();
    const scalarField oldDiag = ((const lduMatrix&)(*this)).diag();

    // Get matrix addressing
    const unallocLabelList& L = lduAddr().lowerAddr();
    const unallocLabelList& U = lduAddr().upperAddr();

    {
        scalarField matrixSumOffDiag(diag().size(), 0.0);
        scalarField dummyInternal(diag().size(), 1.0);

        forAll (L, face)
        {
            matrixSumOffDiag[L[face]] += oldLower[face];
            matrixSumOffDiag[U[face]] += oldUpper[face];
        }

        Info
            << "raw matrix difference: "
            << sum(mag(matrixSumOffDiag + oldDiag)) << endl;
    }

    // Add the coupling coefficients
    addCouplingCoeffs();

    addCouplingSource(dummySource);

    // Prepare for coupled interface update
    FieldField<Field, scalar> coupledBouCoeffs(psi_.boundaryField().size());
    FieldField<Field, scalar> coupledIntCoeffs(psi_.boundaryField().size());
    lduCoupledInterfacePtrsList interfaces(psi_.boundaryField().size());

    forAll(psi_.boundaryField(), patchI)
    {
        const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        coupledBouCoeffs.hook
        (
            new scalarField(ptf.cutBouCoeffs(*this))
        );

        coupledIntCoeffs.hook
        (
            new scalarField(ptf.cutIntCoeffs(*this))
        );

        interfaces[patchI] = &psi_.boundaryField()[patchI];
    }

    eliminateCouplingCoeffs();

    Info << "second check diagonal dominance" << endl;
    // Calculate local matrix off-diag sum
    {
        scalarField matrixSumOffDiag(diag().size(), 0.0);
        scalarField dummyInternal(diag().size(), 1.0);

        const scalarField& lower = ((const lduMatrix&)(*this)).lower();
        const scalarField& upper = ((const lduMatrix&)(*this)).upper();

        forAll (L, face)
        {
            matrixSumOffDiag[L[face]] += lower[face];
            matrixSumOffDiag[U[face]] += upper[face];
        }

        // Initialise coupling
        forAll (interfaces, interfaceI)
        {
            if (interfaces[interfaceI]->coupled())
            {
                interfaces[interfaceI]->initInterfaceMatrixUpdate
                (
                    dummyInternal,
                    matrixSumOffDiag,
                    *this,
                    coupledBouCoeffs[interfaceI],
                    0
                );
            }
        }

        // Update coupled interface
        forAll (interfaces, interfaceI)
        {
            if (interfaces[interfaceI]->coupled())
            {
                interfaces[interfaceI]->updateInterfaceMatrix
                (
                    dummyInternal,
                    matrixSumOffDiag,
                    *this,
                    coupledBouCoeffs[interfaceI],
                    0
                );
            }
        }

        Info
            << "parallel matrix difference: "
            << sum(mag(matrixSumOffDiag + diag())) << endl;
    }

    // restore the matrix in the original state
    if (thereIsUpper())
    {
        upper() = oldUpper;
    }

    if (thereIsLower())
    {
        lower() = oldLower;
    }

    if (thereIsDiag())
    {
        diag() = oldDiag;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
