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

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void lduMatrix::initMatrixInterfaces
(
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt
) const
{
    forAll (lduCoupledInterfaceSchedule_, i)
    {
        label interfaceI = lduCoupledInterfaceSchedule_[i].patch;

        if (interfaces[interfaceI]->coupled())
        {
            if
            (
                lduCoupledInterfaceSchedule_[i].init
             && lduCoupledInterfaceSchedule_[i].bufferedTransfer
            )
            {
                interfaces[interfaceI]->initInterfaceMatrixUpdate
                (
                    psiif,
                    result,
                    *this,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    lduCoupledInterfaceSchedule_[i].bufferedTransfer
                );
            }
        }
    }

    for
    (
        label interfaceI=lduCoupledInterfaceSchedule_.size()/2;
        interfaceI<interfaces.size();
        interfaceI++
    )
    {
        if (interfaces[interfaceI]->coupled())
        {
            interfaces[interfaceI]->initInterfaceMatrixUpdate
            (
                psiif,
                result,
                *this,
                coupleCoeffs[interfaceI],
                cmpt,
                true
            );
        }
    }
}


void lduMatrix::updateMatrixInterfaces
(
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt
) const
{
    forAll (lduCoupledInterfaceSchedule_, i)
    {
        label interfaceI = lduCoupledInterfaceSchedule_[i].patch;

        if (interfaces[interfaceI]->coupled())
        {
            if (lduCoupledInterfaceSchedule_[i].init)
            {
                if (!lduCoupledInterfaceSchedule_[i].bufferedTransfer)
                {
                    interfaces[interfaceI]->initInterfaceMatrixUpdate
                    (
                        psiif,
                        result,
                        *this,
                        coupleCoeffs[interfaceI],
                        cmpt,
                        lduCoupledInterfaceSchedule_[i].bufferedTransfer
                    );
                }
            }
            else
            {
                interfaces[interfaceI]->updateInterfaceMatrix
                (
                    psiif,
                    result,
                    *this,
                    coupleCoeffs[interfaceI],
                    cmpt
                );
            }
        }
    }

    for
    (
        label interfaceI=lduCoupledInterfaceSchedule_.size()/2;
        interfaceI<interfaces.size();
        interfaceI++
    )
    {
        if (interfaces[interfaceI]->coupled())
        {
            interfaces[interfaceI]->updateInterfaceMatrix
            (
                psiif,
                result,
                *this,
                coupleCoeffs[interfaceI],
                cmpt
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
