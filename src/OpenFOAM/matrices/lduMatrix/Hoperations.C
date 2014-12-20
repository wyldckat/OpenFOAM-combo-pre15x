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
    lduMatrix member H operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// H operator
template<class Type>
tmp<Field<Type> > lduMatrix::H(const Field<Type>& sf) const
{
    tmp<Field<Type> > tHphi
    (
        new Field<Type>(lduAddr_.size(), pTraits<Type>::zero)
    );

    if (lowerPtr_ || upperPtr_)
    {
        Field<Type> & Hphi = tHphi();

        const scalarField& Lower = lower();
        const scalarField& Upper = upper();

        // Take refereces to addressing
        const unallocLabelList& l = lduAddr_.lowerAddr();
        const unallocLabelList& u = lduAddr_.upperAddr();

        for (register label face=0; face<l.size(); face++)
        {
            Hphi[u[face]] -= Lower[face]*sf[l[face]];
            Hphi[l[face]] -= Upper[face]*sf[u[face]];
        }
    }

    return tHphi;
}

template<class Type>
tmp<Field<Type> > lduMatrix::H(const tmp<Field<Type> >& tsf) const
{
    tmp<Field<Type> > tHphi(H(tsf()));
    tsf.clear();
    return tHphi;
}


// face H operator
template<class Type>
tmp<Field<Type> > lduMatrix::faceH(const Field<Type>& sf) const
{
    const scalarField& Lower = ((const lduMatrix&)(*this)).lower();
    const scalarField& Upper = ((const lduMatrix&)(*this)).upper();

    // Take refereces to addressing
    const unallocLabelList& l = lduAddr_.lowerAddr();
    const unallocLabelList& u = lduAddr_.upperAddr();

    tmp<Field<Type> > tfaceHphi(new Field<Type> (Lower.size()));
    Field<Type> & faceHphi = tfaceHphi();

    for (register label face=0; face<l.size(); face++)
    {
        faceHphi[face] = Upper[face]*sf[u[face]] - Lower[face]*sf[l[face]];
    }

    return tfaceHphi;
}

template<class Type>
tmp<Field<Type> > lduMatrix::faceH(const tmp<Field<Type> >& tsf) const
{
    tmp<Field<Type> > tfaceHphi(faceH(tsf()));
    tsf.clear();
    return tfaceHphi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
