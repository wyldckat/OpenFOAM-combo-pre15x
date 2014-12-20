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
    cuttingPlane

Description
    Cutting plane sampling functionality

\*---------------------------------------------------------------------------*/

#include "cuttingPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > cuttingPlane::sample
(
    const Field<Type>& sf
) const
{
    if (sf.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > cuttingPlane::sample"
            "(const Field<Type>& sf) const"
        )   << "The argument field does not correspond to the mesh. "
            << "Field size: " << sf.size()
            << " mesh size: " << mesh_.nCells()
            << abort(FatalError);
    }

    return tmp<Field<Type> >(new Field<Type>(sf, cells()));
}


template<class Type>
tmp<Field<Type> > cuttingPlane::sample
(
    const tmp<Field<Type> >& tsf
) const
{
    tmp<Field<Type> > tint = sample(tsf());
    tsf.clear();
    return tint;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
