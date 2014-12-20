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

\*---------------------------------------------------------------------------*/

#include "directPolyTopoChange.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
void Foam::directPolyTopoChange::reorder
(
    const labelList& oldToNew,
    DynamicList<T>& lst
)
{
    // Create copy
    DynamicList<T> oldLst(lst);

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI] = oldLst[elemI];
        }
    }    
}


template <class T>
void Foam::directPolyTopoChange::reorder
(
    const labelList& oldToNew,
    List<DynamicList<T> >& lst
)
{
    // Create copy
    List<DynamicList<T> > oldLst(lst);

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI].transfer(oldLst[elemI]);
        }
    }    
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
