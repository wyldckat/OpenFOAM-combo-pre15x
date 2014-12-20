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

#include "raw.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::raw<Type>::raw()
:
    writer<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::raw<Type>::~raw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::raw<Type>::getFileName
(
    const coordSet& points,
    const Foam::HashTable<Field<Type>*>& valueSets
) const
{
    return getBaseName(points, valueSets) + ".xy";
}


template<class Type>
void Foam::raw<Type>::write
(
    const coordSet& points,
    const HashTable<Field<Type>*>& valueSets,
    Ostream& os
) const
{
    // Collect sets into columns
    List<List<Type>*> columns(valueSets.size());

    label fieldI = 0;

    typename HashTable<Field<Type>*>::const_iterator iter = valueSets.begin();

    for(; iter != valueSets.end(); ++iter)
    {
        columns[fieldI++] = iter();
    }

    writeTable(points, columns, os);
}


// ************************************************************************* //
