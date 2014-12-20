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

#include "xmgr.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::xmgr<Type>::xmgr()
:
    writer<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::xmgr<Type>::~xmgr()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::xmgr<Type>::getFileName
(
    const coordSet& points,
    const Foam::HashTable<Field<Type>*>& valueSets
) const
{
    return getBaseName(points, valueSets) + ".agr";
}


template<class Type>
void Foam::xmgr<Type>::write
(
    const coordSet& points,
    const HashTable<Field<Type>*>& valueSets,
    Ostream& os
) const
{
    os  << "@title \"" << points.name() << '"' << endl
        << "@xaxis label " << '"' << points.axis() << '"' << endl;

    label fieldI = 0;

    typename HashTable<Field<Type>*>::const_iterator iter = valueSets.begin();

    for(; iter != valueSets.end(); ++iter)
    {
        os  << "@s" << fieldI << " legend " << '"'
            << iter.key() << '"' << endl
            << "@target G0.S" << fieldI << endl
            << "@type xy" << endl;

        writeTable(points, *iter(), os);

        os << endl;

        fieldI++;
    }
}


// ************************************************************************* //
