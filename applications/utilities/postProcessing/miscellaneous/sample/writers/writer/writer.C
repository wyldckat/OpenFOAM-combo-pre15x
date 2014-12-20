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

#include "writer.H"
#include "coordSet.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
autoPtr<writer<Type> > writer<Type>::New(const word& writeType)
{
    typename wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_
            ->find(writeType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "writer::New(const word&)"
        )   << "Unknown write type " << writeType
            << endl << endl
            << "Valid write types : " << endl
            << wordConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<writer<Type> >(cstrIter()());
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
fileName writer<Type>::getBaseName
(
    const coordSet& points,
    const HashTable<Field<Type>*>& valueSets
) const
{
    fileName fName(points.name());

    typename HashTable<Field<Type>*>::const_iterator iter =
        valueSets.begin();

    for(; iter != valueSets.end(); ++iter)
    {
        fName += '_' + iter.key();
    }

    return fName;
}


template<class Type>
void writer<Type>::writeCoord
(
    const coordSet& points,
    const label pointI,
    Ostream& os
) const
{
    if (points.hasVectorAxis())
    {
        write(points.vectorCoord(pointI), os);
    }
    else
    {
        write(points.scalarCoord(pointI), os);
    }
}


template<class Type>
void writer<Type>::writeTable
(
    const coordSet& points,
    const List<Type>& values,
    Ostream& os
) const
{
    forAll(points, pointI)
    {
        writeCoord(points, pointI, os);

        os << token::SPACE;

        write(values[pointI], os);

        os << endl;
    }
}


template<class Type>
void writer<Type>::writeTable
(
    const coordSet& points,
    const List<List<Type>*>& valuesPtrList,
    Ostream& os
) const
{
    forAll(points, pointI)
    {
        writeCoord(points, pointI, os);

        forAll(valuesPtrList, i)
        {
            os << token::SPACE;

            const List<Type>& values = *valuesPtrList[i];

            write(values[pointI], os);
        }
        os << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template<class Type>
writer<Type>::writer()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
writer<Type>::~writer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::writer<Type>::write(const scalar value, Ostream& os) const
{
    return os << value;
}


template<class Type>
Foam::Ostream& Foam::writer<Type>::write(const vector& value, Ostream& os) const
{
    for (direction d=0; d<vector::nComponents; d++)
    {
        os << value.component(d);

        if (d <= vector::nComponents-1)
        {
            os << token::TAB;
        }
    }
    return os;
}


template<class Type>
Foam::Ostream& Foam::writer<Type>::write(const tensor& value, Ostream& os) const
{
    for (direction d=0; d<tensor::nComponents; d++)
    {
        os << value.component(d);

        if (d <= tensor::nComponents-1)
        {
            os << token::TAB;
        }
    }
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
