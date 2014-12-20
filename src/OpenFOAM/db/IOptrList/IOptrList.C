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
    An IOptrList of a given type is a ptrList of that type which supports
    automated input and output.

\*---------------------------------------------------------------------------*/

#include "IOptrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
template<class INew>
IOptrList<T>::IOptrList(const IOobject& io, const INew& inewt)
:
    regIOobject(io),
    ptrList<T>(readStream(typeName), inewt)
{
    close();
}


template<class T>
IOptrList<T>::IOptrList(const IOobject& io)
:
    regIOobject(io),
    ptrList<T>(readStream(typeName))
{
    close();
}


template<class T>
IOptrList<T>::IOptrList(const IOobject& io, const ptrList<T>& list)
:
    regIOobject(io),
    ptrList<T>(list)
{}


template<class T>
void IOptrList<T>::operator=(const IOptrList<T>& rhs)
{
    ptrList<T>::operator=(rhs);
}


template<class T>
bool IOptrList<T>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
