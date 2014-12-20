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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct given size
template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const Foam::UList<T>& completeList,
    const Foam::List<label>& addr
)
:
    completeList_(completeList),
    addressing_(addr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::IndirectList<T>::size() const
{
    return addressing_.size();
}


template<class T>
inline const Foam::UList<T>& Foam::IndirectList<T>::
completeList() const
{
    return completeList_;
}


template<class T>
inline const Foam::List<Foam::label>& Foam::IndirectList<T>::addressing() const
{
    return addressing_;
}


template<class T>
inline void Foam::IndirectList<T>::resetAddressing
(
    const List<label>& addr
)
{
    addressing_ = addr;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::List<T> Foam::IndirectList<T>::operator()() const
{
    List<T> result(size());

    forAll(*this, i)
    {
        result[i] = operator[](i);
    }

    return result;
}


template<class T>
inline const T& Foam::IndirectList<T>::operator[](const Foam::label i) const
{
    return completeList_[addressing_[i]];
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const IndirectList<T>& L)
{
    // Write size of list and start contents delimiter
    os << nl << L.size() << nl << token::BEGIN_LIST;

    // Write list contents
    forAll(L, i)
    {
        os << nl << L[i];
    }

    // Write end of contents delimiter
    os << nl << token::END_LIST << nl;

    return os;
}

// ************************************************************************* //
