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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList()
:
    List<T>(SizeInc),
    nextFree_(0)
{}


//- Construct given size
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const label s
)
:
    List<T>(s),
    nextFree_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::label Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::size()
const
{
    return nextFree_;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label s
)
{
    if (s < nextFree_)
    {
        nextFree_ = s;
    }
    else
    {
        List<T>::setSize(s);
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::clear()
{
    nextFree_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::shrink()
{
    List<T>::setSize(nextFree_);
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
    const T& e
)
{
    operator()(size()) = e;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T& Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator()
(
    const Foam::label i
)
{
    nextFree_ = max(nextFree_, i + 1);

    if (nextFree_ > List<T>::size())
    {
        List<T>::setSize
        (
            max
            (
                nextFree_,
                label(SizeMult*List<T>::size()/SizeDiv + SizeInc)
            )
        );
    }

    return this->operator[](i);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const T& t
)
{
    List<T>::operator=(t);
}


// ************************************************************************* //
