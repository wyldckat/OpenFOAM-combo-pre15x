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
    A dynamic list is a 1-D vector of objects of type T which resizes
    itself as necessary to accept the new objects.  Internal storage
    is a compact array and the list can be shrunk to compact storage.
    The increment of list size is controlled by two template parameters,
    which allows the list storage to either increate by the given increment
    or a given multiplier.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::DynamicList()
:
    List<T>(SizeIncrement),
    nextFree_(0)
{}


//- Construct given size
template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::DynamicList
(
    const label s
)
:
    List<T>(s),
    nextFree_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline Foam::label Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::size()
const
{
    return nextFree_;
}

template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline void Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::setSize
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

template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline void Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::clear()
{
    nextFree_ = 0;
}


template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline Foam::DynamicList<T, SizeIncrement, SizeMultiplier>&
Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::shrink()
{
    List<T>::setSize(nextFree_);
    return *this;
}


template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline void Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::append
(
    const T& e
)
{
    operator()(size()) = e;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline T& Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::operator()
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
                label(SizeMultiplier*List<T>::size() + SizeIncrement)
            )
        );
    }

    return this->operator[](i);
}


template<class T, unsigned SizeIncrement, unsigned SizeMultiplier>
inline void Foam::DynamicList<T, SizeIncrement, SizeMultiplier>::operator=
(
    const T& t
)
{
    List<T>::operator=(t);
}


// ************************************************************************* //
