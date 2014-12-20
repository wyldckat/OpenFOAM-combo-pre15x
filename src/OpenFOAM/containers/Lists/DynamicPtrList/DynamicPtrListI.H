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

//- Construct null
template<class T>
inline Foam::DynamicList<T>::DynamicList()
:
    List<T>(sizeIncrement_)
{}


//- Construct given size
template<class T>
inline Foam::DynamicList<T>::DynamicList(const label s)
:
    List<T>(s)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::DynamicList<T>::size() const
{
    return nextFree();
}

template<class T>
inline Foam::label Foam::DynamicList<T>::setIncrement
(
    const Foam::label sizeIncrement
)
{
    label oldSizeIncrement = sizeIncrement_;
    sizeIncrement_ = sizeIncrement;
    sizeMultiplier_ = 1;
    return oldSizeIncrement;
}

template<class T>
inline Foam::label Foam::DynamicList<T>::setMultiplier
(
    const Foam::label sizeMultiplier
)
{
    label oldSizeMultiplier = sizeMultiplier_;
    sizeMultiplier_ = sizeMultiplier;
    sizeIncrement_ = 0;
    return oldSizeMultiplier;
}


template<class T>
inline Foam::DynamicList<T>& Foam::DynamicList<T>::shrink()
{
    setSize(nextFree());
    return *this;
}


template<class T>
inline void Foam::DynamicList<T>::append(const T& e)
{
    operator()(size()) = e;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::DynamicList<T>::operator()(const Foam::label i)
{
    label used = max(nextFree(), i + 1);

    if (used > ptrList<T>::size())
    {
        setSize
        (
            max(used, sizeMultiplier_*ptrList<T>::size() + sizeIncrement_)
        );
    }

    return operator[](i);
}


template<class T>
inline void Foam::DynamicList<T>::operator=(const T& t)
{
    List<T>::operator=(t);
}


// ************************************************************************* //
