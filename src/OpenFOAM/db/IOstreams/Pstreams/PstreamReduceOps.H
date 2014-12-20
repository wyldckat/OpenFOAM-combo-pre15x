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
    Pstream

Description

\*---------------------------------------------------------------------------*/

#ifndef PstreamReduceOps_H
#define PstreamReduceOps_H

#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Reduce operation with user specified communication schedule
template <class T, class BinaryOp>
void reduce
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    BinaryOp bop
)
{
    Pstream::gather(comms, Value, bop);
    Pstream::scatter(comms, Value);
}    


// Reduce using either linear or tree communication schedule
template <class T, class BinaryOp>
void reduce
(
    T& Value,
    BinaryOp bop
)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum)
    {
        reduce(Pstream::linearCommunication(), Value, bop);
    }
    else
    {
        reduce(Pstream::treeCommunication(), Value, bop);
    }
}    


template <class T>
class sumOp
{

public:

    T operator()(const T& x, const T& y) const
    {
        return x + y;
    }
};


template <class T>
class productOp
{

public:

    T operator()(const T& x, const T& y) const
    {
        return x * y;
    }
};


template <class T>
class maxOp
{

public:

    T operator()(const T& x, const T& y) const
    {
        return max(x, y);
    }
};


template <class T>
class minOp
{

public:

    T operator()(const T& x, const T& y) const
    {
        return min(x, y);
    }
};


template <class T>
class andOp
{

public:

    T operator()(const T& x, const T& y) const
    {
        return x && y;
    }
};


template <class T>
class orOp
{

public:

    T operator()(const T& x, const T& y) const
    {
        return x || y;
    }
};


// Insist there is a specialisation for the reduction of a scalar
void reduce(scalar& Value, sumOp<scalar> bop);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
