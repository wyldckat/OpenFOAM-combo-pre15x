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

#include "processorLduCoupledInterface.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void processorLduCoupledInterface::send
(
    const tmp<Field<Type> >& tf,
    const bool bufferdTransfer
) const
{
    OPstream::write
    (
        neighbProcNo(),
        reinterpret_cast<const char*>(tf().begin()),
        tf().size()*sizeof(Type),
        bufferdTransfer
    );
    tf.clear();
}

template<class Type>
tmp<Field<Type> > processorLduCoupledInterface::receive
(
    const label size
) const
{
    tmp<Field<Type> > tf(new Field<Type>(size));
    Field<Type>& f = tf();

    IPstream::read
    (
        neighbProcNo(), 
        reinterpret_cast<char*>(f.begin()),
        size*sizeof(Type)
    );

    return tf;
}


template<class Type>
void processorLduCoupledInterface::compressedSend
(
    const tmp<Field<Type> >& tf,
    const bool bufferdTransfer
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer)
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (tf().size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        const Field<Type>& f = tf();
        const scalar *sArray = reinterpret_cast<const scalar*>(f.begin());
        const scalar *slast = &sArray[nm1];
        float fArray[nFloats];

        for (register label i=0; i<nm1; i++)
        {
            fArray[i] = sArray[i] - slast[i%nCmpts];
        }

        reinterpret_cast<Type&>(fArray[nm1]) = f[tf().size() - 1];

        OPstream::write
        (
            neighbProcNo(),
            reinterpret_cast<const char*>(fArray),
            nBytes,
            bufferdTransfer
        );
        tf.clear();
    }
    else
    {
        this->send(tf, bufferdTransfer);
    }
}

template<class Type>
tmp<Field<Type> > processorLduCoupledInterface::compressedReceive
(
    const label size
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer)
    {
        tmp<Field<Type> > tf(new Field<Type>(size));
        Field<Type>& f = tf();

        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (size - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        float fArray[nFloats];

        IPstream::read(neighbProcNo(), reinterpret_cast<char*>(fArray), nBytes);

        f[size - 1] = reinterpret_cast<const Type&>(fArray[nm1]);
        scalar *sArray = reinterpret_cast<scalar*>(f.begin());
        const scalar *slast = &sArray[nm1];

        for (register label i=0; i<nm1; i++)
        {
            sArray[i] = fArray[i] + slast[i%nCmpts];
        }

        return tf;
    }
    else
    {
        return this->receive<Type>(size);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
