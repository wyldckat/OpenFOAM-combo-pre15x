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
    emptyFaPatch

Description
    A patch which will not exist in the faMesh. Typical example is a front and
    back plane of a 2-D geometry

SourceFiles
    emptyFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef emptyFaPatch_H
#define emptyFaPatch_H

#include "faPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class emptyFaPatch Declaration
\*---------------------------------------------------------------------------*/

class emptyFaPatch
:
    public faPatch
{

public:

    //- Runtime type information
    TypeName("empty");


    // Constructors

        //- Construct from polyPatch
        emptyFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            faPatch(name, dict, index, bm)
        {}

    // Member Functions

        virtual label size() const
        {
            return 0;
        }

        //- Return face normals. Over-riding base class return to get zero size
        // 
//         virtual const vectorField& edgeNormals() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
