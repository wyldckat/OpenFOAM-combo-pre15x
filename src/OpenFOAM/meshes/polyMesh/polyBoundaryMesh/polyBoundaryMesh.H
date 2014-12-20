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
    polyBoundaryMesh

Description

SourceFiles
    polyBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef polyBoundaryMesh_H
#define polyBoundaryMesh_H

#include "polyPatchList.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class polyBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class polyBoundaryMesh
:
    public polyPatchList,
    public regIOobject
{
    // private data

        //- Reference to mesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- Create identity map
        static labelList ident(const label len);

        //- Disallow construct as copy
        polyBoundaryMesh(const polyBoundaryMesh&);

        //- Disallow assignment
        void operator=(const polyBoundaryMesh&);


public:

    //- Runtime type information
    TypeName("polyBoundaryMesh");


    // Constructors

        //- Read constructor given IOobject and a polyMesh reference
        //  Note point pointers are unset, only used in copying meshes
        polyBoundaryMesh
        (
            const IOobject&,
            const polyMesh&
        );

        //- Construct given size
        polyBoundaryMesh
        (
            const IOobject&,
            const polyMesh&,
            const label size
        );


    // Destructor - default


    // Member functions

        //- Calculate the geometry for the patches (transformation tensors etc.)
        void calcGeometry();

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return a list of patch names
        wordList names() const;

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of physical types
        wordList physicalTypes() const;

        //- Find patch index given a name
        label findPatchID(const word& patchName) const;

        //- Return patch index for a given face label
        label whichPatch(const label faceIndex) const;

        //- Check boundary definition
        bool checkDefinition(const bool report = false) const;

        //- Check face order on all patches. Note: does parallel communication.
        bool checkFaceOrder(const bool report = false) const;

        //- Correct polyBoundaryMesh after moving points
        void movePoints(const pointField&);

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const polyBoundaryMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
