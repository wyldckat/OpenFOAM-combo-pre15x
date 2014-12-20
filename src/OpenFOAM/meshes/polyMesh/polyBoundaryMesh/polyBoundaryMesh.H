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
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef List<labelPair> labelPairList;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMesh;
class mapPolyMesh;

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

        //- Whether extended (mostly parallel) addressing needs to be
        //  calculated and maintained.
        bool extendedAddressing_;

        //- Edges of neighbouring patches
        mutable List<labelPairList>* neighbourEdgesPtr_;


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


    // Destructor

        ~polyBoundaryMesh();

        //- Clear geometry at this level and at patches
        void clearGeom();

        //- Clear addressing at this level and at patches
        void clearAddressing();


    // Member functions

        //- Calculate the geometry for the patches (transformation tensors etc.)
        void calcGeometry();

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Per patch the edges on the neighbouring patch. Is for every external
        //  edge the neighbouring patch and neighbouring (external) patch edge
        //  label. Note that edge indices are offset by nInternalEdges to keep
        //  it as much as possible consistent with coupled patch addressing
        //  (where coupling is by local patch face index).
        //  Only valid for singly connected polyBoundaryMesh and not parallel
        const List<labelPairList>& neighbourEdges() const;

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

        //- Correct polyBoundaryMesh after topology update
        void updateTopology();

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;

        //- Write using given format, version and form uncompression
        bool write
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;

        //- Calculate extended addressing yes or no.
        void setExtendedAddressing(const bool extended)
        {
            extendedAddressing_ = extended;
        }

        bool extendedAddressing() const
        {
            return extendedAddressing_;
        }


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const polyBoundaryMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
