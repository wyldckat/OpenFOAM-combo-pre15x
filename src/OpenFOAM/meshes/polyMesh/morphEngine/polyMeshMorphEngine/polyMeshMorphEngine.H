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
    polyMeshMorphEngine

Description
    List of mesh modifiers defining the mesh dynamics.

SourceFiles
    polyMeshMorphEngine.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshMorphEngine_H
#define polyMeshMorphEngine_H

#include "regIOobject.H"
#include "PtrList.H"
#include "polyMeshModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                      Class polyMeshMorphEngine Declaration
\*---------------------------------------------------------------------------*/

class polyMeshMorphEngine
:
    public PtrList<polyMeshModifier>,
    public regIOobject
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        polyMeshMorphEngine(const polyMeshMorphEngine&);

        //- Disallow default bitwise assignment
        void operator=(const polyMeshMorphEngine&);


public:

    //- Runtime type information
    TypeName("meshMorphEngine");


    // Constructors

        //- Read constructor given IOobject and a polyMesh reference
        polyMeshMorphEngine
        (
            const IOobject&,
            const polyMesh&
        );

        //- Construct given size
        polyMeshMorphEngine
        (
            const IOobject&,
            const polyMesh&,
            const label size
        );


    // Destructor - default


    // Member functions

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of patch names
        wordList names() const;

        //- Is topology change required
        bool changeTopology() const;

        //- Return topology change request
        tmp<polyTopoChange> topoChangeRequest() const;

        //- Modify point motion
        void modifyMotionPoints(pointField&) const;

        //- Force recalculation of locally stored data on topological change
        void updateTopology(const mapPolyMesh& m);


        //- Find modifier given a name
        label findModifierID(const word& modName) const;


        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const polyMeshMorphEngine&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
