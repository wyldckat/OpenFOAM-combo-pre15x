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
    staticTopoFvMesh

Description

SourceFiles
    linearValveFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef linearValveLayersFvMesh_H
#define linearValveLayersFvMesh_H

#include "topoFvMesh.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                       Class linearValveLayersFvMesh Declaration
\*---------------------------------------------------------------------------*/

class linearValveLayersFvMesh
:
    public topoFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearValveLayersFvMesh(const linearValveLayersFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const linearValveLayersFvMesh&);


        //- Add valve and piston zones and modifiers
        void addZonesAndModifiers();

        //- Make layering modifiers live
        void makeLayersLive();

        //- Make sliding modifiers live
        void makeSlidersLive();

        //- Return true if sliders are attached
        bool attached() const;

        //- Return new point field
        tmp<pointField> newPoints() const;


public:

    //- Runtime type information
    TypeName("linearValveLayersFvMesh");


    // Constructors

        //- Construct from database
        explicit linearValveLayersFvMesh(const IOobject& io);


    // Destructor

        virtual ~linearValveLayersFvMesh();


    // Member Functions

        //- Move and morph
        virtual void moveAndMorph();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
