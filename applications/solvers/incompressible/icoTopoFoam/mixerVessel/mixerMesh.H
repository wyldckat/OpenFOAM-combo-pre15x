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
    mixerMesh

Description

SourceFiles
    mixerMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mixerMesh_H
#define mixerMesh_H

#include "fvMesh.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class mixerMesh Declaration
\*---------------------------------------------------------------------------*/

class mixerMesh
:
    public fvMesh
{
    // Private data

        //- Coordinate system
        cylindricalCS cs_;

        // - Rotational speed in rotations per minute (rpm)
        scalar rpm_;

        //- Name of cell zone containing moving cells
        word movingCellsName_;

        //- Name of face zone containing inner sliding surface
        word innerSliderName_;

        //- Name of face zone containing outer sliding surface
        word outerSliderName_;

        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixerMesh(const mixerMesh&);

        //- Disallow default bitwise assignment
        void operator=(const mixerMesh&);

        //- Calculate moving masks
        void calcMovingMasks() const;


public:

    // Constructors

        //- Construct from IOobject and dictionary
        mixerMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


    // Destructor

        ~mixerMesh();


    // Member Functions

        //- Return coordinate system
        const cylindricalCS& cs() const
        {
            return cs_;
        }

        //- Return moving points mask.  Moving cells marked with 1
        const scalarField& movingPointsMask() const;

        //- Slide
        void slide();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
