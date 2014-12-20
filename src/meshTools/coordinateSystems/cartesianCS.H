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
    cartesianCS

Description
    Cartesian coordinate system.

SourceFiles
    cartesianCS.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianCS_H
#define cartesianCS_H

#include "coordinateSystem.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cartesianCS Declaration
\*---------------------------------------------------------------------------*/

class cartesianCS
:
    public coordinateSystem
{
    // Private data

        //- Local-to-global transformation tensor
        tensor R_;

        //- Global-to-local transformation tensor
        tensor Rtrans_;

        static const scalar nonOrthogonalError;


    // Private member functions

        //- Calculate transformation tensors
        void calcTransformations();


public:

    //- Runtime type information
    TypeName("Cartesian");


    // Constructors

        //- Construct from origin and 2 axes
        cartesianCS
        (
            const word& name,
            const vector& origin,
            const vector& axis,
            const vector& dir
        );

        //- Construct from origin and rotation angles
        cartesianCS
        (
            const word& name,
            const vector& origin,
            const coordinateRotation& cr
        );

        //- Construct from dictionary
        cartesianCS(const word& name, const dictionary& dict);


    // Destructor

        virtual ~cartesianCS()
        {}


    // Member Functions

        //- Convert from local coordinate system to the global Cartesian system
        virtual vector toGlobal(const vector& localV) const;

        virtual tmp<vectorField> toGlobal(const vectorField& localV) const;

        //- Convert from global Cartesian system to the local coordinate system 
        virtual vector toLocal(const vector& globalV) const;

        virtual tmp<vectorField> toLocal(const vectorField& globalV) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
