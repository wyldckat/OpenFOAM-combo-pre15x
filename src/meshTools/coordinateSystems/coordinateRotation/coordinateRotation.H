/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    coordinateRotation

Description
    Virtual base class for coordinate rotation.
    Allows multiple definitions of rotations of coordinate sytems.

\*---------------------------------------------------------------------------*/

#ifndef coordinateRotation_H
#define coordinateRotation_H

#include "vector.H"
#include "tensor.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                       Class coordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class coordinateRotation
{

protected:

    // Protected data members

        //- Local-to-global transformation tensor
        tensor R_;


    // Protected Member Functions

        //- Calculate transformation tensors
        virtual void calcTransformations
        (
            const scalar,
            const scalar,
            const scalar,
            const bool inDegrees
        ) = 0;


public:

    //- Runtime type information
    TypeName("coordinateRotation");


    // Static data members

        //- Axis
        static const vector z;

        //- Direction
        static const vector x;


    // Constructors

        //- Construct null
        coordinateRotation();


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            coordinateRotation,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Select constructed from Istream
        static autoPtr<coordinateRotation> New
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~coordinateRotation()
        {}


    // Member Functions

        //- Return local-to-global transformation tensor
        const tensor& R() const
        {
            return R_;
        }

        //- Return local axis
        vector axis() const
        {
            return R_ & z;
        }

        //- Return local direction
        vector direction() const
        {
            return R_ & x;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
