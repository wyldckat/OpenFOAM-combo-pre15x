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
    coordinateSystem

Description
    Virtual base class for coordinate systems.

    All systems are defined in two ways:
    1) In terms of a position (the coordinate origin), and axis
    (usually the z axis) and a direction (effectively the x-axis).
    The code requires that the axis and the direction are orthogonal.
    2) Given the origin and a coordinate rotation.

SourceFiles
    coordinateSystem.C

\*---------------------------------------------------------------------------*/

#ifndef coordinateSystem_H
#define coordinateSystem_H

#include "vector.H"
#include "tensor.H"
#include "vectorField.H"
#include "tmp.H"
#include "coordinateRotation.H"
#include "physicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coordinateSystem Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystem
{
    // Private data

        //- Name of patch
        word name_;

        //- Origin
        vector origin_;

        //- Axis
        vector axis_;

        //- Origin
        vector dir_;


public:

    //- Runtime type information
    TypeName("coordinateSystem");


    // Constructors

        //- Construct from origin and 2 axes
        coordinateSystem
        (
            const word& name,
            const vector& origin,
            const vector& axis,
            const vector& dir
        );

        //- Construct from origin and rotation angles
        coordinateSystem
        (
            const word& name,
            const vector& origin,
            const coordinateRotation& cr
        );

        //- Construct from dictionary
        coordinateSystem(const word& name, const dictionary& dict);


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            coordinateSystem,
            origAxisDir,
            (
                const word& name,
                const vector& origin,
                const vector& axis,
                const vector& dir
            ),
            (name, origin, axis, dir)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            coordinateSystem,
            origRotation,
            (
                const word& name,
                const vector& origin,
                const coordinateRotation& cr
            ),
            (name, origin, cr)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            coordinateSystem,
            dictionary,
            (
                const word& name,
                const dictionary& dict
            ),
            (name, dict)
        );


    // Selectors

        //- Select constructed from origin and 2 axes
        static autoPtr<coordinateSystem> New
        (
            const word& coordType,
            const word& name,
            const vector& origin,
            const vector& axis,
            const vector& dir
        );

        //- Select constructed from origin and rotation
        static autoPtr<coordinateSystem> New
        (
            const word& coordType,
            const word& name,
            const vector& origin,
            const coordinateRotation& cr
        );

        //- Select constructed from Istream
        static autoPtr<coordinateSystem> New
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        virtual ~coordinateSystem();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return origin
        const vector& origin() const
        {
            return origin_;
        }

        //- Return axis
        const vector& axis() const
        {
            return axis_;
        }

        //- Return direction
        const vector& direction() const
        {
            return dir_;
        }


        //- Convert from local coordinate system to the global Cartesian system
        virtual vector toGlobal(const vector& localV) const = 0;

        virtual tmp<vectorField> toGlobal(const vectorField& localV) const = 0;

        //- Convert from global Cartesian system to the local coordinate system 
        virtual vector toLocal(const vector& globalV) const = 0;

        virtual tmp<vectorField> toLocal(const vectorField& globalV) const = 0;

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const coordinateSystem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
