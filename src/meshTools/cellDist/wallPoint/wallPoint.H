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
    wallPoint

Description
    Holds information regarding nearest wall point. Used in wall distance
    calculation.

SourceFiles
    wallPointI.H
    wallPoint.C

\*---------------------------------------------------------------------------*/

#ifndef wallPoint_H
#define wallPoint_H

#include "point.H"
#include "label.H"
#include "scalar.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class wallPoint Declaration
\*---------------------------------------------------------------------------*/

class wallPoint
{
    // Private data

        //- position of nearest wall center
        point origin_;

        //- normal distance (squared) from cellcenter to origin
        scalar distSqr_;

    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        inline bool update
        (
            const point&,
            const wallPoint& w2,
            const scalar tol
        );

public:
    // Static data members

        //- initial point far away.
        static point greatPoint;

    // Constructors

        //- Construct null
        inline wallPoint();

        //- Construct from origin, distance
        inline wallPoint
        (
            const point& origin,
            const scalar distSqr
        );

        //- Construct as copy
        inline wallPoint
        (
            const wallPoint&
        );

    // Member Functions

        // Access

            inline const point& origin() const;

            inline point& origin();

            inline scalar distSqr() const;

            inline scalar& distSqr();


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid() const;

            //- Check for identical geometrical data. Used for cyclics checking.
            inline bool sameGeometry(const wallPoint&, const scalar tol) const;

            //- Convert origin to relative vector to leaving point
            //  (= face centre)
            inline void leaveDomain
            (
                const polyPatch& patch,
                const label patchFaceI,
                const point& faceCentre
            );

            //- Convert relative origin to absolute by adding entering point
            inline void enterDomain
            (
                const polyPatch& patch,
                const label patchFaceI,
                const point& faceCentre
            );

            //- Apply rotation matrix to origin
            inline void transform(const tensor& rotTensor);

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of cellI
            inline bool updateCell
            (
                const polyMesh& mesh,
                const label thisCellI,
                const label neighbourFaceI,
                const wallPoint& neighbourWallInfo,
                const scalar tol
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of faceI
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const label neighbourCellI,
                const wallPoint& neighbourWallInfo,
                const scalar tol
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of faceI
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const wallPoint& neighbourWallInfo,
                const scalar tol
            );

    // Member Operators

        //Note: Used to determine whether to call update.
        inline bool operator==(const wallPoint&) const;

        inline bool operator!=(const wallPoint&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const wallPoint&);
        friend Istream& operator>>(Istream&, wallPoint&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
