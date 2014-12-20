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
    wallNormalInfo

Description
    Holds information regarding nearest wall point. Used in wall refinement.

SourceFiles
    wallNormalInfoI.H
    wallNormalInfo.C

\*---------------------------------------------------------------------------*/

#ifndef wallNormalInfo_H
#define wallNormalInfo_H

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
                           Class wallNormalInfo Declaration
\*---------------------------------------------------------------------------*/

class wallNormalInfo
{
    // Private data

        //- Normal at nearest wall point
        vector normal_;

    // Private Member Functions

        //- Evaluate distance to point. Update normal_
        inline bool update(const wallNormalInfo& w2);

public:

    // Static data members

        //- initial point far away.
        static point greatVector;

    // Constructors

        //- Construct null
        inline wallNormalInfo();

        //- Construct from normal
        inline wallNormalInfo(const vector& normal);

        //- Construct as copy
        inline wallNormalInfo(const wallNormalInfo&);

    // Member Functions

        // Access

            inline const vector& normal() const;


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid() const;

            //- Check for identical geometrical data. Used for cyclics checking.
            inline bool sameGeometry(const wallNormalInfo&, const scalar) const;

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
                const wallNormalInfo& neighbourWallInfo,
                const scalar tol
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of faceI
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const label neighbourCellI,
                const wallNormalInfo& neighbourWallInfo,
                const scalar tol
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of faceI
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const wallNormalInfo& neighbourWallInfo,
                const scalar tol
            );

    // Member Operators

        //Note: Used to determine whether to call update.
        inline bool operator==(const wallNormalInfo&) const;

        inline bool operator!=(const wallNormalInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const wallNormalInfo&);
        friend Istream& operator>>(Istream&, wallNormalInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallNormalInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
