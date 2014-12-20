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
    plane

Description
    Geometric class that creates a 2D plane and can return the intersection
     point between a line and the plane.

SourceFiles
    plane.C

\*---------------------------------------------------------------------------*/

#ifndef plane_H
#define plane_H

#include "point.H"
#include "scalarList.H"
#include "dictionary.H"
#include "line.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class plane Declaration
\*---------------------------------------------------------------------------*/

class plane
{
    // Private data

        //- Plane normal
        vector unitVector_;

        //- Base point
        point basePoint_;


    // Private Member Functions

        //- Calculates basePoint and normal vector given plane coefficients
        void calcPntAndVec(const scalarList& C);

        //- Calculates basePoint and normal vector given three points
        //- Normal vector determined using right hand rule
        void calcPntAndVec
        (
            const point& point1,
            const point& point2,
            const point& point3
        );

public:

    // Constructors

        //- Construct null
        plane();
 
        //- Construct from normal vector and point in plane
        plane(const point& basePoint, const vector& normalVector);
 
        //- Construct from three points in plane
        plane(const point& point1, const point& point2, const point& point3);

        //- Construct from coefficients for the
        //  plane equation: ax + by + cz + d = 0
        plane(const scalarList& C);

        //- Construct from dictionary
        plane(const dictionary& planeDict);


    // Member Functions

        //- Return plane normal
        const vector& normal() const;

        //- Return or return plane base point
        const point& refPoint() const;

        //- Return coefficients for the
        //  plane equation: ax + by + cz + d = 0
        scalarList planeCoeffs() const;

        //- Return nearest point in the plane for the given point
        point nearestPoint(const point& p) const;

        //- Return distance from the given point to the plane
        scalar distance(const point& p) const;

        //- Return cut coefficient for plane and line defined by
        //  origin and direction
        scalar normalIntersect(const point& pnt0, const vector& dir) const;

        //- Return the cutting point between the plane and
        // a line passing through the supplied points
        template<class Point, class PointRef>
        scalar lineIntersect(const line<Point, PointRef>& l) const
        {
            return normalIntersect(l.start(), l.vec());
        }

        //- Return the cutting point between this plane and two other planes
        point planePlaneIntersect(const plane&, const plane&) const;

        //- Return the cutting line between this plane and another
        //line planeIntersect(const plane&) const;

    // IOstream Operators

        //- Write plane properties
        friend Ostream& operator<<(Ostream&, const plane&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
