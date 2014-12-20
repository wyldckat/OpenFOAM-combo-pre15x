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

Description
    Geometric class that creates a 2D plane and can return the cutting point 
    between a line and the plane.

\*---------------------------------------------------------------------------*/

#include "plane.H"
#include "tensor.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate base point and unit normal vector from plane equation
void plane::calcPntAndVec(const scalarList& C)
{
    if(mag(C[0]) > VSMALL)
    {
        basePoint_ = vector((-C[3]/C[0]), 0, 0);
    }
    else
    {
        if(mag(C[1]) > VSMALL)
        {
            basePoint_ = vector(0, (-C[3]/C[1]), 0);
        }
        else
        {
            if(mag(C[2]) > VSMALL)
            {
                basePoint_ = vector(0, 0, (-C[3]/C[2]));
            }
            else
            {
                FatalErrorIn("void plane::calcPntAndVec(const scalarList& C)")
                    << "At least one plane coefficient must have a value"
                    << abort(FatalError);
            }
        }
    }

    unitVector_ = vector(C[0], C[1], C[2]);
    scalar magUnitVector(mag(unitVector_));

    if (magUnitVector < VSMALL)
    {
        FatalErrorIn("void plane::calcPntAndVec(const scalarList& C)")
            << "Plane normal defined with zero length"
            << abort(FatalError);
    }

    unitVector_ /= magUnitVector;
}


void plane::calcPntAndVec
(
    const point& point1,
    const point& point2,
    const point& point3
)
{
    basePoint_ = (point1 + point2 + point3)/3;
    vector line12 = point1 - point2;
    vector line23 = point2 - point3;

    if
    (
        mag(line12) < VSMALL
     || mag(line23) < VSMALL
     || mag(point3-point1) < VSMALL
    )
    {
        FatalErrorIn
        (
            "void plane::calcPntAndVec\n"
            "(\n"
            "    const point point1,\n"
            "    const point point2,\n"
            "    const point point3\n"
            ")\n"
        ) << "Bad points." << abort(FatalError);
    }

    unitVector_ = line12 ^ line23;
    scalar magUnitVector(mag(unitVector_));

    if (magUnitVector < VSMALL)
    {
        FatalErrorIn
        (
            "void plane::calcPntAndVec\n"
            "(\n"
            "    const point point1,\n"
            "    const point point2,\n"
            "    const point point3\n"
            ")\n"
        )   <<  "Plane normal defined with zero length"
            << abort(FatalError);
    }

    unitVector_ /= magUnitVector;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
plane::plane()
:
    unitVector_(GREAT, GREAT, GREAT),
    basePoint_(GREAT, GREAT, GREAT)
{}


// Construct from plane equation
plane::plane(const scalarList& C) 
{
    calcPntAndVec(C);
}


// Construct from point and normal vector
plane::plane(const point& basePoint, const vector& normalVector)
:
    unitVector_(normalVector),
    basePoint_(basePoint)
{
    scalar magUnitVector(mag(unitVector_));

    if (magUnitVector > VSMALL)
    {
        unitVector_ /= magUnitVector;
    }
    else
    {
        FatalErrorIn
        (
            "plane::plane(const point& basePoint, const vector& normalVector)"
        )   
        << "plane normal has got zero length"
        << abort(FatalError);
    }
}


// Construct from three points
plane::plane
(
    const point& a,
    const point& b,
    const point& c
)
{
    calcPntAndVec(a, b, c);
}


// Construct from dictionary
plane::plane(const dictionary& planeToInletDict)
:
    unitVector_(0, 0, 0),
    basePoint_(0, 0, 0)
{
    word planeTypeDesignator
    (
        planeToInletDict.lookup("planeType")
    );

    if(planeTypeDesignator == word("planeEquation"))
    {
        scalarList C(4);

        dictionary planeEquationCoeffs
        (
            planeToInletDict.subDict("planeEquationDict")
        );

        C[0] = readScalar(planeEquationCoeffs.lookup("a"));
        C[1] = readScalar(planeEquationCoeffs.lookup("b"));
        C[2] = readScalar(planeEquationCoeffs.lookup("c"));
        C[3] = readScalar(planeEquationCoeffs.lookup("d"));

        calcPntAndVec(C);

    }
    else
    {
        if(planeTypeDesignator == word("embeddedPoints"))
        {
            dictionary embeddedPoints
            (
                planeToInletDict.subDict("embeddedPointsDict")
            );

            point point1
            (
                embeddedPoints.lookup("point1")
            );
            point point2
            (
                embeddedPoints.lookup("point2")
            );
            point point3
            (
                embeddedPoints.lookup("point3")
            );

            calcPntAndVec(point1, point2, point3);
        }
        else
        {
            if(planeTypeDesignator == word("pointAndNormal"))
            {

                dictionary pointAndVector
                (
                    planeToInletDict.subDict("pointAndNormalDict")
                );

                basePoint_ = pointAndVector.lookup("basePoint");
                unitVector_ = pointAndVector.lookup("normalVector");
                unitVector_ /= mag(unitVector_);
            }
            else
            {
                FatalIOErrorIn
                (
                    "plane::plane(IOdictionary& planeToInletDict)",
                    planeToInletDict
                )   
                << "Invalid plane type: " << planeTypeDesignator
                << abort(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return plane normal vector
const vector& plane::normal() const
{ 
    return unitVector_;
}


// Return plane base point
const point& plane::refPoint() const
{ 
    return basePoint_;
}


// Return coefficcients for plane equation: ax + by + cz + d = 0
scalarList plane::planeCoeffs() const
{
    scalarList C(4);

    if(unitVector_.component(vector::X))
    {
        C[0] = 1;
        C[1] = unitVector_.component(vector::Y)
                /unitVector_.component(vector::X);
        C[2] = unitVector_.component(vector::Z)
                /unitVector_.component(vector::X);
    }
    else
    {
        if(unitVector_.component(vector::Y))
        {
            C[0] = 0;
            C[1] = 1;
            C[2] = unitVector_.component(vector::Z)
                   /unitVector_.component(vector::Y);
            
        }
        else
        {
            C[0] = 0;
            C[1] = 0;
            C[2] = 1;
        }
    }

    C[3] = - C[0] * basePoint_.component(point::X)
           - C[1] * basePoint_.component(point::Y)
           - C[2] * basePoint_.component(point::Z);

    return C;
}


// Return nearest point in the plane for the given point
point plane::nearestPoint(const point& p) const
{
    return p - unitVector_*((p - basePoint_) & unitVector_);
}


// Return distance from the given point to the plane
scalar plane::distance(const point& p) const
{
    return mag((p - basePoint_) & unitVector_);
}


// Cutting point for plane and line defined by origin and direction
scalar plane::normalIntersect(const point& pnt0, const vector& dir) const
{
    scalar denom = stabilise((dir & unitVector_), VSMALL);

    return ((basePoint_ - pnt0) & unitVector_)/denom;
}


// Cutting point of three planes
point plane::planePlaneIntersect(const plane& plane2, const plane& plane3) const
{
    List<scalarList> pcs(3);
    pcs[0]= planeCoeffs();
    pcs[1]= plane2.planeCoeffs();
    pcs[2]= plane3.planeCoeffs();
    
    tensor a
    (
        pcs[0][0],pcs[0][1],pcs[0][2],
        pcs[1][0],pcs[1][1],pcs[1][2],
        pcs[2][0],pcs[2][1],pcs[2][2]
    );
    
    vector b(pcs[0][3],pcs[1][3],pcs[2][3]);
    
    return (inv(a) & (-b));
    
}


Ostream& operator<<(Ostream& os, const plane& a)
{
    os  << "Plane definition:  normal = " << a.unitVector_
        << " base point = " << a.basePoint_ << endl;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
