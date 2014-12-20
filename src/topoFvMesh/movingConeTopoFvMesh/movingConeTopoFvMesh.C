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

Description

\*---------------------------------------------------------------------------*/

#include "movingConeTopoFvMesh.H"
#include "mapPolyMesh.H"
#include "layerAdditionRemoval.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingConeTopoFvMesh, 0);

    addToRunTimeSelectionTable(topoFvMesh, movingConeTopoFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::movingConeTopoFvMesh::vertexMarkup
(
    const pointField& p,
    const scalar& curLeft,
    const scalar& curRight
) const
{
    Info<< "Updating vertex markup.  curLeft: "
        << curLeft << " curRight: " << curRight << endl;

    tmp<scalarField> tvertexMarkup(new scalarField(p.size()));
    scalarField& vertexMarkup = tvertexMarkup();

    forAll (p, pI)
    {
        if (p[pI].x() < curLeft - SMALL)
        {
            vertexMarkup[pI] = -1;
        }
        else if (p[pI].x() > curRight + SMALL)
        {
            vertexMarkup[pI] = 1;
        }
        else
        {
            vertexMarkup[pI] = 0;
        }
    }

    return tvertexMarkup;
}


void Foam::movingConeTopoFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
     || morphEngine().size() > 0
    )
    {
        Info<< "void movingConeTopoFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    const vectorField& fc = faceCentres();
    const vectorField& fa = faceAreas();

    labelList zone1(fc.size());
    boolList flipZone1(fc.size(), false);
    label nZoneFaces1 = 0;

    labelList zone2(fc.size());
    boolList flipZone2(fc.size(), false);
    label nZoneFaces2 = 0;

    forAll (fc, faceI)
    {
        if
        (
            fc[faceI].x() > -0.003501
         && fc[faceI].x() < -0.003499
        )
        {
            if ((fa[faceI] & vector(1, 0, 0)) < 0)
            {
                flipZone1[nZoneFaces1] = true;
            }

            zone1[nZoneFaces1] = faceI;
            Info<< "face " << faceI << " for zone 1.  Flip: "
                << flipZone1[nZoneFaces1] << endl;
            nZoneFaces1++;
        }
        else if 
        (
            fc[faceI].x() > -0.00701
         && fc[faceI].x() < -0.00699
        )
        {
            zone2[nZoneFaces2] = faceI;

            if ((fa[faceI] & vector(1, 0, 0)) > 0)
            {
                flipZone2[nZoneFaces2] = true;
            }

            Info << "face " << faceI << " for zone 2.  Flip: "
                << flipZone2[nZoneFaces2] << endl;
            nZoneFaces2++;
        }
    }

    zone1.setSize(nZoneFaces1);
    flipZone1.setSize(nZoneFaces1);

    zone2.setSize(nZoneFaces2);
    flipZone2.setSize(nZoneFaces2);

    Info << "zone: " << zone1 << endl;
    Info << "zone: " << zone2 << endl;

    List<pointZone*> pz(0);
    List<faceZone*> fz(2);
    List<cellZone*> cz(0);

    label nFz = 0;

    fz[nFz] =
        new faceZone
        (
            "rightExtrusionFaces",
            zone1,
            flipZone1,
            nFz,
            faceZones()
        );
    nFz++;

    fz[nFz] =
        new faceZone
        (
            "leftExtrusionFaces",
            zone2,
            flipZone2,
            nFz,
            faceZones()
        );
    nFz++;

    fz.setSize(nFz);

    Info << "Adding mesh zones." << endl;
    addZones(pz, fz, cz);


    // Add layer addition/removal interfaces

    List<polyMeshModifier*> tm(2);
    label nMods = 0;

    tm[nMods] =
        new layerAdditionRemoval
        (
            "right",
            nMods,
            *this,
            "rightExtrusionFaces",
            readScalar
            (
                motionDict_.subDict("left").lookup("minThickness")
            ),
            readScalar
            (
                motionDict_.subDict("left").lookup("maxThickness")
            )
        );
    nMods++;

    tm[nMods] =
        new layerAdditionRemoval
        (
            "left",
            nMods,
            *this,
            "leftExtrusionFaces",
            readScalar
            (
                motionDict_.subDict("left").lookup("minThickness")
            ),
            readScalar
            (
                motionDict_.subDict("left").lookup("maxThickness")
            )
        );
    nMods++;

    Info << "Adding " << nMods << " mesh modifiers" << endl;
    tm.setSize(nMods);
    addTopologyModifiers(tm);

    write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::movingConeTopoFvMesh::movingConeTopoFvMesh(const IOobject& io)
:
    topoFvMesh(io),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "meshMotionDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    motionVelAmplitude_(motionDict_.lookup("motionVelAmplitude")),
    motionVelPeriod_(readScalar(motionDict_.lookup("motionVelPeriod"))),
    curMotionVel_
    (
        motionVelAmplitude_*
        Foam::sin(time().value()*M_PI/motionVelPeriod_)
    ),
    leftEdge_(readScalar(motionDict_.lookup("leftEdge"))),
    curLeft_(readScalar(motionDict_.lookup("leftObstacleEdge"))),
    curRight_(readScalar(motionDict_.lookup("rightObstacleEdge"))),
    motionMask_
    (
        vertexMarkup
        (
            points(),
            curLeft_,
            curRight_
        )
    )
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingConeTopoFvMesh::~movingConeTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingConeTopoFvMesh::moveAndMorph()
{
    updateTopology();

    // Calculate the new point positions depending on whether the
    // topological change has happened or not
    pointField newPoints;

    vector curMotionVel_ =
        motionVelAmplitude_*
        Foam::sin(time().value()*M_PI/motionVelPeriod_); 

    if (morphing())
    {
        Info << "Topology change. Calculating motion points" << endl;

        motionMask_ =
            vertexMarkup
            (
                morphMap().preMotionPoints(),
                curLeft_,
                curRight_
            );

        // Create new points by moving old points but using the
        // pre-motion points at the motion selector for the moving
        // region
        newPoints =
            points()
          + (
                pos(0.5 - mag(motionMask_)) // cells above the body
//               + pos(motionMask_ - 0.5)*      // cells in front of the body
//                 (
//                     points().component(vector::X)/curRight
//                 )
//               + pos(-motionMask_ - 0.5)*      // cells behind the body
//                 (
//                     (
//                         points().component(vector::X)
//                       - leftEdge
//                     )/
//                     (curLeft_ - leftEdge_)
//                 )
            )*curMotionVel_*time().deltaT().value();
    }
    else
    {
        Info << "No topology change" << endl;
        // Set the mesh motion
        newPoints =
            points()
          + (
                pos(0.5 - mag(motionMask_)) // cells above the body
//               + pos(motionMask_ - 0.5)*  // cells in front of the body
//                 (
//                     points().component(vector::X)/curRight
//                 )
//               + pos(-motionMask_ - 0.5)*  // cells behind the body
//                 (
//                     (
//                         points().component(vector::X)
//                       - leftEdge
//                     )/
//                     (curLeft_ - leftEdge_)
//                 )
            )*curMotionVel_*time().deltaT().value();
    }

    curLeft_ += curMotionVel_.x()*time().deltaT().value();
    curRight_ += curMotionVel_.x()*time().deltaT().value();

    // The mesh now contains the cells with zero volume

    Info << "Executing mesh motion" << endl;
    movePoints(newPoints);

    //  The mesh now has got non-zero volume cells
}


// ************************************************************************* //
