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

#include "dx.H"
#include "fileName.H"
#include "OFstream.H"
#include "faceList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::dx<Type>::writeDXGeometry
(
    const pointField& points,
    const faceList& faces,
    Ostream& os
) const
{
    // Write vertex coordinates

    os  << "# The irregular positions" << endl
        << "object 1 class array type float rank 1 shape 3 items "
        << points.size() << " data follows" << endl;

    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }
    os  << endl;

    // Write triangles

    os  << "# The irregular connections (triangles)" << endl
        << "object 2 class array type int rank 1 shape 3 items "
        << faces.size() << " data follows" << endl;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        if (f.size() != 3)
        {
            FatalErrorIn
            (
                "writeDXGeometry(Ostream&, const pointField&, const faceList&)"
            )   << "Face " << faceI << " vertices " << f
                << " is not a triangle."
                << exit(FatalError);
        }

        os << f[0] << ' ' << f[1] << ' ' << f[2] << endl;
    }
    os << "attribute \"element type\" string \"triangles\"" << endl
       << "attribute \"ref\" string \"positions\"" << endl << endl;
}


// Write scalarField in DX format
template<class Type>
void Foam::dx<Type>::writeDXData
(
    const pointField& points,
    const scalarField& values,
    Ostream& os
) const
{
    // Write data
    os  << "object 3 class array type float rank 0 items "
        << values.size()
        << " data follows" << endl;

    forAll(values, elemI)
    {
        os << values[elemI] << endl;
    }

    if (values.size() == points.size())
    {
        os  << endl << "attribute \"dep\" string \"positions\""
            << endl << endl;
    }
    else
    {
        os  << endl << "attribute \"dep\" string \"connections\""
            << endl << endl;
    }
}


// Write vectorField in DX format
template<class Type>
void Foam::dx<Type>::writeDXData
(
    const pointField& points,
    const vectorField& values,
    Ostream& os
) const
{
    // Write data
    os  << "object 3 class array type float rank 1 shape 3 items "
        << values.size()
        << " data follows" << endl;

    forAll(values, elemI)
    {
        os  << values[elemI].x() << ' ' << values[elemI].y()
            << ' ' << values[elemI].z()
            << endl;
    }

    if (values.size() == points.size())
    {
        os  << endl << "attribute \"dep\" string \"positions\""
            << endl << endl;
    }
    else
    {
        os  << endl << "attribute \"dep\" string \"connections\""
            << endl << endl;
    }
}


// Write tensorField in DX format
template<class Type>
void Foam::dx<Type>::writeDXData
(
    const pointField& points,
    const tensorField& values,
    Ostream& os
) const
{
    // Write data
    os  << "object 3 class array type float rank 2 shape 3 items "
        << values.size()
        << " data follows" << endl;

    forAll(values, elemI)
    {
        const tensor& t = values[elemI];

        os  << t.xx() << ' ' << t.xy() << ' ' << t.xz()
            << t.yx() << ' ' << t.yy() << ' ' << t.yz()
            << t.zx() << ' ' << t.zy() << ' ' << t.zz()
            << endl;
    }

    if (values.size() == points.size())
    {
        os  << endl << "attribute \"dep\" string \"positions\""
            << endl << endl;
    }
    else
    {
        os  << endl << "attribute \"dep\" string \"connections\""
            << endl << endl;
    }
}


// Write trailer in DX format
template<class Type>
void Foam::dx<Type>::writeDXTrailer(Ostream& os) const
{
    os  << "# the field, with three components: \"positions\","
        << " \"connections\", and \"data\"" << endl
        << "object \"irregular positions irregular "
        << "connections\" class field"
        << endl
        << "component \"positions\" value 1" << endl
        << "component \"connections\" value 2" << endl
        << "component \"data\" value 3" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::dx<Type>::dx()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::dx<Type>::~dx()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::dx<Type>::write
(
    const fileName& samplePath,
    const fileName& timeDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const fileName& fieldName,
    const Field<Type>& values
) const
{
    fileName surfaceDir(samplePath/timeDir);

    if (!exists(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    fileName planeFName(surfaceDir/fieldName + '_' + surfaceName + ".dx");

    Info<< "Writing field " << fieldName << " to " << planeFName << endl;

    OFstream dxFile(planeFName);

    writeDXGeometry(points, faces, dxFile);

    writeDXData(points, values, dxFile);

    writeDXTrailer(dxFile);

    dxFile << "end" << endl;
}


// ************************************************************************* //
