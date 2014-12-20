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

\*---------------------------------------------------------------------------*/

#include "writeLagrangian.H"
#include "OFstream.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "floatScalar.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writeLagrangian
(
    const bool binary,
    const vtkMesh& vMesh,
    const fileName& lagrFileName,       // where to create files
    const wordList& scalarNames,
    const wordList& vectorNames
)
{
    std::ofstream pStream(lagrFileName.c_str());

    pStream
        << "# vtk DataFile Version 2.0" << std::endl
        << "lagrangian" << std::endl;
    if (binary)
    {
        pStream << "BINARY" << std::endl;
    }
    else
    {
        pStream << "ASCII" << std::endl;
    }
    pStream << "DATASET POLYDATA" << std::endl;

    // Read particles
    const fvMesh& mesh = vMesh.mesh();

    Cloud<passiveParticle> parcels(mesh);

    pStream << "POINTS " << parcels.size() << " float" << std::endl;

    DynamicList<floatScalar> partField(3*parcels.size());

    for
    (
        Cloud<passiveParticle>::const_iterator elmnt = parcels.begin();
        elmnt != parcels.end();
        ++elmnt
    )
    {
        writeFuns::insert(elmnt().position(), partField);
    }
    writeFuns::write(pStream, binary, partField);

    pStream
        << "POINT_DATA " << parcels.size() << std::endl
        << "FIELD attributes "
        << scalarNames.size() + vectorNames.size() << std::endl;

    forAll(scalarNames, i)
    {
        IOobject header
        (
            scalarNames[i],
            mesh.time().timeName(),
            "lagrangian",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        IOField<scalar> fld(header);

        pStream << scalarNames[i] << " 1 "
            << parcels.size() << " float" << std::endl;

        DynamicList<floatScalar> fField(parcels.size());

        writeFuns::insert(fld, fField);

        writeFuns::write(pStream, binary, fField);
    }
    forAll(vectorNames, i)
    {
        IOobject header
        (
            vectorNames[i],
            mesh.time().timeName(),
            "lagrangian",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        IOField<vector> fld(header);

        pStream
            << vectorNames[i] << " 3 "
            << parcels.size() << " float" << std::endl;

        DynamicList<floatScalar> fField(3*parcels.size());

        writeFuns::insert(fld, fField);

        writeFuns::write(pStream, binary, fField);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
