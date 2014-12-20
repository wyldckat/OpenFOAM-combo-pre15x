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

    See README.

\*---------------------------------------------------------------------------*/

#include "dxFvMesh.H"
#include "foamValid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

DXField dxFvMesh::createDxMesh()
{
    DXArray dxPatchNames = NULL;

    if (boundaryMesh().size())
    {
        const char** patchNameStrings =
            new const char*
            [
                boundaryMesh().size()
              + faceZones().size()
            ];

        label nValidPatches = 0;
        forAll (boundaryMesh(), patchi)
        {
            if (valid(boundaryMesh()[patchi]))
            {
                patchNameStrings[nValidPatches++]
                    = boundaryMesh()[patchi].name().c_str();
            }
        }
        forAll(faceZones(), zoneI)
        {
            patchNameStrings[nValidPatches++]
                = faceZones()[zoneI].name().c_str();
        }

        dxPatchNames
            = DXMakeStringListV(nValidPatches, (char**)patchNameStrings);
        delete[] patchNameStrings;

        if (!dxPatchNames)
        {
            DXSetError
            (
                ERROR_INTERNAL,
                "dxFoamExec:dxFvMesh::createDxMesh : Error in DXMakeStringListV"
            );
            return DXERROR;
        }
    }
    else
    {
        Info<< "dxFoamExec:dxFvMesh::createDxMesh : "
            << "no valid patches" << endl;
    }
    
    DXField dxMesh = DXNewField();

    if (!dxMesh)
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:dxFvMesh::createDxMesh : Error creating dxMesh"
        );
        return DXERROR;
    }

    if (!DXSetStringAttribute((Object)dxMesh, "name", "mesh"))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:dxFvMesh::createDxMesh : "
            "Error in DXSetStringAttribute name"
        );
        return DXERROR;
    }

    if (!DXSetComponentValue(dxMesh, "positions", (Object)dxPositions_))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:dxFvMesh::createDxMesh : "
            "Error in DXSetComponentValue positions"
        );
        return DXERROR;
    }
    
    if (!DXSetComponentValue(dxMesh, "connections", (Object)dxConnections_))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:dxFvMesh::createDxMesh : "
            "Error in DXSetComponentValue connections"
        );
        return DXERROR;
    }

    if (!DXSetComponentValue(dxMesh, "BFconnections", (Object)dxBfaceCons_))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:dxFvMesh::createDxMesh : "
            "Error in DXSetComponentValue connections"
        );
        return DXERROR;
    }

    if (!DXSetComponentValue(dxMesh, "patchNames", (Object)dxPatchNames))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:dxFvMesh::createDxMesh : "
            "Error in DXSetComponentValue patchNames"
        );
        return DXERROR;
    }

    if (!DXEndField(dxMesh))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:dxFvMesh::createDxMesh : Error in DXEndField"
        );
        return DXERROR;
    }

    return dxMesh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
