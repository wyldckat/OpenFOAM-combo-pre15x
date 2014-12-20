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

#include "dxFoamExec.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

DXField createDxField
(
    const word& fieldName,
    const DXArray dxArray,
    const DXArray dxPositions,
    const DXArray dxConnections,
    const char* dependency
)
{
    DXField dxField = DXNewField();

    if (!dxField)
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:createDxField : Error creating dxField"
        );
        return DXERROR;
    }

    if (!DXSetStringAttribute((Object)dxField, "name",(char*)fieldName.c_str()))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:createDxField : Error in DXSetStringAttribute name"
        );
        return DXERROR;
    }

    if (!DXSetStringAttribute((Object)dxArray, "dep", (char*)dependency))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:createDxField : Error in DXSetStringAttribute dep"
        );
        return DXERROR;
    }

    if (!DXSetComponentValue(dxField, "data", (Object)dxArray))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:createDxField : Error in DXSetComponentValue data"
        );
        return DXERROR;
    }

    if (!DXSetComponentValue(dxField, "positions", (Object)dxPositions))
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:createDxField : "
            "Error in DXSetComponentValue positions"
        );
        return DXERROR;
    }
    
    if
    (
        dxConnections
     && !DXSetComponentValue(dxField, "connections", (Object)dxConnections)
    )
    {
        DXSetError
        (
            ERROR_INTERNAL,
            "dxFoamExec:createDxField : "
            "Error in DXSetComponentValue connections"
        );
        return DXERROR;
    }

    return dxField;
}


// ************************************************************************* //
