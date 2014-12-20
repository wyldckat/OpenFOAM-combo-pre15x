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
    Select cells based on field. Give min and max to select.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cellSet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void createSet
(
    const volScalarField& field,
    const scalar minVal,
    const scalar maxVal
)
{
    cellSet cells(field.mesh(), "subset", field.mesh().nCells()/100 + 1);

    Info<< "Field min:" << min(field) << " max:" << max(field) << endl;
    Info<< "Selecting cells within range " << minVal
        << " .. " << maxVal << endl;

    forAll(field, cellI)
    {
        if (field[cellI] >= minVal && field[cellI] <= maxVal)
        {
            cells.insert(cellI);
        }
    }

    Info<< "Writing " << cells.size() << " selected cells to "
        << cells.localPath(field.mesh(), cells.name())
        << endl;

    cells.write();
}
    


// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("fieldName");
    argList::validArgs.append("(min max)");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    fileName fieldName(args.args()[3]);
    scalarField bounds(IStringStream(args.args()[4])());

    if (bounds.size() != 2)
    {
        FatalErrorIn(args.executable())
            << "Bounds should be two numbers, min and max" << exit(FatalError);
    }

    Info<< "Reading field " << fieldName << endl << endl;


    IOobject fieldObject
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!fieldObject.headerOk())
    {
        FatalErrorIn(args.executable())
            << "Cannot read field " << fieldName << exit(FatalError);
    }

    if (fieldObject.headerClassName() == volScalarField::typeName)
    {
        volScalarField field(fieldObject, mesh);

        createSet(field, bounds[0], bounds[1]);
    }
    else if (fieldObject.headerClassName() == volVectorField::typeName)
    {
        volVectorField field(fieldObject, mesh);

        createSet(Foam::mag(field), bounds[0], bounds[1]);
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Cannot handle fields of type " << fieldObject.headerClassName()
            << exit(FatalError);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
