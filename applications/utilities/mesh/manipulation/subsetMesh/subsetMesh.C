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
    Selects a section of mesh based on a cellSet.

    The utility sub-sets the mesh to choose only a part of interest. Check
    the setSet/cellSet utilities to see how to select cells based on various.

    The mesh will subset all points, faces and cells needed to make a sub-mesh
    but will not preserve attached boundary types.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshSubset.H"
#include "cellSet.H"
#include "IOobjectList.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void subsetFields
(
    const meshSubset& mesh,
    const wordList& fieldNames,
    ptrList<GeometricField<Type, fvPatchField, volMesh> >& subFields
)
{
    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting field " << fieldName << endl;

        GeometricField<Type, fvPatchField, volMesh> volField    
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        subFields.hook(mesh.interpolate(volField));
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("set");
    argList::validOptions.insert("patch", "patch name");

#   include "setRootCase.H"
#   include "createTime.H"


    word setName(args.args()[3]);


    Info<< "Reading cell set from " << setName << endl << endl;

    Info<< "Reading mesh for time = " << runTime.value() << endl;

    Info<< "Create mesh\n" << endl;
    meshSubset mesh
    (
        IOobject
        (
            meshSubset::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    label patchI = -1;

    if (args.options().found("patch"))
    {
        word patchName(args.options()["patch"]);

        patchI = mesh.boundaryMesh().findPatchID(patchName);

        if (patchI == -1)
        {
            FatalErrorIn(args.executable()) << "Illegal patch " << patchName
                << nl << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        Info<< "Adding exposed internal faces to patch " << patchName << endl
            << endl;
    }
    else
    {
        Info<< "Adding exposed internal faces to a patch called"
            << " \"oldInternalFaces\" (created if nessecary)" << endl
            << endl;
    }


    cellSet currentSet(mesh, setName);

    mesh.setCellSubset(currentSet, patchI);

    IOobjectList objects(mesh, runTime.timeName());

    // Read fields and subset.

    wordList scalarNames(objects.names(volScalarField::typeName));
    ptrList<volScalarField> scalarFlds(scalarNames.size());

    subsetFields(mesh, scalarNames, scalarFlds);

    wordList vectorNames(objects.names(volVectorField::typeName));
    ptrList<volVectorField> vectorFlds(vectorNames.size());

    subsetFields(mesh, vectorNames, vectorFlds);

    wordList tensorNames(objects.names(volTensorField::typeName));
    ptrList<volTensorField> tensorFlds(tensorNames.size());

    subsetFields(mesh, tensorNames, tensorFlds);

    runTime++;

    Info << "Writing polyMesh to time " << runTime.value() << endl;
    mesh.subMesh().write();


    // Subsetting adds 'subset' prefix. Rename field to be like original.    
    forAll(scalarFlds, i)
    {
        scalarFlds[i].rename(scalarNames[i]);

        scalarFlds[i].write();
    }
    forAll(vectorFlds, i)
    {
        vectorFlds[i].rename(vectorNames[i]);

        vectorFlds[i].write();
    }
    forAll(tensorFlds, i)
    {
        tensorFlds[i].rename(tensorNames[i]);

        tensorFlds[i].write();
    }

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
