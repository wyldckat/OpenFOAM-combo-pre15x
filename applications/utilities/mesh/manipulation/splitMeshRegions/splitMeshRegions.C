/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Splits mesh into multiple regions and writes them to consecutive
    time directories. Each region is defined as a domain
    whose cells can all be reached by cell-face-cell walking.
    Uses MeshWave.
    Could work in parallel but never tested.

\*---------------------------------------------------------------------------*/

#include "SortableList.H"
#include "argList.H"
#include "regionSplit.H"
#include "fvMeshSubset.H"
#include "IOobjectList.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void subsetVolFields
(
    const fvMeshSubset& subsetter,
    const wordList& fieldNames,
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& subFields
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting field " << fieldName << endl;

        GeometricField<Type, fvPatchField, volMesh> volField    
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, subsetter.interpolate(volField));
    }
}


template<class Type>
void subsetSurfaceFields
(
    const fvMeshSubset& subsetter,
    const wordList& fieldNames,
    PtrList<GeometricField<Type, fvPatchField, surfaceMesh> >& subFields
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting field " << fieldName << endl;

        GeometricField<Type, fvPatchField, surfaceMesh> volField    
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, subsetter.interpolate(volField));
    }
}


bool writeSubMesh(const fvMeshSubset& subsetter)
{
    bool writeOk = subsetter.subMesh().write();

    writeOk = writeOk && labelIOList
    (
        IOobject
        (
            "pointMap",
            subsetter.subMesh().facesInstance(),
            polyMesh::meshSubDir,
            subsetter.subMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        subsetter.pointMap()
    ).write();

    writeOk = writeOk && labelIOList
    (
        IOobject
        (
            "faceMap",
            subsetter.subMesh().facesInstance(),
            polyMesh::meshSubDir,
            subsetter.subMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        subsetter.faceMap()
    ).write();

    writeOk = writeOk && labelIOList
    (
        IOobject
        (
            "cellMap",
            subsetter.subMesh().facesInstance(),
            polyMesh::meshSubDir,
            subsetter.subMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        subsetter.cellMap()
    ).write();

    return writeOk;
}


// Subset mesh and fields. 
void subsetMeshAndFields
(
    Time& runTime,
    fvMeshSubset& subsetter,
    const labelList& cellRegion,
    const label regionI
)
{
    subsetter.setLargeCellSubset(cellRegion, regionI, -1);

    Info<< "Mesh subset in = "
        << runTime.cpuTimeIncrement() << " s\n" << endl;


    IOobjectList objects(subsetter.baseMesh(), runTime.timeName());

    // Read vol fields and subset.

    wordList scalarNames(objects.names(volScalarField::typeName));
    PtrList<volScalarField> scalarFlds(scalarNames.size());
    subsetVolFields(subsetter, scalarNames, scalarFlds);

    wordList vectorNames(objects.names(volVectorField::typeName));
    PtrList<volVectorField> vectorFlds(vectorNames.size());
    subsetVolFields(subsetter, vectorNames, vectorFlds);

    wordList sphericalTensorNames
    (
        objects.names(volSphericalTensorField::typeName)
    );
    PtrList<volSphericalTensorField> sphericalTensorFlds
    (
        sphericalTensorNames.size()
    );
    subsetVolFields(subsetter, sphericalTensorNames, sphericalTensorFlds);

    wordList symmTensorNames(objects.names(volSymmTensorField::typeName));
    PtrList<volSymmTensorField> symmTensorFlds(symmTensorNames.size());
    subsetVolFields(subsetter, symmTensorNames, symmTensorFlds);

    wordList tensorNames(objects.names(volTensorField::typeName));
    PtrList<volTensorField> tensorFlds(tensorNames.size());
    subsetVolFields(subsetter, tensorNames, tensorFlds);

    // Read surface fields and subset.

    wordList surfScalarNames(objects.names(surfaceScalarField::typeName));
    PtrList<surfaceScalarField> surfScalarFlds(surfScalarNames.size());
    subsetSurfaceFields(subsetter, surfScalarNames, surfScalarFlds);

    wordList surfVectorNames(objects.names(surfaceVectorField::typeName));
    PtrList<surfaceVectorField> surfVectorFlds(surfVectorNames.size());
    subsetSurfaceFields(subsetter, surfVectorNames, surfVectorFlds);

    wordList surfSphericalTensorNames
    (
        objects.names(surfaceSphericalTensorField::typeName)
    );
    PtrList<surfaceSphericalTensorField> surfSphericalTensorFlds
    (
        surfSphericalTensorNames.size()
    );
    subsetSurfaceFields
    (
        subsetter,
        surfSphericalTensorNames,
        surfSphericalTensorFlds
    );

    wordList surfSymmTensorNames
    (
        objects.names(surfaceSymmTensorField::typeName)
    );
    PtrList<surfaceSymmTensorField> surfSymmTensorFlds
    (
        surfSymmTensorNames.size()
    );
    subsetSurfaceFields(subsetter, surfSymmTensorNames, surfSymmTensorFlds);

    wordList surfTensorNames(objects.names(surfaceTensorField::typeName));
    PtrList<surfaceTensorField> surfTensorFlds(surfTensorNames.size());
    subsetSurfaceFields(subsetter, surfTensorNames, surfTensorFlds);

    runTime++;

    Info<< "Writing region " << regionI << " mesh to time "
        << runTime.timeName() << endl;

    // Subsetting adds 'subset' prefix. Rename field to be like original.    
    forAll(scalarFlds, i)
    {
        scalarFlds[i].rename(scalarNames[i]);
        scalarFlds[i].writeOpt() = IOobject::AUTO_WRITE;
    }
    forAll(vectorFlds, i)
    {
        vectorFlds[i].rename(vectorNames[i]);
        vectorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
    }
    forAll(tensorFlds, i)
    {
        tensorFlds[i].rename(tensorNames[i]);
        tensorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
    }

    forAll(surfScalarFlds, i)
    {
        surfScalarFlds[i].rename(surfScalarNames[i]);
        surfScalarFlds[i].writeOpt() = IOobject::AUTO_WRITE;
    }
    forAll(surfVectorFlds, i)
    {
        surfVectorFlds[i].rename(surfVectorNames[i]);
        surfVectorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
    }
    forAll(surfTensorNames, i)
    {
        surfTensorFlds[i].rename(surfTensorNames[i]);
        surfTensorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
    }

    writeSubMesh(subsetter);

    Info<< endl;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("largestOnly", "");
    argList::validOptions.insert("insidePoint", "point");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    bool largestOnly = args.options().found("largestOnly");
    bool insidePoint = args.options().found("insidePoint");

    if (insidePoint && largestOnly)
    {
        FatalErrorIn(args.executable())
            << "You cannot specify both -largestOnly"
            << " (keep region with most cells)"
            << " and -insidePoint (keep region containing point)"
            << exit(FatalError);
    }

    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    // Create mesh subsetting engine
    fvMeshSubset subsetter(mesh);


    // Determine connected regions. regionSplit is the labelList with the
    // region per cell.
    regionSplit cellRegion(mesh);

    Info<< endl << "Number of regions:" << cellRegion.nRegions() << nl << endl;

    if (cellRegion.nRegions() == 1)
    {
        Info<< "Only one region. Doing nothing." << endl;
    }
    else if (insidePoint)
    {
        point insidePoint(IStringStream(args.options()["insidePoint"])());

        Pout<< "Keeping region containing point " << insidePoint << endl;

        label regionI = -1;

        label cellI = mesh.findCell(insidePoint);

        Pout<< "Found point " << insidePoint << " in cell " << cellI << endl;

        if (cellI != -1)
        {
            regionI = cellRegion[cellI];
        }

        reduce(regionI, maxOp<label>());

        Pout<< "Subsetting region " << regionI << endl;

        if (regionI == -1)
        {
            FatalErrorIn(args.executable())
                << "Point " << insidePoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh.globalData().bb()
                << exit(FatalError);
        }

        // Subset. Create "oldInternalFaces" patch.
        subsetMeshAndFields
        (
            runTime,
            subsetter,
            cellRegion,
            regionI
        );
    }
    else
    {
        SortableList<label> regionSizes(cellRegion.nRegions(), 0);

        forAll(cellRegion, cellI)
        {
            regionSizes[cellRegion[cellI]]++;
        }
        forAll(regionSizes, regionI)
        {
            reduce(regionSizes[regionI], sumOp<label>());
        }

        // So now all processors have same data. Sort in ascending order
        // so largest is last. Walk through in reverse order so largest is
        // first (cannot reverse list since would not reverse
        // SortableList::indices)

        regionSizes.sort();

        label regionI = 0;

        for (label i = cellRegion.nRegions()-1; i >= 0; i--)
        {
            label oldRegionI = regionSizes.indices()[i];

            Info<< "Region " << regionI << nl
                << "-------- " << endl;

            Info<< "Region " << regionI
                << " mesh has " << regionSizes[i]
                << " cells" << endl;

            // Subset. Create "oldInternalFaces" patch.
            subsetMeshAndFields
            (
                runTime,
                subsetter,
                cellRegion,
                oldRegionI
            );

            if (largestOnly)
            {
                break;
            }
            regionI++;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
