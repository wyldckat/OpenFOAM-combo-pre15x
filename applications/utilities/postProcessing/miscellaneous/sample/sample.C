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
    Sample field data with a choice of interpolation schemes, sampling options
    and write formats.

    interpolationScheme : choice of
        cell            : use cell-centre value only; constant over cells
        cellPoint       : use cell-centre and vertex values
        cellPointFace   : use cell-centre, vertex and face values.

    sample: choice of
        uniform             evenly distributed points on line
        face                one point per face intersection
        midPoint            one point per cell, inbetween two face intersections
        midPointAndFace     combination of face and midPoint

        curve               specified points, not nessecary on line, uses
                            tracking
        cloud               specified points, uses findCell

    writeFormat : choice of
        xmgr
        jplot
        gnuplot
        raw

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "argList.H"
#include "OSspecific.H"

#include "Cloud.H"
#include "passiveParticle.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "volPointInterpolation.H"

#include "writer.H"
#include "sampleSet.H"
#include "volFieldSampler.H"
#include "dictionaryEntry.H"

#include "combineSampleSets.H"
#include "combineSampleValues.H"

using namespace Foam;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    //
    // Hack: initialize Cloud to initialize the processor table so from
    // now on we can use cloud on single processors only.
    //
    Cloud<passiveParticle> dummyCloud(mesh, IDLList<passiveParticle>());

    // Read control dictionary
    IOdictionary sampleDict
    (
        IOobject
        (
            "sampleDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word interpolationScheme(sampleDict.lookup("interpolationScheme"));
    const wordList fieldNames = sampleDict.lookup("fields");

    //
    // Construct writers
    //

    word writeFormat(sampleDict.lookup("writeFormat"));

    autoPtr<writer<scalar> > scalarFormatter
    (
        writer<scalar>::New(writeFormat)
    );
    autoPtr<writer<vector> > vectorFormatter
    (
        writer<vector>::New(writeFormat)
    );
    autoPtr<writer<tensor> > tensorFormatter
    (
        writer<tensor>::New(writeFormat)
    );

    //
    // Construct interpolation dictionary (same interpolation for all fields)
    //

    dictionary interpolationSchemes;

    forAll(fieldNames, fieldI)
    {
        interpolationSchemes.add
        (
            fieldNames[fieldI],
            interpolationScheme
        );
    }

    // Set up interpolation
    autoPtr<pointMesh> pMeshPtr(new pointMesh(mesh));
    autoPtr<volPointInterpolation> pInterpPtr
    (
        new volPointInterpolation(mesh, pMeshPtr())
    );

    // Set up mesh searching
    meshSearch searchEngine(mesh, true);


    fileName samplePath;

    if (Pstream::master())
    {
        if (Pstream::parRun())
        {
            samplePath = runTime.path()/".."/"samples";
        }
        else
        {
            samplePath = runTime.path()/"samples";
        }

        if (exists(samplePath))
        {
            Info<< "Deleting samples/ directory" << endl << endl;
            rmDir(samplePath);
        }
    }

    fileName oldPointsDir("constant");

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        //
        // Handle geometry/topology changes
        //
        polyMesh::readUpdateState state = mesh.readUpdate();

        if
        (
            state == polyMesh::POINTS_MOVED
         || state == polyMesh::TOPO_CHANGE
        )
        {
            // Geometry and topology changes            
            searchEngine.correct();

            pMeshPtr.reset(new pointMesh(mesh));

            pInterpPtr.reset(new volPointInterpolation(mesh, pMeshPtr()));
        }


        //
        // Construct sampling point generators
        //

        ptrList<sampleSet> sampleSets
        (
            sampleDict.lookup("sampleSets"),
            sampleSet::iNew(mesh, searchEngine)
        );
        if (sampleSets.size() < 1)
        {
            FatalErrorIn(args.executable())
                << "No sampleSets provided in sampleDict"
                << exit(FatalError);
        }

        // Storage for interpolated values
        ptrList<volFieldSampler<scalar> > sampledScalarFields
        (
            fieldNames.size()
        );
        ptrList<volFieldSampler<vector> > sampledVectorFields
        (
            fieldNames.size()
        );
        ptrList<volFieldSampler<tensor> > sampledTensorFields
        (
            fieldNames.size()
        );

        //
        // Do actual interpolation
        //

        label nScalarFields = 0;
        label nVectorFields = 0;
        label nTensorFields = 0;

        forAll(fieldNames, fieldI)
        {
            const word& fieldName = fieldNames[fieldI];

            IOobject fieldHeader
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            // Determine number of processor actually having this field
            label fieldFound = (fieldHeader.headerOk() ? 1 : 0);
            reduce(fieldFound, sumOp<label>());

            if (fieldFound == Pstream::nProcs())
            {
                if 
                (
                    fieldHeader.headerClassName() == volScalarField::typeName
                )
                {
                    Info<< "Sampling " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volScalarField sField(fieldHeader, mesh);

                    sampledScalarFields.hook
                    (
                        new volFieldSampler<scalar>
                        (
                            pInterpPtr(),
                            interpolationSchemes,
                            sField,
                            sampleSets
                        )
                    );

                    nScalarFields++;
                }
                else if
                (
                    fieldHeader.headerClassName() == volVectorField::typeName
                )
                {
                    Info<< "Sampling " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volVectorField vField(fieldHeader, mesh);

                    sampledVectorFields.hook
                    (
                        new volFieldSampler<vector>
                        (
                            pInterpPtr(),
                            interpolationSchemes,
                            vField,
                            sampleSets
                        )
                    );

                    nVectorFields++;
                }
                else if
                (
                    fieldHeader.headerClassName() == volTensorField::typeName
                )
                {
                    Info<< "Sampling " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volTensorField tField(fieldHeader, mesh);

                    sampledTensorFields.hook
                    (
                        new volFieldSampler<tensor>
                        (
                            pInterpPtr(),
                            interpolationSchemes,
                            tField,
                            sampleSets
                        )
                    );

                    nTensorFields++;
                }
            }
            else if (fieldFound != 0)
            {
                FatalErrorIn(args.executable())
                    << "Did not find field " << fieldName
                    << " on all processors" << exit(FatalError);
            }
            else if (fieldName.find('(') != string::npos)
            {
                if (fieldName.find("component") != string::npos)
                {
                    string baseFieldName(fieldName(0, fieldName.find('.')));

                    IOobject fieldHeader
                    (
                        baseFieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    );

                    // Determine number of processor actually having this field
                    label fieldFound = (fieldHeader.headerOk() ? 1 : 0);
                    reduce(fieldFound, sumOp<label>());

                    if (fieldFound == Pstream::nProcs())
                    {
                        if
                        (
                            fieldHeader.headerClassName()
                         == volVectorField::typeName
                        )
                        {
                            size_t cmptPos(fieldName.find_last_of("012"));

                            if (cmptPos == string::npos)
                            {
                                FatalErrorIn(args.executable())
                                    << "cannot find component index for "
                                    << fieldName << exit(FatalError);
                            }

                            direction cmpt =
                                atoi(string(fieldName[cmptPos]).c_str());

                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volVectorField vField(fieldHeader, mesh);

                            sampledScalarFields.hook
                            (
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    vField.component(cmpt)(),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else if
                        (
                            fieldHeader.headerClassName()
                         == volTensorField::typeName
                        )
                        {
                            size_t cmptPos(fieldName.find_last_of("012345678"));

                            if (cmptPos == string::npos)
                            {
                                FatalErrorIn(args.executable())
                                    << "cannot find component index for "
                                    << fieldName
                                    << exit(FatalError);
                            }

                            direction cmpt =
                                atoi(string(fieldName[cmptPos]).c_str());

                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volTensorField tField(fieldHeader, mesh);

                            sampledScalarFields.hook
                            (
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    tField.component(cmpt)(),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else
                        {
                            FatalErrorIn(args.executable())
                                << "component function not supported for field "
                                << fieldName << " of type "
                                << fieldHeader.headerClassName()
                                << exit(FatalError);
                        }
                    }
                    else if (fieldFound != 0)
                    {
                        FatalErrorIn(args.executable())
                            << "Did not find field " << baseFieldName
                            << " on all processors" << exit(FatalError);
                    }

                }
                else if (fieldName.find("mag") != string::npos)
                {
                    string baseFieldName
                    (
                        fieldName(fieldName.find('(') + 1,
                        fieldName.find(')') - fieldName.find('(') - 1)
                    );

                    IOobject fieldHeader
                    (
                        baseFieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    );

                    // Determine number of processor actually having this field
                    label fieldFound = (fieldHeader.headerOk() ? 1 : 0);
                    reduce(fieldFound, sumOp<label>());

                    if (fieldFound == Pstream::nProcs())
                    {
                        if
                        (
                            fieldHeader.headerClassName()
                         == volVectorField::typeName
                        )
                        {
                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volVectorField vField(fieldHeader, mesh);

                            sampledScalarFields.hook
                            (
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    mag(vField),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else if
                        (
                            fieldHeader.headerClassName()
                         == volTensorField::typeName
                        )
                        {
                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volTensorField tField(fieldHeader, mesh);

                            sampledScalarFields.hook
                            (
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    mag(tField),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else
                        {
                            FatalErrorIn(args.executable())
                                << "mag function not supported for field "
                                << fieldName << " of type "
                                << fieldHeader.headerClassName()
                                << exit(FatalError);
                        }
                    }
                    else if (fieldFound != 0)
                    {
                        FatalErrorIn(args.executable())
                            << "Did not find field " << baseFieldName
                            << " on all processors" << exit(FatalError);
                    }
                }
            }
        }

        // Set the sampledFields to the correct size
        sampledScalarFields.setSize(nScalarFields);
        sampledVectorFields.setSize(nVectorFields);
        sampledTensorFields.setSize(nTensorFields);

        //
        // Now we have all results
        //   - sampleSets       : list of all sampling sets
        //   - sampledXXXFields : list of all sampled fields
        //

        // Combine sampleSets from processors. Sort by curveDist. Return
        // ordering in indexSets.
        // Note: only master results are valid

        ptrList<coordSet> masterSampleSets(sampleSets.size());
        labelListList indexSets(sampleSets.size());
        combineSampleSets(sampleSets, masterSampleSets, indexSets);


        // Combine sampled fields from processors.
        // Note: only master results are valid

        ptrList<volFieldSampler<scalar> > masterScalarFields
        (
            sampledScalarFields.size()
        );
        combineSampleValues(sampledScalarFields, indexSets, masterScalarFields);

        ptrList<volFieldSampler<vector> > masterVectorFields
        (
            sampledVectorFields.size()
        );
        combineSampleValues(sampledVectorFields, indexSets, masterVectorFields);

        ptrList<volFieldSampler<tensor> > masterTensorFields
        (
            sampledTensorFields.size()
        );
        combineSampleValues(sampledTensorFields, indexSets, masterTensorFields);


        //
        // Write each set, each Field type (scalar/vector/tensor) to separate
        // file.
        //

        if (Pstream::master())
        {
            fileName timeDir(samplePath/runTime.timeName());

            // Mirror the time structure under "samples"
            mkDir(timeDir);

            forAll(masterSampleSets, setI)
            {
                // ScalarFields

                if (masterScalarFields.size() > 0)
                {
                    Foam::HashTable<scalarField*> scalarValueSets;
                    forAll(masterScalarFields, fieldI)
                    {
                        scalarValueSets.insert
                        (
                            masterScalarFields[fieldI].name(),
                            &(masterScalarFields[fieldI][setI])
                        );
                    }

                    fileName scalarFName
                    (
                        timeDir
                      / scalarFormatter().getFileName
                        (
                            masterSampleSets[setI],
                            scalarValueSets
                        )
                    );

                    Info<< "Writing scalarFields to " << scalarFName << endl;

                    scalarFormatter().write
                    (
                        masterSampleSets[setI],
                        scalarValueSets,
                        OFstream(scalarFName)()
                    );
                }

                // VectorFields

                if (masterVectorFields.size() > 0)
                {
                    Foam::HashTable<vectorField*> vectorValueSets;
                    forAll(masterVectorFields, fieldI)
                    {
                        vectorValueSets.insert
                        (
                            masterVectorFields[fieldI].name(),
                            &(masterVectorFields[fieldI][setI])
                        );
                    }

                    fileName vectorFName
                    (
                        timeDir
                      / vectorFormatter().getFileName
                        (
                            masterSampleSets[setI],
                            vectorValueSets
                        )
                    );

                    Info<< "Writing vectorFields to " << vectorFName << endl;

                    vectorFormatter().write
                    (
                        masterSampleSets[setI],
                        vectorValueSets,
                        OFstream(vectorFName)()
                    );
                }

                // TensorFields

                if (masterTensorFields.size() > 0)
                {
                    Foam::HashTable<tensorField*> tensorValueSets;
                    forAll(masterTensorFields, fieldI)
                    {
                        tensorValueSets.insert
                        (
                            masterTensorFields[fieldI].name(),
                            &(masterTensorFields[fieldI][setI])
                        );
                    }

                    fileName tensorFName
                    (
                        timeDir
                      / tensorFormatter().getFileName
                        (
                            masterSampleSets[setI],
                            tensorValueSets
                        )
                    );

                    Info<< "Writing tensorFields to " << tensorFName << endl;

                    tensorFormatter().write
                    (
                        masterSampleSets[setI],
                        tensorValueSets,
                        OFstream(tensorFName)()
                    );
                }
            }

            Info<< endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
