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

\*---------------------------------------------------------------------------*/

#include "vtkFoam.H"

#include "Time.H"
#include "polyBoundaryMeshEntries.H"
#include "IOobjectList.H"
#include "wordList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"

#include "vtkFoamReader.h"
#include "vtkDataArraySelection.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkCharArray.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::vtkFoam, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "vtkFoamConvertFields.H"

void Foam::vtkFoam::SetName
(
    vtkUnstructuredGrid* vtkMesh,
    const char* name
)
{
    vtkCharArray* nmArray =  vtkCharArray::New();
    nmArray->SetName("Name");
    size_t len = strlen(name);
    nmArray->SetNumberOfTuples(static_cast<vtkIdType>(len)+1);
    char* copy = nmArray->GetPointer(0);
    memcpy(copy, name, len);
    copy[len] = '\0';
    vtkMesh->GetFieldData()->AddArray(nmArray);
    nmArray->Delete();
}


Foam::string Foam::vtkFoam::padTimeString(const string& ts)
{
    return ts + string("            ", max(label(12 - ts.size()), 0));
}


void Foam::vtkFoam::setSelectedTime
(
    Time& runTime,
    vtkFoamReader* reader
)
{
    // Get times list
    instantList Times = runTime.times();
    int timeIndex = min(max(reader->GetTimeStep() + 1, 0), Times.size()-1);

    // If this is the first call timeIndex will be 0 ("constant")
    // so reset to the first time step if one exists and deselect every
    // element of the selection array
    if (timeIndex == 0)
    {
        timeIndex = min(1, Times.size()-1);
        reader->GetTimeSelection()->DisableAllArrays();
    }

    label selectedTimeIndex = -1;
    label nSelectedTimes = reader->GetTimeSelection()->GetNumberOfArrays();

    for (label i=nSelectedTimes-1; i>=0; i--)
    {
        if(reader->GetTimeSelection()->GetArraySetting(i))
        {
            word timeName = reader->GetTimeSelection()->GetArrayName(i);
            forAll(Times, j)
            {
                if (Times[j].name() == timeName)
                {
                    selectedTimeIndex = j;
                    break;
                }
            }
            break;
        }
    }

    if (selectedTimeIndex != -1)
    {
        timeIndex = min(selectedTimeIndex, Times.size()-1);
    }

    if (debug)
    {
        Info<< "Selecting time " << Times[timeIndex].name() << endl;
    }

    runTime.setTime(Times[timeIndex], timeIndex);

    Times = runTime.times();

    reader->SetTimeStepRange(0, max(Times.size()-2, 0));

    // reset the time steps ...
    reader->GetTimeSelection()->RemoveAllArrays();

    int* TimeStepLimits = reader->GetTimeStepLimits();
    label maxStartTimes = min(Times.size(), TimeStepLimits[0]);
    label maxNTimes = min(Times.size() - maxStartTimes, TimeStepLimits[1]);

    for (label i=0; i<maxStartTimes; i++)
    {
        reader->GetTimeSelection()
            ->AddArray(padTimeString(Times[i].name()).c_str());
    }

    if (Times.size() > TimeStepLimits[0] + TimeStepLimits[1])
    {
        reader->GetTimeSelection()->AddArray(padTimeString("...").c_str());
    }

    for (label i=Times.size() - maxNTimes; i<Times.size(); i++)
    {
        reader->GetTimeSelection()
            ->AddArray(padTimeString(Times[i].name()).c_str());
    }

    // Disable all the time selections (which are all selected by default) ...
    reader->GetTimeSelection()->DisableAllArrays();

    // But maintain the selections made previously
    if (selectedTimeIndex != -1 && selectedTimeIndex < Times.size())
    {
        reader->GetTimeSelection()->EnableArray
            (padTimeString(Times[selectedTimeIndex].name()).c_str());
    }
}


void Foam::vtkFoam::updateSelectedRegions()
{
    if (debug)
    {
        Info<< "Foam::vtkFoam::updateSelectedRegions()" << endl;
    }

    label nRegions = reader_->GetRegionSelection()->GetNumberOfArrays();

    // Read the selected patches and add to the region list
    for (int i=0; i<nRegions; i++)
    {
        selectedRegions_[i] = 
            reader_->GetRegionSelection()->GetArraySetting(i);
    }
}


void Foam::vtkFoam::convertMesh()
{
    if (debug)
    {
        Info<< "Foam::vtkFoam::convertMesh()" << endl;
    }

    const fvMesh& mesh = *meshPtr_;

    label nRegions = reader_->GetRegionSelection()->GetNumberOfArrays();

    // Read the internal mesh as region 0 if selected
    if (reader_->GetRegionSelection()->GetArraySetting(0))
    {
        selectedRegions_[0] = true;
        addInternalMesh
        (
            mesh,
            vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(0))
        );
    }
    else
    {
        selectedRegions_[0] = false;

        vtkUnstructuredGrid *vtkMesh =
            vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(0));

        vtkMesh->Initialize();
        SetName(vtkMesh, "(Internal Mesh)");
    }


    // Read the selected patches and add to the region list
    for (label i=1; i<nRegions; i++)
    {
        vtkUnstructuredGrid *vtkMesh =
            vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(i));

        label patchi = i - 1;

        if (reader_->GetRegionSelection()->GetArraySetting(i))
        {
            if (mesh.boundaryMesh()[patchi].size())
            {
                selectedRegions_[i] = true;
                addPatch(mesh.boundaryMesh()[patchi], vtkMesh);
            }
            else
            {
                selectedRegions_[i] = false;
            }
        }
        else
        {
            selectedRegions_[i] = false;
            vtkMesh->Initialize();
            SetName
            (
                vtkMesh,
                ('(' + mesh.boundaryMesh()[patchi].name() + ')').c_str()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkFoam::vtkFoam(const char* const FileName, vtkFoamReader* reader)
:
    reader_(reader),
    dbPtr_(NULL),
    meshPtr_(NULL)
{
    rootPath_ = fileName(FileName).path();
    casePath_ = rootPath_.name();
    rootPath_ = rootPath_.path();

    if (!dir(rootPath_/casePath_))
    {
        return;
    }

    dbPtr_ = new Time(Time::controlDictName, rootPath_, casePath_);
    setSelectedTime(*dbPtr_, reader_);

    if (debug)
    {
        Info<< "vtkFoam::ExecuteInformation: Initialising outputs" << endl;
    }

    reader_->GetRegionSelection()->AddArray("Internal Mesh");

    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
    reader_->SetNthOutput(0, ugrid);
    ugrid->Delete();
    reader_->GetOutput(0)->Initialize();

    polyBoundaryMeshEntries patchEntries
    (
        IOobject
        (
            "boundary",
            dbPtr_->findInstance(polyMesh::meshSubDir, "boundary"),
            polyMesh::meshSubDir,
            *dbPtr_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(patchEntries, entryi)
    {
        reader_->GetRegionSelection()->AddArray
        (
            patchEntries[entryi].keyword().c_str()
        );

        vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
        reader_->SetNthOutput(entryi + 1, ugrid);
        ugrid->Delete();
        reader_->GetOutput(entryi + 1)->Initialize();
    }

    selectedRegions_.setSize(patchEntries.size() + 1);
    selectedRegions_ = true;

    UpdateInformation();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkFoam::~vtkFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "vtkFoamAddFields.H"

void Foam::vtkFoam::UpdateInformation()
{
    if (debug)
    {
        Info<< "TimeStep = " << reader_->GetTimeStep() << endl;
    }

    setSelectedTime(*dbPtr_, reader_);

    // Search for list of objects for this time
    IOobjectList objects(*dbPtr_, dbPtr_->timeName());

    addFields<volScalarField>(reader_->GetVolFieldSelection(), objects);
    addFields<volVectorField>(reader_->GetVolFieldSelection(), objects);
    addFields<volTensorField>(reader_->GetVolFieldSelection(), objects);

    addFields<pointScalarField>(reader_->GetPointFieldSelection(), objects);
    addFields<pointVectorField>(reader_->GetPointFieldSelection(), objects);
    addFields<pointTensorField>(reader_->GetPointFieldSelection(), objects);
}


void Foam::vtkFoam::Update()
{
    // Clear the current set of selected fields

    for (label i=0; i<reader_->GetNumberOfOutputs(); i++)
    {
        vtkUnstructuredGrid *vtkMesh =
            vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(i));

        vtkCellData* cellData = vtkMesh->GetCellData();
        int numberOfCellArrays = cellData->GetNumberOfArrays();

        wordList cellFieldNames(numberOfCellArrays);
        for (int j=0; j<numberOfCellArrays; j++)
        {
            cellFieldNames[j] = cellData->GetArrayName(j);
        }

        for (int j=0; j<numberOfCellArrays; j++)
        {
            cellData->RemoveArray(cellFieldNames[j].c_str());
        }

        vtkPointData* pointData = vtkMesh->GetPointData();
        int numberOfPointArrays = pointData->GetNumberOfArrays();

        wordList pointFieldNames(numberOfPointArrays);
        for (int j=0; j<numberOfPointArrays; j++)
        {
            pointFieldNames[j] = pointData->GetArrayName(j);
        }

        for (int j=0; j<numberOfPointArrays; j++)
        {
            pointData->RemoveArray(pointFieldNames[j].c_str());
        }
    }

    // Check to see if the mesh has been created

    if (!meshPtr_)
    {
        if (debug)
        {
            Info<< "Reading Mesh" << endl;
        }
        meshPtr_ = new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                dbPtr_->timeName(),
                *dbPtr_
            )
        );
        convertMesh();
    }
    else
    {
        boolList oldSelectedRegions = selectedRegions_;
        updateSelectedRegions();
        if
        (
            meshPtr_->readUpdate() != fvMesh::UNCHANGED
         || oldSelectedRegions != selectedRegions_
        )
        {
            convertMesh();
        }
    }

    if (debug)
    {
        Info<< "converting fields" << endl;
    }

    const fvMesh& mesh = *meshPtr_;

    // Construct interpolation on the raw mesh
    Foam::pointMesh pMesh(mesh);

    Foam::volPointInterpolation pInterp(mesh, pMesh);

    // Search for list of objects for this time
    Foam::IOobjectList objects(mesh, dbPtr_->timeName());

    convertVolFields<Foam::scalar>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection()
    );
    convertVolFields<Foam::vector>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection()
    );
    convertVolFields<Foam::tensor>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection()
    );

    convertPointFields<Foam::scalar>
    (
        mesh, objects, reader_->GetPointFieldSelection()
    );
    convertPointFields<Foam::vector>
    (
        mesh, objects, reader_->GetPointFieldSelection()
    );
    convertPointFields<Foam::tensor>
    (
        mesh, objects, reader_->GetPointFieldSelection()
    );

    if (!reader_->GetCacheMesh())
    {
        delete meshPtr_;
        meshPtr_ = NULL;
    }

    if (debug)
    {
        Info<< "done" << endl;
    }
}


// ************************************************************************* //
