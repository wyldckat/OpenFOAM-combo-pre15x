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
    Generic Geometric Type field class

\*---------------------------------------------------------------------------*/

#include "GeometricField.H"
#include "Time.H"
#include "demandDrivenData.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// check mesh for two fields

#define checkField(gf1, gf2, op)                                    \
if ((gf1).mesh() != (gf2).mesh())                                   \
{                                                                   \
    FatalErrorIn("checkField(gf1, gf2, op)")                        \
        << "different mesh for fields "                             \
        << (gf1).name() << " and " << (gf2).name()                  \
        << " during operatrion " <<  op                             \
        << abort(FatalError);                                       \
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
bool GeometricField<Type, PatchField, GeoMesh>::readIfPresent()
{
    if (readOpt() == IOobject::MUST_READ)
    {
        Warning
            << "GeometricField<Type, PatchField, GeoMesh>::readIfPresent() : "
            << endl
            << "    read option IOobject::MUST_READ "
            << "suggests that a read constuctor for field " << name()
            << " would be more appropriate." << endl;
    }
    else if (readOpt() == IOobject::READ_IF_PRESENT)
    {
        if (headerOk())
        {
            GeometricField<Type, PatchField, GeoMesh> gfRead(*this, mesh_);
            operator==(gfRead);

            readOldTimeIfPresent();

            return true;
        }
    }

    return false;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool GeometricField<Type, PatchField, GeoMesh>::readOldTimeIfPresent()
{
    // Read the old time field if present
    IOobject field0
    (
        name()  + "_0",
        time().timeName(),
        db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
        
    if (field0.headerOk())
    {
        if (debug)
        {
            Info<< "Reading old time level for field"
                << endl << info() << endl;
        }

        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            field0,
            mesh_
        );

        field0Ptr_->timeIndex_ = timeIndex_ - 1;
        
        if (!field0Ptr_->readOldTimeIfPresent())
        {
            field0Ptr_->oldTime();
        }
        
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Constructor given a GeometricField and dimensionSet
// This allocates storage for the field but not values.
// Note : This constructor should only be used to
//       construct TEMPORARY variables
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const word& patchFieldType
)
:
    regIOobject(io),
    Field<Type>(GeoMesh::size(mesh)),
    mesh_(mesh),
    dimensions_(ds),
    timeIndex_(time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << info() << endl;
    }

    readIfPresent();
}


// Constructor given a GeometricField and dimensionSet
// This allocates storage for the field but not values.
// Note : This constructor should only be used to
//       construct TEMPORARY variables
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const wordList& patchFieldTypes
)
:
    regIOobject(io),
    Field<Type>(GeoMesh::size(mesh)),
    mesh_(mesh),
    dimensions_(ds),
    timeIndex_(time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << info() << endl;
    }

    readIfPresent();
}


// Constructor given a GeometricField and dimensioned<Type>
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
:
    regIOobject(io),
    Field<Type>(GeoMesh::size(mesh), dt.value()),
    mesh_(mesh),
    dimensions_(dt.dimensions()),
    timeIndex_(time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << info() << endl;
    }

    boundaryField_ == dt.value();

    readIfPresent();
}


// Constructor given a GeometricField and dimensioned<Type>
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes
)
:
    regIOobject(io),
    Field<Type>(GeoMesh::size(mesh), dt.value()),
    mesh_(mesh),
    dimensions_(dt.dimensions()),
    timeIndex_(time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << info() << endl;
    }

    boundaryField_ == dt.value();

    readIfPresent();
}


//  construct from components
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& iField,
    const ptrList<PatchField<Type> >& ptfl
)
:
    regIOobject(io),
    Field<Type>(iField),
    mesh_(mesh),
    dimensions_(ds),
    timeIndex_(time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh_.boundary(), *this, ptfl)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing from components"
            << endl << info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<typename GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField>
GeometricField<Type, PatchField, GeoMesh>::readField()
{
    Istream& is = readStream(typeName);

    if (is.version() < 2.0)
    {
        FatalIOErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::readField()",
            is
        )   << "IO versions < 2.0 are not supported."
            << exit(FatalIOError);
    }

    dictionary fieldDict(is);

    dimensions_.reset(dimensionSet(fieldDict.lookup("dimensions")));

    Type fieldAverage = pTraits<Type>::zero;

    if (fieldDict.found("referenceLevel"))
    {
        fieldAverage = pTraits<Type>(fieldDict.lookup("referenceLevel"));
    }

    Field<Type>::operator=
    (
        tmp<Field<Type> >
        (
            new Field<Type>("internalField", fieldDict, GeoMesh::size(mesh_))
        )
    );

    tmp<GeometricBoundaryField> tboundaryField
    (
        new GeometricBoundaryField
        (
            mesh_.boundary(),
            *this,
            fieldDict.subDict("boundaryField")
        )
    );

    GeometricBoundaryField& boundaryField = tboundaryField();

    // Check fieldAverage_, if not zero then switch field referencing on
    if (mag(fieldAverage) > SMALL)
    {
        Field<Type>::operator+=(fieldAverage);

        forAll(boundaryField, patchi)
        {
            boundaryField[patchi] == boundaryField[patchi] + fieldAverage;
        }
    }

    return tboundaryField;
}


template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh
)
:
    regIOobject(io),
    Field<Type>(0),
    mesh_(mesh),
    dimensions_(dimless),
    timeIndex_(time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, readField())
{
    close();

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(mesh_))
    {
        FatalIOErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::GeometricField"
            "(const IOobject&, const Mesh&)",
            readStream(typeName)
        )   << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(mesh_)
            << exit(FatalIOError);
    }

    readOldTimeIfPresent();

    if (debug)
    {
        Info<< "Finishing read-construct of "
               "GeometricField<Type, PatchField, GeoMesh>"
            << endl << info() << endl;
    }
}


// construct as copy
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
#   ifdef ConstructFromTmp
    regIOobject(gf),
#   else
    regIOobject(gf, true),
#   endif
    Field<Type>(gf),
    mesh_(gf.mesh_),
    dimensions_(gf.dimensions_),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy"
            << endl << info() << endl;
    }

    if (gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            *gf.field0Ptr_
        );
    }

    writeOpt() = IOobject::NO_WRITE;
}

// construct as copy of tmp<GeometricField> deleting argument
#ifdef ConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
:
    regIOobject(tgf(), true),
    Field<Type>((Field<Type>&)tgf(), tgf.isTmp()),
    mesh_(tgf().mesh_),
    dimensions_(tgf().dimensions_),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy"
            << endl << info() << endl;
    }

    writeOpt() = IOobject::NO_WRITE;

    tgf.clear();
}
#endif


// construct as copy resetting IO parameters
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    regIOobject(io),
    Field<Type>(gf),
    mesh_(gf.mesh_),
    dimensions_(gf.dimensions_),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy resetting IO params"
            << endl << info() << endl;
    }

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


// construct as copy resetting name
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    regIOobject(IOobject(newName, gf.time().timeName(), gf.db())),
    Field<Type>(gf),
    mesh_(gf.mesh_),
    dimensions_(gf.dimensions_),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy resetting name"
            << endl << info() << endl;
    }

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            newName + "_0",
            *gf.field0Ptr_
        );
    }
}


// construct as copy resetting name
#ifdef ConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
:
    regIOobject(IOobject(newName, tgf().time().timeName(), tgf().db())),
    Field<Type>((Field<Type>&)tgf(), tgf.isTmp()),
    mesh_(tgf().mesh_),
    dimensions_(tgf().dimensions_),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing from tmp resetting name"
            << endl << info() << endl;
    }

    tgf.clear();
}
#endif

// construct as copy resetting IO parameters and patch type
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const word& patchFieldType
)
:
    regIOobject(io),
    Field<Type>(gf),
    mesh_(gf.mesh_),
    dimensions_(gf.dimensions_),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh_.boundary(), *this, patchFieldType)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy resetting IO params"
            << endl << info() << endl;
    }

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


// construct as copy resetting IO parameters and boundary types
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const wordList& patchFieldTypes
)
:
    regIOobject(io),
    Field<Type>(gf),
    mesh_(gf.mesh_),
    dimensions_(gf.dimensions_),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh_.boundary(), *this, patchFieldTypes)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy resetting IO params"
            << endl << info() << endl;
    }

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


// * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * * //

//- Return a pointer to a GeometricField created by re-using the given
//  tmp<GeometricField> and resetting the patch types to the given type
template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> >
GeometricField<Type, PatchField, GeoMesh>::New
(
    const IOobject& io,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf,
    const dimensionSet& ds
)
{
    GeometricField<Type, PatchField, GeoMesh>* gfPtr = tgf.ptr();
    GeometricField<Type, PatchField, GeoMesh>& gf = *gfPtr;

    deleteDemandDrivenData(gf.field0Ptr_);
    deleteDemandDrivenData(gf.fieldPrevIterPtr_);

    gf.regIOobject::operator=(io);
    gf.dimensions().reset(ds);

    return tmp<GeometricField<Type, PatchField, GeoMesh> >(gfPtr);
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>::~GeometricField()
{
    deleteDemandDrivenData(field0Ptr_);
    deleteDemandDrivenData(fieldPrevIterPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return internal field
template<class Type, template<class> class PatchField, class GeoMesh>
Field<Type>& GeometricField<Type, PatchField, GeoMesh>::internalField()
{
    storeOldTimes();
    return *this;
}


// Return reference to GeometricBoundaryField
template<class Type, template<class> class PatchField, class GeoMesh>
typename GeometricField<Type, PatchField, GeoMesh>::
GeometricBoundaryField&
GeometricField<Type, PatchField, GeoMesh>::boundaryField()
{
    storeOldTimes();
    return boundaryField_;
}


// Store old-time field
template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::storeOldTimes() const
{
    if
    (
        field0Ptr_
     && timeIndex_ != time().timeIndex()
     && !(name().size() > 2 && name()(name().size()-2, 2) == "_0")
    )
    {
        storeOldTime();

        // Correct time index
        timeIndex_ = time().timeIndex();
    }
}

// Store old-time field
template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::storeOldTime() const
{
    if (field0Ptr_)
    {
        field0Ptr_->storeOldTime();

        if (debug)
        {
            Info<< "Storing old time field for field" << endl << info() << endl;
        }

        *field0Ptr_ == *this;
        field0Ptr_->timeIndex_ = timeIndex_;

        if (field0Ptr_->field0Ptr_)
        {
            field0Ptr_->writeOpt() = writeOpt();
        }
    }
}

// Return the number of old time fields stored
template<class Type, template<class> class PatchField, class GeoMesh>
label GeometricField<Type, PatchField, GeoMesh>::nOldTimes() const
{
    if (field0Ptr_)
    {
        return field0Ptr_->nOldTimes() + 1;
    }
    else
    {
        return 0;
    }
}

// Return old time internal field
template<class Type, template<class> class PatchField, class GeoMesh>
const GeometricField<Type, PatchField, GeoMesh>&
GeometricField<Type, PatchField, GeoMesh>::oldTime() const
{
    if (!field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                name() + "_0",
                time().timeName(),
                db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this
        );
    }
    else
    {
        storeOldTimes();
    }

    return *field0Ptr_;
}

// Return old time internal field
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>&
GeometricField<Type, PatchField, GeoMesh>::oldTime()
{
    ((const GeometricField<Type, PatchField, GeoMesh>&)(*this)).oldTime();
    return *field0Ptr_;
}


// Store previous iteration field
template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::storePrevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        if (debug)
        {
            Info<< "Allocating previous iteration field" << endl
                << info() << endl;
        }

        fieldPrevIterPtr_ =
            new GeometricField<Type, PatchField, GeoMesh>(*this);
    }
    else
    {
        *fieldPrevIterPtr_ == *this;
    }
}


// Return previous iteration field
template<class Type, template<class> class PatchField, class GeoMesh>
const GeometricField<Type, PatchField, GeoMesh>&
GeometricField<Type, PatchField, GeoMesh>::prevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        FatalErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::prevIter() const"
        )   << "previous iteration field" << endl << info() << endl
            << "  not stored."
            << "  Use field.storePrevIter() at start of iteration."
            << abort(FatalError);
    }

    return *fieldPrevIterPtr_;
}


// Calculate and return sum value average
template<class Type, template<class> class PatchField, class GeoMesh>
dimensioned<Type> GeometricField<Type, PatchField, GeoMesh>::average() const
{
    dimensioned<Type> Average
    (
        name() + ".average()",
        dimensions_,
        gAverage(internalField())
    );

    return Average;
}


// Calculate and return weighted average
template<class Type, template<class> class PatchField, class GeoMesh>
dimensioned<Type> GeometricField<Type, PatchField, GeoMesh>::weightedAverage
(
    const scalarField& weightField
) const
{
    return
    (
        dimensioned<Type>
        (
            name() + ".weightedAverage(weights)",
            dimensions_,
            gSum(weightField*internalField())/
            gSum(weightField)
        )
    );
}


// Correct the boundary conditions
template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::correctBoundaryConditions()
{
    // Correct time index
    timeIndex_ = time().timeIndex();

    boundaryField_.evaluate();
}


// Does the field need a reference level for solution
template<class Type, template<class> class PatchField, class GeoMesh>
bool GeometricField<Type, PatchField, GeoMesh>::needReference() const
{
    // Search all boundary conditions, if any are
    // fixed-value or mixed (Robin) do not set reference level for solution.

    static bool searched = false;
    static bool needRef = true;

    if (!searched)
    {
        forAll(boundaryField_, patchi)
        {
            if (boundaryField_[patchi].fixesValue())
            {
                needRef = false;
            }
        }

        reduce(needRef, andOp<bool>());

        searched = true;
    }

    return needRef;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::relax(const scalar alpha)
{
    operator==(prevIter() + alpha*(*this - prevIter()));
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::relax()
{
    scalar alpha = 0;

    if (mesh().relax(name()))
    {
        alpha = mesh().relaxationFactor(name());
    }

    if (alpha > 0)
    {
        relax(alpha);
    }
}


// writeData member function required by regIOobject
template<class Type, template<class> class PatchField, class GeoMesh>
bool GeometricField<Type, PatchField, GeoMesh>::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>&
GeometricField<Type, PatchField, GeoMesh>::null()
{
    GeometricField<Type, PatchField, GeoMesh>* nullPtr = 
        (GeometricField<Type, PatchField, GeoMesh>*)NULL;
    return *nullPtr;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> >
GeometricField<Type, PatchField, GeoMesh>::T() const
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > result
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                name() + ".T()",
                instance(),
                db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensions()
        )
    );

    ::Foam::T(result().internalField(), internalField());
    ::Foam::T(result().boundaryField(), boundaryField());

    return result;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
GeometricField<Type, PatchField, GeoMesh>::component
(
    const direction d
) const
{
    tmp<GeometricField<cmptType, PatchField, GeoMesh> > Component
    (
        new GeometricField<cmptType, PatchField, GeoMesh>
        (
            IOobject
            (
                name() + ".component(" + ::Foam::name(d) + ')',
                instance(),
                db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensions()
        )
    );

    ::Foam::component(Component().internalField(), internalField(), d);
    ::Foam::component(Component().boundaryField(), boundaryField(), d);

    return Component;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
     >& gcf
)
{
    internalField().replace(d, gcf.internalField());
    boundaryField().replace(d, gcf.boundaryField());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::max
(
    const dimensioned<Type>& dt
)
{
    ::Foam::max(internalField(), internalField(), dt.value());
    ::Foam::max(boundaryField(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::min
(
    const dimensioned<Type>& dt
)
{
    ::Foam::min(internalField(), internalField(), dt.value());
    ::Foam::min(boundaryField(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::negate()
{
    internalField().negate();
    boundaryField().negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    if (this == &gf)
    {
        FatalErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::operator="
            "(const GeometricField<Type, PatchField, GeoMesh>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, gf, "=");

    // only equate field contents not ID

    dimensions() = gf.dimensions();
    internalField() = gf.internalField();
    boundaryField() = gf.boundaryField();

    // Correct time index
    timeIndex_ = time().timeIndex();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
{
    if (this == &(tgf()))
    {
        FatalErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::operator="
            "(const tmp<GeometricField<Type, PatchField, GeoMesh> >&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    const GeometricField<Type, PatchField, GeoMesh>& gf = tgf();

    checkField(*this, gf, "=");

    // only equate field contents not ID

    dimensions() = gf.dimensions();

    // This is dodgy stuff, don't try it at home.
    internalField().transfer((Field<Type>&)gf.internalField());

    boundaryField() = gf.boundaryField();

    // Correct time index
    timeIndex_ = time().timeIndex();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const dimensioned<Type>& dt
)
{
    dimensions() = dt.dimensions();
    internalField() = dt.value();
    boundaryField() = dt.value();

    // Correct time index
    timeIndex_ = time().timeIndex();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::operator==
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
{
    const GeometricField<Type, PatchField, GeoMesh>& gf = tgf();

    checkField(*this, gf, "==");

    // only equate field contents not ID

    dimensions() = gf.dimensions();
    internalField() = gf.internalField();
    boundaryField() == gf.boundaryField();

    // Correct time index
    timeIndex_ = time().timeIndex();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GeometricField<Type, PatchField, GeoMesh>::operator==
(
    const dimensioned<Type>& dt
)
{
    dimensions() = dt.dimensions();
    internalField() = dt.value();
    boundaryField() == dt.value();

    // Correct time index
    timeIndex_ = time().timeIndex();
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                         \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void GeometricField<Type, PatchField, GeoMesh>::operator op                   \
(                                                                             \
    const GeometricField<TYPE, PatchField, GeoMesh>& gf                       \
)                                                                             \
{                                                                             \
    checkField(*this, gf, #op);                                               \
                                                                              \
    dimensions() op gf.dimensions();                                          \
    internalField() op gf.internalField();                                    \
    boundaryField() op gf.boundaryField();                                    \
                                                                              \
    /* Correct time index */                                                  \
    timeIndex_ = time().timeIndex();                                         \
}                                                                             \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void GeometricField<Type, PatchField, GeoMesh>::operator op                   \
(                                                                             \
    const tmp<GeometricField<TYPE, PatchField, GeoMesh> >& tgf                \
)                                                                             \
{                                                                             \
    operator op(tgf());                                                       \
    tgf.clear();                                                              \
}                                                                             \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void GeometricField<Type, PatchField, GeoMesh>::operator op                   \
(                                                                             \
    const dimensioned<TYPE>& dt                                               \
)                                                                             \
{                                                                             \
    dimensions() op dt.dimensions();                                          \
    internalField() op dt.value();                                            \
    boundaryField() op dt.value();                                            \
                                                                              \
    /* Correct time index */                                                  \
    timeIndex_ = time().timeIndex();                                         \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream& os,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    os.writeKeyword("dimensions") << gf.dimensions() << token::END_STATEMENT
        << nl << nl;

    gf.internalField().writeEntry("internalField", os);
    os  << nl;
    gf.boundaryField().writeEntry("boundaryField", os);

    // Check state of IOstream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const GeometricField<Type, PatchField, GeoMesh>&)"
    );

    return (os);
}


template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream& os,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
{
    os << tgf();
    tgf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "GeometricBoundaryField.C"
#   include "GeometricFieldFunctions.C"

// ************************************************************************* //
