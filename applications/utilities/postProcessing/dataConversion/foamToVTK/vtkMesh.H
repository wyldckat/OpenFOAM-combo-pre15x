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

Class
    vtkMesh

Description
    Encapsulation of VTK mesh data. Holds mesh or meshsubset and
    polyhedral-cell decomposition on it.

SourceFiles
    vtkMesh.C

\*---------------------------------------------------------------------------*/

#ifndef vtkMesh_H
#define vtkMesh_H

#include "fvMeshSubset.H"
#include "vtkTopo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class Time;

/*---------------------------------------------------------------------------*\
                           Class vtkMesh Declaration
\*---------------------------------------------------------------------------*/

class vtkMesh
:
    public fvMeshSubset
{
    // Private data

        //- Current cellSet (or empty)
        const word setName_;

        //- cellMap, pointMap for if not running subset
        mutable labelList* cMapPtr_;
        mutable labelList* pMapPtr_;

        //- Current decomposition of topology
        mutable vtkTopo* topoPtr_;



    // Private Member Functions

        //- Recalculate the constant mapping arrays
        void calcMapping() const;

        //- Disallow default bitwise copy construct
        vtkMesh(const vtkMesh&);

        //- Disallow default bitwise assignment
        void operator=(const vtkMesh&);


public:

    // Constructors

        //- Construct from components
        vtkMesh(const IOobject& io, const word& setName = "");


    // Destructor

        ~vtkMesh();


    // Member Functions

        // Access

            //- Check if running subMesh
            bool useSubMesh() const
            {
                return setName_.size() > 0;
            }

            //- topology
            const vtkTopo& topo() const
            {
                if (!topoPtr_)
                {
                    topoPtr_ = new vtkTopo(mesh());
                }
                return *topoPtr_;
            }

            //- Access either mesh or submesh
            const fvMesh& mesh() const
            {
                if (setName_.size() == 0)
                {
                    return *this;
                }
                else
                {
                    return subMesh();
                }
            }

            //- Cell mapping array.
            const labelList& cMap() const
            {
                if (useSubMesh())
                {
                    return cellMap();
                }
                else
                {
                    if (!cMapPtr_)
                    {
                        calcMapping();
                    }
                    return *cMapPtr_;
                }
            }

            //- Point mapping array.
            const labelList& pMap() const
            {
                if (useSubMesh())
                {
                    return pointMap();
                }
                else
                {
                    if (!pMapPtr_)
                    {
                        calcMapping();
                    }
                    return *pMapPtr_;
                }
            }



        // Edit

            //- Read mesh
            polyMesh::readUpdateState readUpdate();


            //- Map volume field (does in fact do very little interpolation;
            //  just copied from fvMeshSubset)
            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>& fld
            ) const
            {
                if (setName_.size() == 0)
                {
                    return fld;
                }
                else
                {
                    tmp<GeometricField<Type, fvPatchField, volMesh> > subFld =
                        fvMeshSubset::interpolate(fld);

                    subFld().rename(fld.name());

                    return subFld;
                }
            }

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
