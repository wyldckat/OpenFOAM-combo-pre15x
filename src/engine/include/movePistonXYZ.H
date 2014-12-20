{
    scalar deltaZ = runTime.pistonDisplacement();
    Info<< "deltaZ = " << deltaZ << endl;

    pointField newPoints = mesh.points();

    if (layeredMesh)
    {
        forAll (newPoints, pointi)
        {
            point& p = newPoints[pointi];

            if (p.z() < pistonPosition)       // In piston bowl
            {
                p.z() += deltaZ;
            }
            else if (p.z() < deckHeight)      // In liner region
            {
                p.z() += 
                    deltaZ
                   *(deckHeight - p.z())/(deckHeight - pistonPosition);
            }
        }
    }
    else
    {
        tetPolyMesh tetMesh(mesh);

        tetPointVectorField motionU
        (
            IOobject
            (
                "motionU",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh,
            dimensionedVector("0", dimLength, vector::zero),
            fixedValueTetPolyPatchVectorField::typeName
        );

        motionU.boundaryField()[pistonIndex] == vector(0, 0, deltaZ);

        {
            scalarField linerPoints =
                tetMesh.boundary()[linerIndex]
                .localPoints().component(vector::Z);

            motionU.boundaryField()[linerIndex] ==
                vector(0, 0, deltaZ)*pos(deckHeight - linerPoints)
                *(deckHeight - linerPoints)/(deckHeight - pistonPosition);
        }

        scalarField diffusion(mesh.nCells(), 1.0);

        const fvPatchList& patches = mesh.boundary();

        forAll (patches, patchI)
        {
            const labelList::subList fc = patches[patchI].faceCells();

            forAll (fc, fcI)
            {
                diffusion[fc[fcI]] = 2;
            }
        }

        solve(tetFem::laplacian(diffusion, motionU), 1e-6, 0);

        newPoints +=
            (vectorField&)vectorField::subField
            (
                motionU.internalField(),
                newPoints.size()
            );
    }

    mesh.movePoints(newPoints);

    pistonPosition += deltaZ;
    scalar pistonSpeed = deltaZ/runTime.deltaT().value();

    Info<< "clearance: " << deckHeight - pistonPosition << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;

#   include "volContinuity.H"

#   ifdef CheckMesh
    mesh.checkMesh();
#   endif

#   ifndef NoFields
        U.boundaryField()[pistonIndex] == vector(0, 0, pistonSpeed);
#   endif

}
