{
    scalar deltaZ = runTime.pistonDisplacement().value();
    Info<< "deltaZ = " << deltaZ << endl;

    // Position of the top of the static mesh layers above the piston
    scalar pistonPlusLayers = pistonPosition.value() + pistonLayers.value();

    pointField newPoints = mesh.points();

    if (layeredMesh)
    {
        forAll (newPoints, pointi)
        {
            point& p = newPoints[pointi];

            if (p.z() < pistonPlusLayers)       // In piston bowl
            {
                p.z() += deltaZ;
            }
            else if (p.z() < deckHeight.value())      // In liner region
            {
                p.z() += 
                    deltaZ
                   *(deckHeight.value() - p.z())
                   /(deckHeight.value() - pistonPlusLayers);
            }
        }
    }
    else
    {
        tetPolyMesh tetMesh(mesh);

        // Select the set of boundary condition types.  For symmetry planes
        // and wedge boundaries, the slip condition should be used;
        // otherwise, use the fixedValue
        wordList boundaryTypes
        (
            mesh.boundary().size(),
            fixedValueTetPolyPatchScalarField::typeName
        );

        forAll (mesh.boundary(), patchI)
        {
            if
            (
                mesh.boundary()[patchI].type() == symmetryFvPatch::typeName
             || mesh.boundary()[patchI].type() == wedgeFvPatch::typeName
             || mesh.boundary()[patchI].type() == emptyFvPatch::typeName
            )
            {
                boundaryTypes[patchI] = slipTetPolyPatchScalarField::typeName;
            }
        }

        tetPointScalarField motionUz
        (
            IOobject
            (
                "motionUz",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh,
            dimensionedScalar("0", dimLength, 0),
            boundaryTypes
        );

        motionUz.boundaryField()[pistonIndex] == deltaZ;

        {
            scalarField linerPoints =
                motionUz.boundaryField()[linerIndex].patch()
                .localPoints().component(vector::Z);

            motionUz.boundaryField()[linerIndex] ==
                deltaZ*pos(deckHeight.value() - linerPoints)
                *(deckHeight.value() - linerPoints)
                /(deckHeight.value() - pistonPlusLayers);
        }

        elementScalarField diffusion
        (
            IOobject
            (
                "motionDiffusion",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh,
            dimensionedScalar("d", dimless, 1.0)
        );

        const fvPatchList& patches = mesh.boundary();

        forAll (patches, patchI)
        {
            const labelList::subList fc = patches[patchI].faceCells();

            forAll (fc, fcI)
            {
                diffusion[fc[fcI]] = 2;
            }
        }

        solve(tetFem::laplacian(diffusion, motionUz));

        newPoints.replace
        (
            vector::Z,
            newPoints.component(vector::Z)
          + scalarField::subField
            (
                motionUz.internalField(),
                newPoints.size()
            )
        );
    }

#   ifndef NoFields
    bool absolutePhi = false;
    if (mesh.moving())
    {
        phi += fvc::interpolate(rho)*fvc::meshPhi(rho, U);
        absolutePhi = true;
    }
#   endif

    mesh.movePoints(newPoints);

#   ifndef NoFields
    if (absolutePhi)
    {
        phi -= fvc::interpolate(rho)*fvc::meshPhi(rho, U);
    }
#   endif

    pistonPosition.value() += deltaZ;
    scalar pistonSpeed = deltaZ/runTime.deltaT().value();

    Info<< "clearance: " << deckHeight.value() - pistonPosition.value() << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;

#   include "volContinuity.H"

#   ifdef CheckMesh
    mesh.checkMesh();
#   endif

}
