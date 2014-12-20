    Info << "Calculating patch traction" << endl;
    PtrList<vectorField> traction(mesh.boundary().size());

    {
        const fvPatchVectorFieldField& Apatches = mesh.Sf().boundaryField();

        const fvPatchScalarFieldField& magApatches =
            mesh.magSf().boundaryField();

        forAll(traction, patchI)
        {
            vectorField n = Apatches[patchI]/magApatches[patchI];

            traction.hook
                (
                    new vectorField
                    (
                        patchTraction[patchI]
                      - patchPressure[patchI]*n
                    )
                );
        }
    }

