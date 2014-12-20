{
    fvPatchVectorFieldField& Upatches = U.boundaryField();

    const fvPatchTensorFieldField& gradUpatches = gradU.boundaryField();

    const fvPatchVectorFieldField& Apatches = mesh.Sf().boundaryField();
    const fvPatchScalarFieldField& magApatches = mesh.magSf().boundaryField();

    forAll(Upatches, patchI)
    {
        if
        (
            Upatches[patchI].type()
         == fixedGradientFvPatchVectorField::typeName
        )
        {
            fixedGradientFvPatchVectorField& Upatch =
                refCast<fixedGradientFvPatchVectorField>(Upatches[patchI]);

            vectorField n = Apatches[patchI]/magApatches[patchI];

            Upatch.gradient() =
#           include "traction.H"
        }

        if
        (
            Upatches[patchI].type()
         == directionMixedFvPatchVectorField::typeName
        )
        {

            directionMixedFvPatchVectorField& Upatch =
                refCast<directionMixedFvPatchVectorField>(Upatches[patchI]);

            vectorField n = Apatches[patchI]/magApatches[patchI];

            Upatch.refGrad() =
#           include "traction.H"
       }
    }
}
