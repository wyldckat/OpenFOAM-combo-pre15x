    dimensionedScalar StCorr("StCorr", dimless, 1.0);

    if (ign.igniting())
    {
        // Calculate volume of ignition kernel
        dimensionedScalar Vk = gSum(c*mesh.V());

        if (Vk.value() > SMALL)
        {
            // Calculate kernel area from it's volume

            // assuming it is part-spherical
            scalar sphereFraction = 1.0;
            dimensionedScalar Ak = 
                sphereFraction*4.0*physicalConstant::pi
               *pow(3.0*Vk/(sphereFraction*4.0*physicalConstant::pi),2.0/3.0);

            // assuming it is part-circular
            //dimensionedScalar thickness("l", dimLength, 0.1);
            //scalar circleFraction = 0.5;
            //dimensionedScalar Ak = 
            //    circleFraction*physicalConstant::pi*thickness
            //   *sqrt(4.0*Vk/(circleFraction*thickness*physicalConstant::pi));

            // Calculate kernel area from b field consistent with the
            // discretisation of the b equation.
            volScalarField mgb =
                fvc::div(nf, b, "div(phiSt,b)") - b*fvc::div(nf) + dMgb;
            dimensionedScalar AkEst = gSum(mgb*mesh.V());

            StCorr.value() = max(min((Ak/AkEst).value(), 10.0), 1.0);

            Info<< "StCorr = " << StCorr.value() << endl;
        }
    }
