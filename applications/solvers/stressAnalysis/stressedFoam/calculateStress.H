    if (runTime.outputTime())
    {
        volTensorField gradU = fvc::grad(U);

        volTensorField sigma = rho*(2.0*mu*symm(gradU) + lambda*I*tr(gradU));

        /*
        volTensorField sigma = rho*fvc::reconstruct
        (
            (2.0*mu + lambda)*fvc::snGrad(U)*mesh.magSf()
          + (
                mesh.Sf() & fvc::interpolate
                (
                    mu*gradU.T()
                  - (mu + lambda)*gradU 
                  + I*
                    (
                        lambda*tr(gradU)

#                       ifdef ThermalStresses
                      - threeKalpha*T
#                       endif
                    ),
                    "interpolate(sigmaRot)"
                )
            )
        );
        */

        if (thermalStress)
        {
            sigma -= rho*fvc::reconstruct
            (
                mesh.Sf() * fvc::interpolate(threeKalpha*T)
            );
        }

        volScalarField sigmaEq
        (
            IOobject
            (
                "sigmaEq",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
        );

        Info<< "Max sigmaEq = " << max(sigmaEq).value()
            << endl;

        volScalarField sigmaxx
        (
            IOobject
            (
                "sigmaxx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::XX)
        );

        volScalarField sigmayy
        (
            IOobject
            (
                "sigmayy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::YY)
        );

        volScalarField sigmazz
        (
            IOobject
            (
                "sigmazz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::ZZ)
        );

        Info<< "Max sigmazz = " << max(sigmazz).value()
            << endl;

        volScalarField sigmaxy
        (
            IOobject
            (
                "sigmaxy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::XY)
        );

        volScalarField sigmaxz
        (
            IOobject
            (
                "sigmaxz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::XZ)
        );

        volScalarField sigmayz
        (
            IOobject
            (
                "sigmayz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::YZ)
        );

        runTime.write();
    }
