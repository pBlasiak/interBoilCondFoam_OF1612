{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            mixture->surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    Pair<tmp<volScalarField> > vDotP = mixture->vDotP();
    const volScalarField& vDotcP = vDotP[0](); // it should be 0 for evaporation and (-) for condensation
    const volScalarField& vDotvP = vDotP[1](); // it should be (+) for evaporation and 0 for condensation
    const volScalarField vDotcmvP = vDotcP - vDotvP;
		Info<< "vDotcP = " << vDotcP << endl;
		Info<< "vDotvP = " << vDotvP << endl;
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
//            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
//		-	fvm::SuSp(vDotcP-vDotvP,p_rgh) + (vDotcP-vDotvP)*(pSat-rho*gh) 
//	//	-	fvm::Sp(vDotcP,p_rgh) - vDotvP + vDotcP*(pSat-rho*gh) 
//		//  -	(vDotcP + vDotvP)
            fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
		==
			//fvm::Sp(vDotcP-vDotvP,p_rgh) - (vDotcP-vDotvP)*(pSat-rho*gh) 
		//-	fvm::SuSp(-vDotcmvP,p_rgh) - (vDotcmvP)*(pSat-rho*gh) 
		//-	fvm::Sp(-vDotcmvP,p_rgh) - (vDotcmvP)*(pSat-rho*gh) // chyba daje te same wyniki co z SuSp i troche szybciej liczy ale trzeba sprawdzic fvm::Sp(vDotcmvP,p_rgh) - (vDotcmvP)*(pSat-rho*gh) // chyba daje te same wyniki co z SuSp i troche szybciej liczy ale trzeba sprawdzic
		    fvm::Sp(vDotcmvP,p_rgh) - (vDotcmvP)*(pSat-rho*gh) // chyba daje te same wyniki co z SuSp i troche szybciej liczy ale trzeba sprawdzic
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
    Info<< "min/max(U) = " << min(U).value() << ", "
        << max(U).value() <<endl;
}
