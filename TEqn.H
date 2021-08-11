{
    const volScalarField& vDotc = mixture->vmCond()/rho1;
    const volScalarField& vDotv = mixture->vmEvap()/rho1;
    const volScalarField vDotNetto(vDotc - vDotv);

	TSat = mixture->TSat();
    //alphat = turbulence->nut()/Prt;
    //alphat.correctBoundaryConditions();
	//alphaEff = mixture->alphaEff() + alphat;
	alphaEff = mixture->alphaEff();
	surfaceScalarField alphaEfff = mixture->alphaEfff();
	
	// trzeba dodac kt dla turbulencji
	//surfaceScalarField kf = mixture->kfHarmonic();
	cp = mixture->cp();
	

	fvScalarMatrix TEqn
	(
		  fvm::ddt(T)
		+ fvm::div(phi, T)
	    - fvm::Sp(fvc::div(phi),T)
		- fvm::laplacian(alphaEff, T,"laplacian(alphaEff,T)")
	 ==
		  vDotNetto*hEvap/rhoCp
//		  fvm::Sp(vDotcmvT/rhoCp, T) - vDotcmvT/rhoCp*TSat	
	);

    TEqn.relax();
    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
