{
	Pair<tmp<volScalarField> > vDotT = mixture->vDotT();
		const volScalarField& vDotcT = vDotT[0]();
		const volScalarField& vDotvT = vDotT[1]();
	    const volScalarField vDotcmvT(vDotcT - vDotvT);

	TSat = mixture->TSat();
    //alphat = turbulence->nut()/Prt;
    //alphat.correctBoundaryConditions();
	//alphaEff = mixture->alphaEff() + alphat;
	//k = alphaEff*rhoCp;
	
	// trzeba dodac kt dla turbulencji
	//k = mixture->k();
	surfaceScalarField kf = mixture->kfHarmonic();
	cp = mixture->cp();

	fvScalarMatrix TEqn
	(
		  rhoCp*(
			  fvm::ddt(T)
			+ fvm::div(phi, T)
			- fvm::Sp(fvc::div(phi), T)
			  )
		- fvm::laplacian(kf, T,"laplacian(kf,T)")
	 ==
		  fvm::Sp(vDotcmvT, T) - vDotcmvT*TSat	
	);

	//fvScalarMatrix TEqn
	//(
	//	  fvm::ddt(rhoCp,T)
	//	+ fvm::div(rhoPhiCpf, T)
	////	- fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCpf), T) // similarly to U in interPhaseChangeFoam
	//	//- fvm::laplacian(k, T,"laplacian(k,T)")
	//	- fvm::laplacian(kf, T,"laplacian(kf,T)")
	// ==
	//	  fvm::Sp(vDotcmvT, T) - vDotcmvT*TSat	
	//);

    TEqn.relax();
    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
