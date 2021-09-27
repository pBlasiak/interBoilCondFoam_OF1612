{
	Pair<tmp<volScalarField> > vDotT = mixture->vDotT();
		const volScalarField& vDotcT = vDotT[0](); // it should be 0 for evaporation and (-) for condensation
		const volScalarField& vDotvT = vDotT[1](); // it should be (+) for evaporation and 0 for condensation
	    const volScalarField vDotcmvT(vDotcT - vDotvT);

	Info<< "vDotcT = " << vDotcT << endl; 
	Info<< "vDotvT = " << vDotvT << endl; 

	TSat = mixture->TSat();
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();
	alphaEff = mixture->k()/rhoCp + alphat;
	
	// trzeba zmienic w BC dla q na rhoCp i wtedy cp nie bedzie potrzebne
	//cp = mixture->cp();
	

	fvScalarMatrix TEqn
	(
		  fvm::ddt(T)
		+ fvm::div(phi, T)
	    - fvm::Sp(fvc::div(phi),T)
		- fvm::laplacian(alphaEff, T,"laplacian(alphaEff,T)")
	 ==
		  fvm::Sp(vDotcmvT/rhoCp, T) - vDotcmvT/rhoCp*TSat	
	);

    TEqn.relax();
    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
