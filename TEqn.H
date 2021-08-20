{
	Pair<tmp<volScalarField> > vDotT = mixture->vDotT();
		const volScalarField& vDotcT = vDotT[0]();
		const volScalarField& vDotvT = vDotT[1]();
	    const volScalarField vDotcmvT(vDotcT - vDotvT);

		//Info<< "vDotcT = " << vDotcT << endl;
		//Info<< "vDotvT = " << vDotvT << endl;
		//Info<< "vDotcmvT = " << vDotcmvT << endl;

	TSat = mixture->TSat();
    //alphat = turbulence->nut()/Prt;
    //alphat.correctBoundaryConditions();
	//alphaEff = mixture->alphaEff() + alphat;
	alphaEff = mixture->alphaEff();
	surfaceScalarField alphaEfff = mixture->alphaEfff();
	
	// trzeba dodac kt dla turbulencji
	//surfaceScalarField kf = mixture->kfHarmonic();
	
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
