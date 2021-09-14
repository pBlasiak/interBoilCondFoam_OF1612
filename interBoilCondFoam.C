/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interBoilCondFoam

Group
    grpMultiphaseSolvers

Description
    Solver for 2 incompressible, immiscible fluids with phase-change
    (boiling, condensation).  Uses a VOF (volume of fluid) phase-fraction based
    interface capturing approach. Solver is based on interFoam.

    The momentum and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    The set of phase-change models provided are designed to simulate boiling
    and condensation.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "getCellDims.H"
	#include "GalusinskiVigneauxNo.H"
    #include "setInitialDeltaT.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "createTimeControls.H"
	    #include "GalusinskiVigneauxNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "createTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "GalusinskiVigneauxNo.H"
            #include "alphaGalusinskiVigneauxNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture->correct();

            #include "UEqn.H"
			#include "TControls.H"
            #include "TEqnSubCycle.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

		// czy to jest tu potrzebne?
		// chyba do aktualizacji Tsat bylo?
		// albo robic korekte na przed petla albo po niej
		// trzeba sprawdziac co robi correct w danej klasie i zastanowic
		// sie czy wszystko ma byc aktualizowane (np. interfejs i nu)
        mixture->correct();

		//if (printWallHeatFluxes)
		//{
		//	#include "calcWallHeatFluxes.H"
		//}

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
