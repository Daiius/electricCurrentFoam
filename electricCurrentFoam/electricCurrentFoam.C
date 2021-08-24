/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    electricCurrentFoam

Group
	grpElectroMagneticsSolvers

Description
	Solver for electric current with distributed conductivity.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	argList::addNote
	(
		"Solver for electric current with distributed conductivity."
	);

	#include "addCheckCaseOptions.H"
	#include "setRootCaseLists.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;


	while (runTime.loop())
	{
		Info<< "Iteration = " << runTime.timeName() << nl << endl;
		
		solve( fvm::laplacian(sigma, voltage) );
		
		current = fvc::reconstruct(
			- fvc::interpolate(sigma) * fvc::snGrad(voltage) * mesh.magSf()
		);

		runTime.write();
		runTime.printExecutionTime(Info);
	}

    Info<< "End\n" << endl;

	// Boundary condition tests
	auto bf = current.boundaryField();
	auto wlist = bf.types();
	Info<< wlist << endl;
	forAll(wlist, wlistI)
	{
		Info<< wlist[wlistI] << endl;
		if (wlist[wlistI] == "fixedValue")
		{
			Info<< "Current condition found!" << endl;
			const fvPatch &patch = current.boundaryField()[wlistI].patch();
			auto boundaryCurrent = patch.lookupPatchField<volVectorField,vector>("current");
			auto internalCurrent = current.internalField()[wlistI];
			forAll(patch.faceCells(), pcid)
			{
				label cid = patch.faceCells()[pcid];
				const cell &c = mesh.cells()[cid];

				Info<< c << "," << current[cid] << endl;
				
			}

			//forAll(boundaryCurrent, bci)
			//{
			//	Info<< bci << "," << boundaryCurrent[bci] << endl;
			//	Info<< internalCurrent[bci] << endl;
			//}
		}
	}
	//

    return 0;
}


// ************************************************************************* //
