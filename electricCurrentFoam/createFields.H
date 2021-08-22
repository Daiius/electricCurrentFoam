
	Info<< "Reading field sigma (conductivity)" << endl;
    volScalarField sigma
    (
		IOobject
		(
			"sigma",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
    );


    Info<< "Reading field voltage\n" << endl;
    volScalarField voltage
    (
        IOobject
        (
            "voltage",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Calculating field current\n" << endl;
    volVectorField current
    (
        IOobject
        (
            "current",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
		mesh
    );

