Info<< "Reading the concentrations \n" << endl;

// prepare list to store the concentrations in
PtrList<volScalarField> MFrList(specieNames.size());

forAll(specieNames,countI)
{
    word currentSp = specieNames[countI];
    Info<< "Reading the concentration of the: " << currentSp << endl;

    MFrList.set
    (
	countI,
	new volScalarField
	(
	    IOobject
	    (
		currentSp,
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	)
    );
}

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
	"U",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

// -- dummy entry to resolve in fvOptions source 
volScalarField dummyConcEqn
(
    IOobject
    (
        "dummyConcEqn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("dummyConcEqn",dimMoles/dimVolume,1.0)
);
//dummy volVectorField for dimension correction in equation
//NOTE: I dont know why this is necessary, dimensions are correct
//but openFOAM doesnt think so :(
volScalarField dummyEEqn
(
    IOobject
    (
        "dummyEEqn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("dummyEEqn",dimMass/dimLength*Foam::pow(dimTime,-2)/dimTemperature,1.0)
    //~ dimensionedScalar ("dummy2",dimMass/dimLength*Foam::pow(dimTime,-2)/dimTemperature,1.0)
);

//~ volScalarField p
//~ (
    //~ IOobject
    //~ (
        //~ "p",
        //~ runTime.timeName(),
        //~ mesh,
        //~ IOobject::MUST_READ,
        //~ IOobject::AUTO_WRITE
    //~ ),
    //~ mesh
//~ );

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field isOnWallCoat\n" << endl;
volScalarField isOnWallCoat
(
    IOobject
    (
        "isOnWallCoat",
        runTime.timeName(0),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field isInWallCoat\n" << endl;
volScalarField isInWallCoat
(
    IOobject
    (
        "isInWallCoat",
        runTime.timeName(0),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


//~ #include "createPhi.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);


fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);


#include "compressibleCreatePhi.H"

pressureControl pressureControla(p, rho, simple.dict());

mesh.setFluxRequired(p.name());


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"




