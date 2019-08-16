    Info<< "Creating twoPhaseSystem\n" << endl;

    twoPhaseSystem fluid(mesh, g);

    phaseModel& phase1 = fluid.phase1();
    phaseModel& phase2 = fluid.phase2();


    volScalarField& alpha1 = phase1;
    volScalarField& alpha2 = phase2;

    volVectorField& U1 = phase1.U();
    surfaceScalarField& phi1 = phase1.phi();
    surfaceScalarField& alphaPhi1 = phase1.phiAlpha();

    volVectorField& U2 = phase2.U();
    surfaceScalarField& phi2 = phase2.phi();
    surfaceScalarField& alphaPhi2 = phase2.phiAlpha();

    surfaceScalarField& phi = fluid.phi();

    dimensionedScalar pMin
    (
        "pMin",
        dimPressure,
        fluid.lookup("pMin")
    );

    rhoThermo& thermo1 = phase1.thermo();
    rhoThermo& thermo2 = phase2.thermo();

    volScalarField& p = thermo1.p();

    volScalarField& rho1 = thermo1.rho();
    const volScalarField& psi1 = thermo1.psi();

    volScalarField& rho2 = thermo2.rho();
    const volScalarField& psi2 = thermo2.psi();

	 volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.U()
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.rho()
    );

    Info<< "Creating field kinetic energy k\n" << endl;
    volScalarField k("k", 0.5*magSqr(U));

    scalar Cmu=0.09;
    scalar dh = 0.1;

    Info<< "Creating field dissipation of kinetic energy epsilon\n" << endl;
	 volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
				"epsilon",
				sqr(dimLength)/pow3(dimTime),
				Foam::pow(Cmu,3./4.)*Foam::pow(max(k).value(), 3./2.)/0.007/dh
		  )
    );

	 // ********************************************************************************************
	 // Ehsan Askari ,
	 // declaration of population balance packages
	 // ********************************************************************************************

	 volScalarField source_ini
	 (
	     IOobject
	     (
			   "source_ini",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
		  ),
	     mesh
	 );

   volScalarField m_ini
   (
       IOobject
       (
         "m_ini",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
       mesh
   );

	 volScalarField f_ini
	 (
	     IOobject
	     (
           "f_ini",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
		  ),
	     mesh
	 );

   volScalarField w_ini
	 (
	     IOobject
	     (
			   "w_ini",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
		  ),
	     mesh
	 );

	 volScalarField L_ini
	 (
	     IOobject
	     (
			   "L_ini",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
		  ),
	     mesh
	 );


    PtrList<volScalarField> source(3000);
    PtrList<volScalarField> f(3000);
    PtrList<volScalarField> w(3000);
    PtrList<volScalarField> L(3000);
    PtrList<volScalarField> m(3000);




    // read PBM class from 'populationModels' packages
	 populationBalanceModel PBM
	 (
		  source, 					// source term
  		  f,							// VF of each class i
  //      m,
		  alpha1						// total volume fraction of dispersed phase
	 );

  //       populationBalanceReturn PBM2;

	 // ************************************************************************************************

    Info<< "Calculating field DDtU1 and DDtU2\n" << endl;

    volVectorField DDtU1
    (
        "DDtU1",
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1
    );

    volVectorField DDtU2
    (
        "DDtU2",
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2
    );

    volScalarField rAU1
    (
        IOobject
        (
            IOobject::groupName("rAU", phase1.name()),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0)
    );

    volScalarField rAU2
    (
        IOobject
        (
            IOobject::groupName("rAU", phase2.name()),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K1(IOobject::groupName("K", phase1.name()), 0.5*magSqr(U1));
    volScalarField K2(IOobject::groupName("K", phase2.name()), 0.5*magSqr(U2));
