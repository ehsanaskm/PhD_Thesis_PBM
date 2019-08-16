
dimensionedScalar residualAlpha 
(
     residualAlphaProperties.lookup ("residualAlpha") 
);

dimensionedScalar Y_O2Star 
(
     residualAlphaProperties.lookup ("Y_O2Star") 
);

dimensionedScalar maxEpsilon 
(
     residualAlphaProperties.lookup ("maxEpsilon") 
);

dimensionedScalar rhoL ("rhoL",dimensionSet(1,-3,0,0,0,0,0),998.2);
dimensionedScalar muL ("muL",dimensionSet(1,-1,-1,0,0,0,0),0.001);  
dimensionedScalar diffCoef ("diffCoef",dimensionSet(0,2,-1,0,0,0,0),2.01*10e-9); 
dimensionedScalar dair ("dair",dimensionSet(0,1,0,0,0,0,0),0.0035);


   volScalarField epsilon(mesh.lookupObject<volScalarField>("epsilon.water"));
         

  KLa=6*(2/M_PI*sqrt(diffCoef)*pow((min(epsilon,maxEpsilon)*rhoL/muL),0.25))*min(alpha1,residualAlpha)/dair;



		
				volScalarField SY_O2
				(
				    IOobject
				    (
					"SY_O2",
					runTime.timeName(),
					mesh,
			/////////////////////////////////////
				 IOobject::NO_READ,
				 IOobject::AUTO_WRITE
			/////////////////////////////////////
				    ),
				    KLa*(Y_O2Star-Y_O2)
				);


				fvScalarMatrix Y_O2Eqn
				(
				      fvm::ddt(Y_O2)
				    + fvm::div(phi2, Y_O2)
				    - fvm::laplacian(diffCoef, Y_O2)
				==
				      fvc::Su(SY_O2,Y_O2)

				);

			   Y_O2Eqn.relax();

			   Y_O2Eqn.solve();

			   Y_O2=max(min(Y_O2,Y_O2Star),1e-15);  
			




 //  Y_O2=min(Y_O2,Y_O2Star); 
 //  Y_O2=max(Y_O2,1e-15);

   Info << nl << " *********************************************************************** "<< endl;
   Info << nl << "           min/avg/max alpha: " << min(alpha1).value()
              << "/" << average(alpha1).value()<< "/" << max(alpha1).value() << endl;
   Info << nl << "           min/avg/max Y_O2: " << min(Y_O2).value()
              << "/" << average(Y_O2).value()<< "/" << max(Y_O2).value() << endl;
   Info << nl << "           min/avg/max KLa: " << min(KLa).value()
              << "/" << average(KLa).value()<< "/" << max(KLa).value() << endl;
   Info << nl << " *********************************************************************** "<< endl;
 	Info << " " << endl;			
			


