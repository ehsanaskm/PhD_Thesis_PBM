
dimensionedScalar residualAlpha 
(
     residualAlphaProperties.lookup ("residualAlpha") 
);
/*
dimensionedScalar Y_O2Star 
(
     residualAlphaProperties.lookup ("Y_O2Star") 
);

dimensionedScalar maxEpsilon 
(
     residualAlphaProperties.lookup ("maxEpsilon") 
);
*/
/*
dimensionedScalar dminForKLa 
(
     residualAlphaProperties.lookup ("dminForKLa") 
);
*/
/*
dimensionedScalar rhoL ("rhoL",dimensionSet(1,-3,0,0,0,0,0),998.2);
dimensionedScalar muL ("muL",dimensionSet(1,-1,-1,0,0,0,0),0.001);  
dimensionedScalar diffCoef ("diffCoef",dimensionSet(0,2,-1,0,0,0,0),2.01*1e-9); 
*/


    volScalarField moment0(mesh.lookupObject<volScalarField>("moment.0.populationBalance"));
    volScalarField moment1(mesh.lookupObject<volScalarField>("moment.1.populationBalance"));
    volScalarField moment2(mesh.lookupObject<volScalarField>("moment.2.populationBalance"));
    volScalarField moment3(mesh.lookupObject<volScalarField>("moment.3.populationBalance"));
    volScalarField moment4(mesh.lookupObject<volScalarField>("moment.4.populationBalance"));
    volScalarField moment5(mesh.lookupObject<volScalarField>("moment.5.populationBalance"));
    volScalarField moment6(mesh.lookupObject<volScalarField>("moment.6.populationBalance"));
    volScalarField moment7(mesh.lookupObject<volScalarField>("moment.7.populationBalance"));
    volScalarField moment8(mesh.lookupObject<volScalarField>("moment.8.populationBalance"));

    volScalarField sigmaNDF(mesh.lookupObject<volScalarField>("node0.populationBalance.sigma"));

    volScalarField weight1(mesh.lookupObject<volScalarField>("node0.populationBalance.weight"));
    volScalarField weight2(mesh.lookupObject<volScalarField>("node1.populationBalance.weight"));
    volScalarField weight3(mesh.lookupObject<volScalarField>("node2.populationBalance.weight"));
    volScalarField weight4(mesh.lookupObject<volScalarField>("node3.populationBalance.weight"));

    volScalarField abscissa1(mesh.lookupObject<volScalarField>("node0.populationBalance.abscissa"));
    volScalarField abscissa2(mesh.lookupObject<volScalarField>("node1.populationBalance.abscissa"));
    volScalarField abscissa3(mesh.lookupObject<volScalarField>("node2.populationBalance.abscissa"));
    volScalarField abscissa4(mesh.lookupObject<volScalarField>("node3.populationBalance.abscissa"));

    volScalarField epsilon(mesh.lookupObject<volScalarField>("epsilon.water"));
         


    d32=max(moment3/(moment2+dimensionedScalar("SMALL", dimless/dimLength, SMALL))/1000, dimensionedScalar("Min", dimLength, 0.0001));
 
       
    sigma=sigmaNDF;


    W1=weight1;
    W2=weight2;
    W3=weight3;
    W4=weight4;
    L1=abscissa1;
    L2=abscissa2;
    L3=abscissa3;
    L4=abscissa4;


  forAll(d32, cellI)
  {

	 if(alpha1[cellI] > residualAlpha.value())
	  {
		
				d32[cellI] = d32[cellI];
			
	  }


	  else
	  {
			       d32[cellI] = 0.0001;
	  }

  }

  
   //    KLa=6*(2/1.772*sqrt(diffCoef)*pow((min(epsilon,maxEpsilon)*rhoL/muL),0.25))*max(alpha1,residualAlpha)/d32;

      
/*		
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
*/

    //   Info<<min(max(alpha1,residualAlpha)).value()<<endl;
    
 //     d32 = max(d32, min(d32, 10*dimensionedScalar("one", dimLength, 1.0)))

 //      d32 = max(1*dimensionedScalar("one", dimLength, 1.0), min(d32, 10*dimensionedScalar("one", dimLength, 1.0)));

    


 /*   

  forAll(sigma, cellI)
  {

	 if(alpha1[cellI] > residualAlpha.value())
	  {
		
				sigma[cellI] = sigma[cellI];

	  }
	  else
	  {
			       sigma[cellI] = 0.0;
	  }

  }
*/ 
/*
  forAll(KLa, cellI)
  {

	 if(alpha1[cellI] > residualAlpha.value())
	  {
		
				KLa[cellI] = KLa[cellI];

	  }
	  else
	  {
			       KLa[cellI] = 0.0;
	  }

  } 
*/
   Info << nl << " *********************************************************************** "<< endl;
   Info << nl << "           Updating Sauter diameter d32 by using EQMOM " << endl;
   Info << nl << "           min/avg/max d32: " << min(d32).value()*1000
              << "/" << average(d32).value()*1000<< "/" << max(d32).value()*1000<< "mm." << endl;
   Info << nl << "           Calculated Sigma by using EQMOM " << endl;
   Info << nl << "           min/avg/max sigma: " << min(sigma).value()
              << "/" << average(sigma).value()<< "/" << max(sigma).value() << endl;
   Info << nl << "           min/avg/max moment 8: " << min(moment8).value()
              << "/" << average(moment8).value()<< "/" << max(moment8).value() << endl;
   Info << nl << "           min/avg/max moment 7: " << min(moment7).value()
              << "/" << average(moment7).value()<< "/" << max(moment7).value() << endl;
   Info << nl << "           min/avg/max moment 6: " << min(moment6).value()
              << "/" << average(moment6).value()<< "/" << max(moment6).value() << endl;
   Info << nl << "           min/avg/max moment 5: " << min(moment5).value()
              << "/" << average(moment5).value()<< "/" << max(moment5).value() << endl;
   Info << nl << "           min/avg/max moment 4: " << min(moment4).value()
              << "/" << average(moment4).value()<< "/" << max(moment4).value() << endl;
   Info << nl << "           min/avg/max moment 3: " << min(moment3).value()
              << "/" << average(moment3).value()<< "/" << max(moment3).value() << endl;
   Info << nl << "           min/avg/max moment 2: " << min(moment2).value()
              << "/" << average(moment2).value()<< "/" << max(moment2).value() << endl;
   Info << nl << "           min/avg/max moment 1: " << min(moment1).value()
              << "/" << average(moment1).value()<< "/" << max(moment1).value() << endl;
   Info << nl << "           min/avg/max moment 0: " << min(moment0).value()
              << "/" << average(moment0).value()<< "/" << max(moment0).value() << endl;
   Info << nl << "           min/avg/max alpha: " << min(alpha1).value()
              << "/" << average(alpha1).value()<< "/" << max(alpha1).value() << endl;
  /* Info << nl << "           min/avg/max KLa: " << min(KLa).value()
              << "/" << average(KLa).value()<< "/" << max(KLa).value() << endl;*/
   Info << nl << " *********************************************************************** "<< endl;
 	Info << " " << endl;			
			


