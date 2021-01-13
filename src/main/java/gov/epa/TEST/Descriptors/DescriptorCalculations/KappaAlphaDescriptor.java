package gov.epa.TEST.Descriptors.DescriptorCalculations;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

import java.util.LinkedList;


public class KappaAlphaDescriptor {
	
	public static void Calculate(IAtomContainer m,LinkedList [] paths,DescriptorData dd) {
		
		
		double [] Radius=CalculateRadii(m);
		
		
		int A=m.getAtomCount();

		// kappa[0]=atomsCount*Math.log10(atomsCount);

		double alpha = 0;

		for (int i = 0; i < A; i++) {
			alpha += Radius[i] / 0.77 - 1;
		}

		if (A == 1) {
			dd.ka1 = 0;
			dd.ka2 = 0;
			dd.ka3 = 0;
		} else {

			double NumSinglePaths = paths[1].size();
			dd.ka1 = (A + alpha) * Math.pow(A + alpha - 1, 2.0);
			dd.ka1 /= Math.pow(NumSinglePaths + alpha, 2.0);

			if (m.getAtomCount() == 2) {
				dd.ka2 = 0;
				dd.ka3 = 0;
			} else {

				double NumDoublePaths = paths[2].size();

				dd.ka2 = (A + alpha - 1)
						* Math.pow(A + alpha - 2, 2.0);
				dd.ka2 /= Math.pow(NumDoublePaths + alpha, 2.0);
				;

//				System.out.println(paths[3].size());
				
				if (A == 3) {
					dd.ka3 = 0;
				} else {

					double NumTriplePaths = paths[3].size();

					if (NumTriplePaths > 0) {
						if (A % 2 != 0) { // odd
							dd.ka3 = (A + alpha - 1)
									* Math.pow(A + alpha - 3, 2.0);
						} else { // even
							dd.ka3 = (A + alpha - 3)
									* Math.pow(A + alpha - 2, 2.0);
						}
						dd.ka3 /= Math.pow(NumTriplePaths + alpha, 2.0);
						;
					} else
						dd.ka3 = 0;

				}
			}
		}

		
		dd.phia=dd.ka1*dd.ka2/m.getAtomCount();
		
		// System.out.println("kappaalpha1 ="+kappaalpha[1]);
		// System.out.println("kappaalpha2 ="+kappaalpha[2]);
		// System.out.println("kappaalpha3 ="+kappaalpha[3]);
	}

	static private double [] CalculateRadii(IAtomContainer m) {

		// this sub determines covalent radii for determination of ka

		// radii taken from pg 250, todeschini book
		// here D=#attached atoms, DV=#bonds attached

		double[] Radius= new double[m.getAtomCount()]; // covalent radius;
		
		for (int I = 0; I <= m.getAtomCount() - 1; I++) {

			double BondCount=m.getConnectedBondsCount(m.getAtom(I));
			double BondSum=m.getBondOrderSum(m.getAtom(I));
						
			
			/*AtomType at=null;
			try {
				at=this.hsam.findMatchingAtomType(m,m.getAtom(I));
				
			} catch (Exception e) {
				//System.out.println(e);
			}
			
			if (at instanceof AtomType) {
				System.out.println(at.getAtomTypeName());
			} else {
				
			}
			*/
			
			String symbol = m.getAtom(I).getSymbol();
			int charge = m.getAtom(I).getFormalCharge();

					
			if (symbol.equals("C")) // 40 (carbon)
			{
				if (BondSum - BondCount == 0) { // sp3
					Radius[I] = 0.77;
				} else if (BondSum - BondCount == 1) { // sp2
					Radius[I] = 0.67;
				} else { // sp
					Radius[I] = 0.6;
				}
			} else if (symbol.equals("N")) {
				if (BondSum - BondCount == 0) { // sp3
					Radius[I] = 0.74;

				} else if (BondSum - BondCount == 1) { // sp2
					Radius[I] = 0.62;
				} else { // sp
					Radius[I] = 0.55;
				}

			} else if (symbol.equals("O")) {
				if (BondSum - BondCount == 0) { // sp3
					Radius[I] = 0.74;
				} else { // sp2
					Radius[I] = 0.62;
				}
			} else if (symbol.equals("P")) {
				// TODO: find better way to determine hybridization- matches dragon though
				
				// for all currently known estate types, radius should be 1.1 according to dragon
				
				if (BondSum <= 2) { 
					if (BondSum - BondCount == 0) { // sp3
						Radius[I] = 1.1;
					} else { // sp2
						Radius[I] = 1.;
					}
				} else {
					Radius[I] = 1.1; // assume like sp3
				}
				
				//System.out.println("RadiusP="+Radius[I]);
				
			} else if (symbol.equals("S")) {

				// TODO: find better way to determine hybridization
				if (BondSum <= 2) {
					if (BondSum - BondCount == 0) { // sp3
						Radius[I] = 1.04;

					} else { // sp2
						Radius[I] = 0.94;
					}
				} else {
					Radius[I] = 1.04; // assume like sp3
				}
				//System.out.println("RadiusS="+Radius[I]);
				// Radius[I]=0.94;

			} else if (symbol.equals("F")) {
				Radius[I] = 0.72;
				// value in general chemistry book = 0.64
			} else if (symbol.equals("Cl")) {
				Radius[I] = 0.99;
				// value in general chemistry book = 0.99
			} else if (symbol.equals("Br")) {
				Radius[I] = 1.14;
				// value in general chemistry book = 1.14
			} else if (symbol.equals("I")) {
				Radius[I] = 1.33;
				// value in general chemistry book = 1.33
			} else if (symbol.equals("Sn")) {
				Radius[I] = 1.42; // is Sn always sp3?
				// back calc from dragon
				
			} else if (symbol.equals("Si")) {
				Radius[I] = 1.17; // from Dragon help manual
				// QSAR
				// value in general chemistry book = 1.17
			} else if (symbol.equals("Hg")) {
				Radius[I] = 1.48; // determined by back calculating from Dragon
				// QSAR
			} else if (symbol.equals("As")) {
				Radius[I] = 1.2166; // determined by back calculating from MDL
				// QSAR
				// value in general chemistry book = 1.21
			} else if (symbol.equals("Pb")) {
				Radius[I] = 1.535; // from back calc dragon
				// http://www.ccdc.cam.ac.uk/products/csd/radii/ gives 1.54

			} else {
				Radius[I] = 1; // should we output error for unknown?
				//System.out.println(symbol);
			}
		}
		return Radius;
	}

}// end class
