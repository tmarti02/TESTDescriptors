package gov.epa.TEST.Descriptors.DescriptorCalculations;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

//import org.openscience.cdk.Molecule;
import org.openscience.cdk.RingSet;

public class MolecularPropertyDescriptors {

	IAtomContainer m;
	DescriptorData dd;
	String[] Fragment;
	
	public void Calculate(IAtomContainer m, DescriptorData dd, String[] Fragment) {
		this.m = m;
		this.dd = dd;
		this.Fragment=Fragment;
		
		this.CalculateHy();
		this.CalculateUI();
	}
	
	/**
	 * Calculate Hydrophilicity descriptor - see pg 225 of the Handbook of Molecular Descriptors
	 */
	private void CalculateHy() {
		int A=m.getAtomCount();
		int NHy=0;
		int NC=0;
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			
			String s=m.getAtom(i).getSymbol();
			if (s.equals("O") || s.equals("S") || s.equals("N")) {
				if (m.getAtom(i).getImplicitHydrogenCount()>=1) {
					NHy+=m.getAtom(i).getImplicitHydrogenCount();// need to read Todeschini et al., 1997b to find out whether to use 2 or 1 for NH2 group
//					NHy+=m.getAtom(i).getHydrogenCount(); // to match dragon- use instead of NHy++
				}
			} else if (s.equals("C")) {
				NC++;
			}
			
		}
		
		//System.out.println(NHy+"\t"+NC+"\t"+A);
		double Hy=0;
		
		double term1=(1+NHy)*Log(2,(1.0+NHy));
		double term2=NC*(1.0/A*Log(2,1.0/A));
		double term3=Math.sqrt((double)NHy/(A*A));
		double term4=Log(2,1.0+A);
		dd.Hy=(term1+term2+term3)/term4;
		//System.out.println(term1+"\t"+term2+"\t"+term3+"\t"+term4);
		//System.out.println("Hy="+dd.Hy);
		
		
	}
	
//	private void CalculateUI() {
//		int NC=dd.nC;
//		int NH=dd.nH;
//		int NX=dd.nX;
//		int NN=dd.nN;
//		int NP=dd.nP;
//		int CountOS=0; // number of oxygens bonded to sulphur
//		int NSO3=0; //number of SO3 groups
//		int C=dd.nCIC;
//		
//		
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			
//			int CountSdO=0;
//			int CountSsO=0;
//			IAtom ai=m.getAtom(i);
//			IAtom[] ca = m.getConnectedAtoms(m.getAtom(i));
//			String s=ai.getSymbol();
//			
//			if (s.equals("S")) {
//				for (int j=0;j<=ca.length-1;j++) {
//					String s2=ca[j].getSymbol();
//					if (s2.equals("O")) {
//						CountOS++;
//
//						if (m.getBond(ai,ca[j]).getOrder()==2) {
//							CountSdO++;
//						} else if (m.getBond(ai,ca[j]).getOrder()==1) {
//							CountSsO++;
//						}
//					}
//					
//				}
//				if (CountSdO==2 && CountSsO==1) {
//					NSO3++;
//				}
//			}
//		}
//
//		//		try {
////	    String strval=(String)dd.FragmentList.get("SO3");
////	    NSO3=Integer.parseInt(strval);
////		} catch(Exception e) {
////			System.out.println("need to find fragments before calc Ui");
////			dd.Ui=-9999;
////			return;
////		}
//	    
//		
//		double term1=2*NC+2-NH-NX+NN+NP;
//		double term2=2*(CountOS-NSO3);
//		double b=(term1+term2)/2.0-C;
//		System.out.println(term1+"\t"+term2);
//		dd.Ui=Log(2,1.0+b);
////		System.out.println("b="+b);
////		System.out.println("NOS="+NOS);
//		//System.out.println("UI="+UI);
//		
//	}
	
	/**
	 * 
	 */
	
	private void CalculateUI() {
		dd.Ui=Log(2,1+dd.nDB+dd.nTB+dd.nAB);
		
	}
	
	
	private static double Log(int base,double x) {
		
		double Logbx=0;
		
		Logbx=Math.log10(x)/Math.log10((double)base);
		
		return Logbx;
		
		
	}
}
