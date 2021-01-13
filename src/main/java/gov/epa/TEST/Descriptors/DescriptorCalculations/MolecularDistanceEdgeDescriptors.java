package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.lang.reflect.Field;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;



public class MolecularDistanceEdgeDescriptors {

	static int Method=1;
	// if method=1, atom degree= # of bonds attached to the atom,
	// if method=2, degree = valence - # hydrogens

	// method #1 seems to give better predicted values for the 
	// boiling point of alkenes (the BP model was fit to alkanes)
	
	
	//TODO: in ADAPT, degree = valence - # hydrogens (IBRID(I)=IVALENCE-NH)
	// currently I am using degree = number of bonds
	
	
//	public static void Calculate (IAtomContainer m,DescriptorData dd,String[] Fragment,int[][] Distance) {
//	
//		// this method determines degree of atoms based on e-state indices:
//		 
//		CalculatePair(m,dd,Fragment,Distance,"SsCH3","SsCH3","MDEC11");
//		CalculatePair(m,dd,Fragment,Distance,"SsCH3","SssCH2","MDEC12");
//		CalculatePair(m,dd,Fragment,Distance,"SsCH3","SsssCH","MDEC13");
//		CalculatePair(m,dd,Fragment,Distance,"SsCH3","SssssC","MDEC14");
//		
//		CalculatePair(m,dd,Fragment,Distance,"SssCH2","SssCH2","MDEC22");
//		CalculatePair(m,dd,Fragment,Distance,"SssCH2","SsssCH","MDEC23");
//		CalculatePair(m,dd,Fragment,Distance,"SssCH2","SssssC","MDEC24");
//		
//		CalculatePair(m,dd,Fragment,Distance,"SsssCH","SsssCH","MDEC33");
//		CalculatePair(m,dd,Fragment,Distance,"SsssCH","SssssC","MDEC34");
//		
//		CalculatePair(m,dd,Fragment,Distance,"SssssC","SssssC","MDEC44");
//		
//		CalculatePair(m,dd,Fragment,Distance,"SsOH","SsOH","MDEO11");
//		CalculatePair(m,dd,Fragment,Distance,"SsOH","SssO","MDEO12");
//		CalculatePair(m,dd,Fragment,Distance,"SssO","SssO","MDEO22");
//		
//		CalculatePair(m,dd,Fragment,Distance,"SsNH2","SsNH2","MDEN11");  
//		CalculatePair(m,dd,Fragment,Distance,"SsNH2","SssNH","MDEN12");
//		CalculatePair(m,dd,Fragment,Distance,"SsNH2","SsssN","MDEN13");
//
//		CalculatePair(m,dd,Fragment,Distance,"SssNH","SssNH","MDEN22");
//		CalculatePair(m,dd,Fragment,Distance,"SssNH","SsssN","MDEN23");
//		
//		CalculatePair(m,dd,Fragment,Distance,"SsssN","SsssN","MDEN33");
//		CalculateBP(dd);
//		
//
//	}
//	
public static void Calculate2 (IAtomContainer m,DescriptorData dd,int[][] Distance) {
	// this method determines degree of atoms (i.e. primary, secondary, etc.)
	// based on the actual degree of each atom the advantage of this is that it 
	// still works for unsaturated compounds
	
		
		for (int i=1;i<=4;i++) {
			for (int j=i;j<=4;j++) {
				CalculatePair2(m,dd,Distance,"C",i,j);
			}
		}
		
		for (int i=1;i<=2;i++) {
			for (int j=i;j<=2;j++) {
				CalculatePair2(m,dd,Distance,"O",i,j);
			}
		}
		
		for (int i=1;i<=3;i++) {
			for (int j=i;j<=3;j++) {
				CalculatePair2(m,dd,Distance,"N",i,j);
			}
		}
	
	
//		CalculateBP(dd);
		
//		for (int i=1;i<=4;i++) {
//			for (int j=i;j<=4;j++) {
//				try {
//					String strvar="MDEC"+i+j;
//					Field myField = dd.getClass().getField(strvar);
//					System.out.println(strvar+"\t"+myField.getDouble(dd));
//					
//					//System.out.println(strvar+"="+myField.getDouble(dd));
//					} catch (Exception e) {
//						e.printStackTrace();
//					}
//			}
//		}
		
	
	}
	

	
//	private static void CalculatePair(IAtomContainer m, DescriptorData dd, String [] Fragment,int[][] Distance,String frag1,String frag2,String strvar) {
//		int nij=0;
//		double dij=1;
//		double lamdaij=0;
//		
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			for (int j=i+1;j<=m.getAtomCount()-1;j++) {				
//				if ((Fragment[i]==frag1 && Fragment[j]==frag2) || (Fragment[i]==frag2 && Fragment[j]==frag1)) {
//					nij++;
//					dij*=Distance[i][j];
//					//System.out.println((i+1)+"\t"+(j+1));
//				}
//				
//			}
//		}
//		
//		if (nij==0) {
//			lamdaij=0;
//		} else {
//			dij=Math.pow(dij,1.0/(2.0*nij));		
//			lamdaij=nij/(dij*dij);
//		}
//		
//		try {
//		Field myField = dd.getClass().getField(strvar);
//		myField.setDouble(dd, lamdaij);
//		
////		System.out.println(strvar+"="+myField.getDouble(dd));
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//		
//		
//		//System.out.println(nij);
//		//System.out.println(dij);
//		//System.out.println(lamdaij);
//	}
	

private static void CalculateBP(DescriptorData dd) {
	double BP=5.983-0.1039*dd.MDEC11-0.1034*dd.MDEC12;
	BP+=-0.06133*dd.MDEC13-0.007161*dd.MDEC14-0.08385*dd.MDEC22;
	BP+=-0.02965*dd.MDEC23-0.01056*dd.MDEC24-0.04224*dd.MDEC33;
	BP+=-0.08573*dd.MDEC34-0.2304*dd.MDEC44;
	BP=266.7-Math.exp(BP);
	System.out.println("BP="+BP);
	
	
}


private static void CalculatePair2(IAtomContainer m, DescriptorData dd, int[][] Distance,String s,int d1,int d2) {
		int nij=0;
		double dij=1;
		double lamdaij=0;
		
		int D1,D2;
		String S1,S2;
		
		String strvar="MDE"+s+d1+d2;
		
		
		
		int valence=0;
		
		if (s.equals("C")) valence=4;
		else if (s.equals("O")) valence=2;
		else if (s.equals("N")) valence=3;
		
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			IAtom atom1=m.getAtom(i);
			
			if (Method==1) D1=m.getConnectedBondsCount(atom1);
			else D1=valence-atom1.getImplicitHydrogenCount();
			
			S1=atom1.getSymbol();
			
			for (int j=i+1;j<=m.getAtomCount()-1;j++) {				
				IAtom atom2=m.getAtom(j);
				
				if (Method==1) D2=m.getConnectedBondsCount(atom2);
				else D2=valence-atom2.getImplicitHydrogenCount();
								
				S2=atom2.getSymbol();
				
				if (S1.equals(s) && S2.equals(s)) {
					if ((d1==D1 && d2==D2) || (d1==D2 && d2==D1)) {
						nij++;
						dij*=Distance[i][j];
						//System.out.println((i+1)+"\t"+(j+1));
					}
				}
				
			}
		}
		
		if (nij==0) {
			lamdaij=0;
		} else {
			dij=Math.pow(dij,1.0/(2.0*nij));		
			lamdaij=nij/(dij*dij);
		}
		
		try {
		Field myField = dd.getClass().getField(strvar);
		myField.setDouble(dd, lamdaij);
		
		//System.out.println(strvar+"="+myField.getDouble(dd));
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		//System.out.println(nij);
		//System.out.println(dij);
		//System.out.println(lamdaij);
	}
	
}
