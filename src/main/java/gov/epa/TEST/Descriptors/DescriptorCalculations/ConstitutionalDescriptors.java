package gov.epa.TEST.Descriptors.DescriptorCalculations;

//import org.openscience.cdk.*;


import java.lang.reflect.Field;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.AtomicProperties;

import org.openscience.cdk.Ring;

//import ToxPredictor.Utilities.CDKUtilities;

public class ConstitutionalDescriptors {

	IAtomContainer m; // hydrogen supressed molecule
	DescriptorData dd;
	
	double [] EState;
	AtomicProperties ap;
	IRingSet rs;
	IRingSet srs;
	
	public void Calculate(IAtomContainer m,DescriptorData dd,double [] EState,IRingSet rs,IRingSet srs) {
		
		
//		try {
//			this.m=(IAtomContainer)m.clone();
//			
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		
		this.m=m;
		
		// get rid of charged versions
//		CDKUtilities.FixNitroGroups2(m); 

		// for some reason the following doesnt work:
		//		CDKUtilities.Normalize(m,"ToxPredictor/system/fixsmiles3.xml");
				
		this.dd=dd;
		this.EState=EState;
		this.rs=rs;
		this.srs=srs;
		
		dd.nSK=m.getAtomCount();		
		dd.nBO=m.getBondCount();
		this.DetermineAtomTypeCounts();
		
		
		this.CalculateMultipleBondConstitutionalDescriptors();
		this.Calculate_mw();
		this.Calculate_TotalProps();
//		this.CalculateAromaticRatio();

		dd.nCIR=rs.getAtomContainerCount();
		dd.nCIC=dd.nBO-dd.nSK+1;
		
		
		//TODO: should we use all cycles and rings (rs) or just simple rings (srs)
		this.DetermineRingCounts(rs);
		
//		CDKUtilities.Normalize(m,"ToxPredictor/system/fixsmiles.xml");
		
		
	}
	
	
	
	private void DetermineRingCounts(IRingSet ringset) {
		Field myField;
		try {
		
			int count;
			String strCount;
			
			dd.nBnz=0;
			
			for (int i=3;i<=12;i++) {
				if (i<10) strCount="0"+i;
				else strCount=i+"";
				
				myField = dd.getClass().getField("nR"+strCount);
				myField.setInt(dd,0);
				
			}
			
			for (int i=0;i<=ringset.getAtomContainerCount()-1;i++) {
				
				Ring r=(Ring)ringset.getAtomContainer(i);
				count=r.getAtomCount();
				
		
				if (count<10) strCount="0"+count;
				else strCount=count+"";
				
//				if (count==9) {
//					System.out.println("\nRing#"+i);
//					for (int j=0;j<r.getAtomCount();j++) {
//						System.out.println(r.getAtomNumber(r.getAtom(j)));
//					}
//				}
				
			
//				System.out.println("ringsize="+count);
//				System.out.println("aromatic="+r.getFlag(CDKConstants.ISAROMATIC));
				
				if (count==6 && r.getFlag(CDKConstants.ISAROMATIC)) {
					boolean AllC=true;
					
					for (int j=0;j<=r.getAtomCount()-1;j++) {
						if (!(r.getAtom(j).getSymbol().equals("C"))) {
							AllC=false;
							break;
						}
					}
					
					if (AllC) {
						dd.nBnz++;
					}
				}
				
				if (count>12) continue;
				
				myField = dd.getClass().getField("nR"+strCount);
				int currentcount=myField.getInt(dd);
				currentcount++;
				
				myField = dd.getClass().getField("nR"+strCount);
				myField.setInt(dd,currentcount);
				
			}
		
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	private void DetermineAtomTypeCounts() {
		
		dd.nH=0;
		dd.nC=0;
		dd.nN=0;
		dd.nO=0;
		dd.nP=0;
		dd.nS=0;
		dd.nF=0;
		dd.nCL=0;
		dd.nBR=0;
		dd.nI=0;
		dd.nB=0;
		//dd.nHM=0;
		dd.nX=0;
		
//		System.out.println("enter determine atom type");
		
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			String s=m.getAtom(i).getSymbol();
			
			dd.nH+=m.getAtom(i).getImplicitHydrogenCount(); // set in EStateFragmentDescriptor
			//System.out.println(m.getAtom(i).getHydrogenCount());
			if (s.equals("C")) dd.nC++;
			else if (s.equals("N")) dd.nN++;
			else if (s.equals("O")) dd.nO++;
			else if (s.equals("P")) dd.nP++;
			else if (s.equals("S")) dd.nS++;
			
			else if (s.equals("F")) {
				dd.nF++;
				dd.nX++;
			}
			else if (s.equals("Cl"))  {
				dd.nCL++;
				dd.nX++;
			}
			else if (s.equals("Br")) {
				dd.nBR++;
				dd.nX++;
			}
			else if (s.equals("I")) {
				dd.nI++;
				dd.nX++;
			}
			
			else if (s.equals("B")) dd.nB++;
						
			
			
			
		}
			
		dd.nAT=m.getAtomCount()+dd.nH;
		dd.nBT=m.getBondCount()+dd.nH;
		
	}
	
	
	
	
//	private void CalculateAromaticRatio() {
//		dd.ARR=0;
//			
//		for (int i=0;i<m.getBondCount();i++) {
//			if (m.getBond(i).getFlag(CDKConstants.ISAROMATIC)) {
//				dd.ARR++;
//			}			
//		}
//
////		dd.ARR/=bonds.length;
//		dd.ARR/=m.getBondCount();
//
//		
//	}
	
	
	private void CalculateMultipleBondConstitutionalDescriptors() {
		//number of multiple bonds

		
		dd.nBM=0;
		dd.SCBO=0;
		dd.nDB=0;
		dd.nTB=0;
		dd.nAB=0;
		
		
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			System.out.println(m.getAtom(i).getFlag(CDKConstants.ISAROMATIC));
//		}
		
		
		for (int i=0;i<m.getBondCount();i++) {
			
			double bondorder=0;
			
			//System.out.println(bonds[i].getOrder());							
			
			IAtom atom0=m.getBond(i).getAtom(0);
			IAtom atom1=m.getBond(i).getAtom(1);
			
			if (atom0.getFlag(CDKConstants.ISAROMATIC)
					&& atom1.getFlag(CDKConstants.ISAROMATIC)) {
				
				
				
				// check if in same ring:
				boolean SameRing = EStateFragmentDescriptor.InSameAromaticRing(m,atom0,atom1, rs);
				
				
				if (SameRing) {
					bondorder=1.5;
					dd.nAB++;
					
				} else {					
					bondorder=m.getBond(i).getOrder().numeric();
				}
			} else {
				bondorder=m.getBond(i).getOrder().numeric();
			}
			
			
			
			if (bondorder>1) {
				dd.nBM++;
			}
			dd.SCBO+=bondorder;
			
			if (bondorder==2) {
				dd.nDB++;
			}
			
			if (bondorder==3) {
				dd.nTB++;
			}
			
		} // end loop over bonds
		
		
		dd.ARR=(double)dd.nAB/(double)m.getBondCount();
		
	}
	
	
	private void Calculate_mw() {
		// tried to use CDK built in methods but they suck
		// alternative method would be to use m2 which includes the hydrogens
		
		try {
			ap=AtomicProperties.getInstance();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		dd.MW=0;
				
		for (int i=0;i<=m.getAtomCount()-1;i++) {			
			IAtom a=m.getAtom(i);
			dd.MW+=ap.GetMass(a.getSymbol());
			dd.MW+=a.getImplicitHydrogenCount()*ap.GetMass("H");
			
		}
		dd.AMW=dd.MW/dd.nAT;
			
	}
	private void Calculate_TotalProps() {
		dd.Sv=0;
		dd.Se=0;
		dd.Sp=0;
		dd.Ss=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {			
			IAtom a=m.getAtom(i);
						
			dd.Sv+=ap.GetNormalizedVdWVolume(a.getSymbol());
			dd.Se+=ap.GetNormalizedElectronegativity(a.getSymbol());
			dd.Sp+=ap.GetNormalizedPolarizability(a.getSymbol());
			
			dd.Sv+=ap.GetNormalizedVdWVolume("H")*a.getImplicitHydrogenCount();
			dd.Se+=ap.GetNormalizedElectronegativity("H")*a.getImplicitHydrogenCount();
			dd.Sp+=ap.GetNormalizedPolarizability("H")*a.getImplicitHydrogenCount();
						
			//dd.Ss+=ap.GetNormalized(a.getSymbol());
		}
		
		dd.Mv=dd.Sv/dd.nAT;
		dd.Me=dd.Se/dd.nAT;
		dd.Mp=dd.Sp/dd.nAT;
		
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			dd.Ss+=EState[i];
		}
		dd.Ms=dd.Ss/m.getAtomCount();
		
		
		
		
	}

}


