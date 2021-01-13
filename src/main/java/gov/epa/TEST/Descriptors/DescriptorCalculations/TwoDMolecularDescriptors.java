package gov.epa.TEST.Descriptors.DescriptorCalculations;
import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

/**
 * This class calculates 2D molecular descriptors including the number of
 * hydrogen bond acceptors, hydrogen bond donors, and the number of weak
 * hydrogen bond donors.  In addition the class calculates the sum of the E-state
 * values for the hydrogen bond donors and acceptors, and three molecular 
 * and group polarity indices (Qs, Qv, and Qsv).
 * 
 * The definitions for what constitutes hydrogen bond donors and acceptors were taken from Dragon.
 * 
 * @author Todd Martin
 * @modified 8-17-06
 */

public class TwoDMolecularDescriptors {

	// for convenience local copies of the variables passed in the calculate routine are created:
	
	/**  local instance of the molecule that descriptor is being calculated for */
	private IAtomContainer m; 
	
	/**  local instance of the data class to store descriptor values */
	private DescriptorData dd;
	
	/**  local instance of the intrinsic state for each atom in the molecule */
	private double [] IS;

	/**  local instance of the E-state for each atom in the molecule */
	private double [] EState;
	
	/**  local instance of the Hydrogen E-state for each atom in the molecule */
	private double [] HEState;
	
	/** local instance of the E-state fragment name for each atom in the molecule  */
	private String [] Fragment;
		
	/**
	 *  Performs all of the calculations of the class
	 *
	 * @param  m IMolecule - molecule that descriptor is being calculated for
	 * @param  dd DescriptorData - data class to store descriptor values
	 * @param  IS double [] - intrinsic state for each atom in the molecule 
	 * @param EState double [] - E-state for each atom in the molecule
	 * @param HEState double [] - Hydrogen E-state for each atom in the molecule
	 * @param Fragment [] String [] - E-state fragment name for each atom in the molecule
	 */
	
	public void calculate(IAtomContainer m,DescriptorData dd,double [] IS,double [] EState,double [] HEState,String [] Fragment) {
		//set local variables = passed variables
		this.m=m;
		this.dd=dd;
		this.EState=EState;
		this.HEState=HEState;
		this.IS=IS;
		this.Fragment=Fragment;
		
		
//		this.Calculate_nvx();
		this.Calculate_numHBa_SHHBa();
		this.Calculate_Qs_Qv_Qsv();
		
		// now add in hydrogens:
		//CDKUtilities.AddHydrogens(m);
		//CDKUtilities.FixSulphurH(m);
				
		this.Calculate_numHBd_SHHBd();
		this.Calculate_numwHBd_SHwHBd();
		
		
	}
	
	/**
	 * Calculates group polarity indices (Qs, Qv, and Qsv). The equations for
	 * these indices were taken from <a
	 * href="http://www.mdl.com/products/predictive/qsar/index.jsp">MDL QSAR</a>.
	 * 
	 */
	
	private void Calculate_Qs_Qv_Qsv() {
		
		// TODO: add reference for equations in this method
		
		int N=2; // for Carbon
		
		double sumIalk=0;
		double sumI=0;
		double sumImax=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			
			double DV=m.getConnectedBondsCount(i);
			double D=DV;
															
			sumIalk+=(DV + 1) / D;
			sumI+=IS[i];
			
			if (DV==1) {
				sumImax+=8.0; //(a)  All terminal atoms (atoms bonded to one other non-hydrogen atom) replaced by –F (I = 8.000)
			} else if (DV==2) {
				sumImax+=3.5;//(b) All divalent atoms (atoms bonded to two other non-hydrogen atoms) replaced by –O- (I = 3.500)
			} else if (DV==3) {
				sumImax+=2.0; //(c) All trivalent atoms (atoms bonded to three other non-hydrogen atoms) replaced by >N- (I = 2.000)
			} else if (DV==4) {
				sumImax+=1.25;//(d) All quaternary atoms (atoms bonded to four other non-hydrogen atoms) replaced by >C< (I = 1.250)
			} else {
				sumImax+=IS[i];// is this right???
			}
			
			
		}
				
		dd.Qs=Math.pow(m.getAtomCount()/sumI,2.0)*sumIalk;
		dd.Qv=sumIalk*sumImax/(sumI*sumI);
		
		double sumIave=(sumIalk+sumImax)/2.0;
		dd.Qsv=sumIave*sumIalk/(sumI*sumI);
				
	}
	
//	/**
//	 * Determines the number of non-hydrogen atoms
//	 * 
//	 */
//	private void Calculate_nvx() {					
//		// this method assumes hydrogens have been removed first
//		dd.nvx=m.getAtomCount(); 
//	}
	
	
	/*private void Calculate_fw() {
		
		// did not use CDK's WeightDescriptor since it uses weight of most abundant isotope instead of atomic weights
		
		
		// now add in hydrogens:
		HydrogenAdder hydrogenAdder = null;
		hydrogenAdder = new HydrogenAdder("org.openscience.cdk.tools.ValencyChecker");
		//hydrogenAdder = new HydrogenAdder("org.openscience.cdk.tools.SaturationChecker");
		
				
		dd.fw = 0;
		
		try {
			hydrogenAdder.addExplicitHydrogensToSatisfyValency(m);
			//hydrogenAdder.addImplicitHydrogensToSatisfyValency(m2);
			//hydrogenAdder.addHydrogensToSatisfyValency(m2);
									
			
			//System.out.println(m2.getAtomCount());
			//System.out.println(m.getAtomCount());
			
			for (int i = 0; i <= m.getAtomCount()-1; i++) {
				
				//System.out.println("WEIGHT: "+container.getAtom(i).getSymbol() +" " +IsotopeFactory.getInstance().getMajorIsotope( container.getAtom(i).getSymbol() ).getExactMass());
				Isotope[] isotopes=IsotopeFactory.getInstance().getIsotopes(m.getAtom(i).getSymbol());
				
				double atomweight=0;
				Atom a=m.getAtom(i);
				
				double totala=0; // total abundance for less abundant isotopes
				for (int j=0;j<=isotopes.length-1;j++) {
					Isotope isotope=isotopes[j];
					double ab=isotope.getNaturalAbundance();
					if (ab!=100) totala+=ab; 
				}
				
				for (int j=0;j<=isotopes.length-1;j++) {
					Isotope isotope=isotopes[j];
					double ab=isotope.getNaturalAbundance();
					
					// for some lame reason abundance for most abundant
					// isotope is listed as 100 in CDK (instead of 100-sum of
					// others)
					if (ab==100) ab-=totala; 					
					
					atomweight+=isotope.getExactMass()*ab/100.0;
										
				}
				
				//System.out.println((i+1)+"\t"+a.getSymbol()+"\t"+atomweight);
				dd.fw += atomweight;								
				
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
			

	}
	*/
	
	
	
	/*
	public static double GetWeight(String s) {
		// CDKs methods for getting atomic weight are sucky
		double weight=0;
		if (s.equals("As")) {
			weight=74.92159;
		}else if (s.equals("Br")) {
			weight=79.904;
		}else if (s.equals("C")) {
			weight=12.011;
		}else if (s.equals("Cl")) {
			weight=35.4527;			
		}else if (s.equals("F")) {
			weight=18.9984032;
		}else if (s.equals("H")) {
			weight=1.00794;
		}else if (s.equals("Hg")) {
			weight=200.59;					
		}else if (s.equals("I")) {
			weight=126.90447;					
		}else if (s.equals("N")) {
			weight=14.00674;							
		} else if (s.equals("O")) {
			weight=15.9994;
		}else if (s.equals("P")) {
			weight=30.973762;							
		}else if (s.equals("Pb")) {
			weight=207.2;			
		}else if (s.equals("S")) {
			weight=32.066;							
		}else if (s.equals("Si")) {
			weight=28.0855;
		}else if (s.equals("Sn")) {
			weight=118.71;					
				
		} else {
			weight=-9999;
			System.out.println("**Warning: Need weight for"+s);
		}
		
		return weight;
	}
	*/
	
	
	/**
	 * Calculates the number and sum of the E-state values for the hydrogen bond
	 * acceptors in the molecule. In addition the atoms which are hydrogen bond
	 * acceptors are flagged.  The definition for what constitutes a hydrogen bond
	 * acceptor was taken from Dragon version 5.4.
	 * 
	 */
	private void Calculate_numHBa_SHHBa() {
		dd.numHBa=0;
		dd.SHHBa=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			
			IAtom a=m.getAtom(i);			
			String s=a.getSymbol();
			
			//if (s.equals("Cl") || s.equals("F") || s.equals("O")) {
			if (s.equals("F") || s.equals("O")) { // Cl not included in dragon
				dd.numHBa++;
				dd.SHHBa+=EState[i];
				a.setFlag(CDKConstants.IS_HYDROGENBOND_ACCEPTOR,true);
				
								
			} else if (s.equals("N")) {
				
				//System.out.println(Fragment[i]);
				
				if (a.getFormalCharge()==0) {
					if (a.getFlag(CDKConstants.ISAROMATIC)) {
						if (Fragment[i] instanceof String) {
							if (Fragment[i].equals("SaaN")) { // dragon doesnt
																// include SaaNH
																// (pyrrolyl
																// form)
								// if (Fragment[i].equals("SaaN") ||
								// Fragment[i].equals("SaaNH") ) {
								// if (a.getHydrogenCount()==1) {
								dd.numHBa++; // need zero charge to have
												// unpaired electrons
								dd.SHHBa += EState[i];
								a.setFlag(CDKConstants.IS_HYDROGENBOND_ACCEPTOR,true);
							}
						}
						
						
					} else {
						dd.numHBa++; // need zero charge to have unpaired electrons
						dd.SHHBa+=EState[i];
						a.setFlag(CDKConstants.IS_HYDROGENBOND_ACCEPTOR,true);
					}
				}
								
			} 



			// dragon 5 doesnt include S as donor or acceptor						
			/*
			else if (s.equals("S")) {
			// according MDL QSAR, SH atoms aren't acceptors but others are (as
			// long as have unpaired electrons)
				if (m.getBondOrderSum(a)>1 && m.getBondOrderSum(a)<4) {
					
					dd.numHBa++;
					dd.SHHBa+=EState[i];
				}
				
			}
			*/
			
		}
		
	}
	
	/**
	 * Calculates the number and sum of the E-state values for the hydrogen bond
	 * donors in the molecule. In addition the atoms which are hydrogen bond donors
	 * are flagged. The definition for what constitutes a hydrogen bond
	 * donor was taken from Dragon version 5.4.
	 * 
	 */
	private void Calculate_numHBd_SHHBd() {
		dd.numHBd=0;
		dd.SHHBd=0;
							
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			IAtom a=m.getAtom(i);			
			String s=a.getSymbol();
			
			// dragon 5 doesnt include S as donor or acceptor:
//			if (s.equals("O")) {
//				if (a.getHydrogenCount()>0) {
//					dd.numHBd++;
//					dd.SHHBd+=HEState[i];
//					a.setFlag(CDKConstants.IS_HYDROGENBOND_DONOR,true);
//				}
//								
//			} else if  (s.equals("N")) { // only non aromatic NH's are donors
//				if (!a.getFlag(CDKConstants.ISAROMATIC)) {
//					if (a.getHydrogenCount()>0) {
//						dd.numHBd+=a.getHydrogenCount();
//						dd.SHHBd+=HEState[i];
//						a.setFlag(CDKConstants.IS_HYDROGENBOND_DONOR,true);
//					}
//				}
//			}
			
// in the new version of dragon all NH's are donors:
			if (s.equals("O") || s.equals("N") ) {
				dd.numHBd+=a.getImplicitHydrogenCount();	
				dd.SHHBd+=HEState[i];
				a.setFlag(CDKConstants.IS_HYDROGENBOND_DONOR,true);
			}								

			
		} // end i for loop
		
	}
	
	private void Calculate_numwHBd_SHwHBd() {
		dd.numwHBd=0;
		dd.SHwHBd=0;
		// need to find count of CHX where X=F or Cl 
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			IAtom a=m.getAtom(i);			
			String s=a.getSymbol();
			
			if (s.equals("C")) {
				
				java.util.List ca=m.getConnectedAtomsList(a);
				
				boolean HaveClorF=false;
				boolean HaveH=false;
				
				if (a.getImplicitHydrogenCount()>0) HaveH=true;
				
				for (int j=0;j<ca.size();j++) {
					IAtom caj=(IAtom)ca.get(j);
					//if (b.getSymbol().equals("H")) HaveH=true;
					if (caj.getSymbol().equals("Cl")
							|| caj.getSymbol().equals("F"))
						HaveClorF = true;
				}

				if (HaveClorF && HaveH) {
					dd.numwHBd++;
					dd.SHwHBd+=HEState[i];
					
				}
				
			}
			
			
		} // end i for loop
		
		
	}
	
	
} // end class
