package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.lang.reflect.Field;
import java.util.List;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

import org.openscience.cdk.Ring;

//import org.openscience.cdk.*;
//import org.openscience.cdk.interfaces.*;
//import org.openscience.cdk.tools.ValencyHybridChecker;

;

public class EStateFragmentDescriptor {

	DescriptorData dd;

	IAtomContainer m;
	IRingSet rs;
	
	String[] Fragment;

	

	public void Calculate(IAtomContainer m, DescriptorData dd, double[] EState,
			double[] HEState, IRingSet rs, String[] Fragment) {
		
		//TODO- code needs to be redone to account for cases in which
		// the hydrogen depleted structure causes the estate fragments to come out wrong
		// i.e. 053429-15-5	c1ccccc1CCCCCN(H)(H)(H)CL

		
//		for (int i=0;i<=m.getAtomCount()-1;i++){
//			System.out.println(i+"\t"+m.getAtom(i).getHydrogenCount());
//		}
		
		
		this.dd = dd;
		this.m = m;
		this.Fragment = Fragment;
		this.rs = rs;

		for (int ii = 0; ii <= dd.strES.length - 1; ii++) {
			CalcFragmentEStates(EState, dd.strES[ii]);
		}

		this.FindMissingFragments();

		// TODO- redesign it so it that determines which fragment each atom is then
		// match to existing list-ie loop over the atoms instead of over the
		// current fragment list

		this.CalcFragmentHEStates(HEState);

		this.CalcMinMax(EState, HEState);
		
//		EStateFragmentType e=new EStateFragmentType();		
//
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			try {
//			IAtomType iat=e.findMatchingAtomType(m,m.getAtom(i));
//			System.out.println(Fragment[i]+"\t"+iat.getAtomTypeName());
//			
//			} catch (Exception ex) {
//				
//			}
//			
//		}
		

	}

	private void CalcMinMax(double[] EState, double[] HEState) {

		dd.Hmax = -999999;
		dd.Gmax = -999999;
		dd.Hmaxpos = 0;

		dd.Hmin = 999999;
		dd.Gmin = 999999;
		dd.Hminneg = 0;

		for (int i = 0; i <= m.getAtomCount() - 1; i++) {

			if (!(Fragment[i] instanceof String)) {
				continue;
			}

			if (m.getAtom(i).getImplicitHydrogenCount() == 0) {
				HEState[i] = 0;
			}

			if (Fragment[i].indexOf("H") > -1
					&& Fragment[i].indexOf("Hg") == -1) {
				{

					// System.out.println(Fragment[i]+"\t"+HEState[i]);
					if (HEState[i] > dd.Hmax) {
						dd.Hmax = HEState[i];
						if (HEState[i] > 0) {
							dd.Hmaxpos = HEState[i];
						}
					}

					if (HEState[i] < dd.Hmin) {
						dd.Hmin = HEState[i];
						if (HEState[i] < 0) {
							dd.Hminneg = HEState[i];
						}
					}
				}
			}

			if (EState[i] > dd.Gmax)
				dd.Gmax = EState[i];
			if (EState[i] < dd.Gmin)
				dd.Gmin = EState[i];

		}// end i for loop

		if (dd.Hmax == -999999)
			dd.Hmax = 0;
		if (dd.Gmax == -999999)
			dd.Gmax = 0;
		if (dd.Hmin == 999999)
			dd.Hmin = 0;
		if (dd.Gmin == 999999)
			dd.Gmin = 0;

	}

	private void CalcFragmentEStates(double[] EState, String strFragment) {

		int FragmentCharge;

		try {

			// initialize values:
			Field myField = dd.getClass().getField(strFragment);
			myField.setDouble(dd, 0.0);

			Field myField2 = dd.getClass().getField(strFragment + "_acnt");
			myField2.setInt(dd, 0);

			String strFragment2 = strFragment
					.substring(1, strFragment.length());
			// System.out.println(strFragment);

			int pos = -1;
			String s = "-";

			int NumSingleBonds = 0; // number of single bonds in search fragment
			int NumDoubleBonds = 0; // number of double bonds in search fragment
			int NumTripleBonds = 0; // number of triple bonds in search fragment
			int NumAromaticBonds = 0; // number of aromatics bonds in search
			// fragment
			int NumAromaticBondsTotal = 0; // need to distinguish between SaaNH
			// and SaaN

			int NumH; // number of hydrogens in search fragment

			do {
				pos++;
				s = strFragment2.substring(pos, pos + 1);

				if (s.equals("s")) {
					NumSingleBonds++;
				}

				if (s.equals("d")) {
					NumDoubleBonds++;
				}

				if (s.equals("t")) {
					NumTripleBonds++;
				}

				if (s.equals("a")) {
					NumAromaticBonds++;
				}

				// System.out.println("s="+s);

			} while (s.equals("s") || s.equals("d") || s.equals("t")
					|| s.equals("a"));

			String element = strFragment2.substring(pos, pos + 1);
			String element2 = "";

			if (strFragment2.length() >= pos + 2)
				element2 = strFragment2.substring(pos + 1, pos + 2);

			if (!element2.equals("H") && !element2.equals("m")
					&& !element2.equals("p"))
				element += element2;

			NumH = 0;

			String bob = strFragment2.substring(pos + 1, strFragment2.length());

			if (bob.length() > 0) {
				if (bob.indexOf("H") > -1) {
					String bob2 = bob.substring(bob.indexOf("H"), bob.length());
					if (bob2.length() == 1)
						NumH = 1;
					else if (bob2.length() == 2) {
						// System.out.println(strFragment2+"\t"+bob2.substring(1,2));

						if (bob2.substring(1, 2).equals("2"))
							NumH = 2;
						else if (bob2.substring(1, 2).equals("3"))
							NumH = 3;
						else if (bob2.substring(1, 2).equals("4"))
							NumH = 4;
						else
							NumH = 1;

					}

				}
			}

			// System.out.println(strFragment2+"\t"+NumH);

			if (strFragment.indexOf("aa") > -1) {
				if (element.equals("N")) {
					if (strFragment.equals("SaaNH")) {
						NumAromaticBondsTotal = 2;
					} else if (strFragment.equals("SaaN")) {
						NumAromaticBondsTotal = 3;
					} else if (strFragment.equals("SsaaN")) {
						NumAromaticBondsTotal = 2;
					} else if (strFragment.equals("SsaaNH")) {
						NumAromaticBondsTotal = 3;
					} else if (strFragment.equals("SssaaNH")) {
						NumAromaticBondsTotal = 2;
					} else if (strFragment.equals("SssaaN")) {
						NumAromaticBondsTotal = 3;
					}

					// else if (strFragment2.equals("daaN")) {
					// NumAromaticBondsTotal = 3; // five bonds to N
					//		
					// } else if (strFragment2.equals("aaaN")) {
					// NumAromaticBondsTotal = 3;
					// }
				} else {

				}
			}

			String last = strFragment2.substring(strFragment2.length() - 1,
					strFragment2.length());

			if (last.equals("m")) {
				FragmentCharge = -1;
			} else if (last.equals("p")) {
				FragmentCharge = 1;
			} else
				FragmentCharge = 0;

			// System.out.println(strFragment+"\t"+FragmentCharge);

			// System.out.println(strFragment+"\t"+element+"\t"+NumSingleBonds+"\t"+NumDoubleBonds+"\t"+NumTripleBonds+"\t"+NumAromaticBonds);
			// System.out.println("\n"+NumSingleBonds +"\t"+NumDoubleBonds
			// +"\t"+NumTripleBonds +"\t"+NumAromaticBonds+"\n" );

			for (int i = 0; i <= m.getAtomCount() - 1; i++) {
				
				if (Fragment[i] instanceof String) continue; // skip calcs if already assigned
				
				int NumSingleBonds2 = 0;
				int NumDoubleBonds2 = 0;
				int NumTripleBonds2 = 0;
				int NumAromaticBonds2 = 0;
				int NumAromaticBondsTotal2 = 0;

				IAtom a = m.getAtom(i);

				// System.out.println(a.getSymbol()+"\t"+a.getFormalCharge());

				if (!a.getSymbol().equals(element))
					continue;

				// if (strFragment.equals("SaaNH"))
				// System.out.println((i+1)+"\t"+a.getSymbol());

//				IBond[] bonds = m.getConnectedBonds(a);

				List ca = m.getConnectedAtomsList(a);

				for (int j = 0; j < ca.size() ; j++) {

					IAtom caj=(IAtom)ca.get(j);
					IBond b = m.getBond(a, caj);

					if (a.getFlag(CDKConstants.ISAROMATIC)
							&& caj.getFlag(CDKConstants.ISAROMATIC)) {

						boolean SameRing = InSameAromaticRing(m, a,caj, rs);

						if (SameRing) {
							NumAromaticBonds2++;
							if (element.equals("N")) {
								if (b.getOrder().equals(IBond.Order.SINGLE))
									NumAromaticBondsTotal2++;
								if (b.getOrder().equals(IBond.Order.DOUBLE))
									NumAromaticBondsTotal2 = NumAromaticBondsTotal2 + 2;
							}
						} else {
							if (b.getOrder().equals(IBond.Order.SINGLE))
								NumSingleBonds2++;
							if (b.getOrder().equals(IBond.Order.DOUBLE))
								NumDoubleBonds2++;
							if (b.getOrder().equals(IBond.Order.TRIPLE))
								NumTripleBonds2++;
						}

					} else {

						if (b.getOrder().equals(IBond.Order.SINGLE))
							NumSingleBonds2++;
						if (b.getOrder().equals(IBond.Order.DOUBLE))
							NumDoubleBonds2++;
						if (b.getOrder().equals(IBond.Order.TRIPLE))
							NumTripleBonds2++;
					}
				}

				/*
				 * if (strFragment.equals("SaaNH")) { if (i+1==1) {
				 * System.out.println(NumSingleBonds2);
				 * System.out.println(a.getFlag(CDKConstants.ISAROMATIC)); } }
				 */

				// if (strFragment.equals("SaaNH"))
				// System.out.println("\n"+NumSingleBonds2);
				// System.out.print(NumSingleBonds2+"\t"+NumDoubleBonds2+"\t"+NumTripleBonds2+"\t"+NumAromaticBonds2);
				// System.out.println(strFragment+"\t"+FragmentCharge+"\t"+a.getFormalCharge());
				/*
				 * if (strFragment.equals("SaaaC")) {
				 * System.out.println("i="+(i+1));
				 * System.out.println(NumSingleBonds+"\t"+NumSingleBonds2);
				 * System.out.println(NumDoubleBonds+"\t"+NumDoubleBonds2);
				 * System.out.println(NumTripleBonds+"\t"+NumTripleBonds2);
				 * System.out.println(NumAromaticBonds+"\t"+NumAromaticBonds2); }
				 */
				
				
//				 System.out.println("NumSingleBonds\t"+NumSingleBonds+"\t"+NumSingleBonds2);
//				 System.out.println("NumDoubleBonds\t"+NumDoubleBonds+"\t"+NumDoubleBonds2);
//				 System.out.println("NumTripleBonds\t"+NumTripleBonds+"\t"+NumTripleBonds2);
//				 System.out.println("NumAromaticBonds\t"+NumAromaticBonds+"\t"+NumAromaticBonds2); 
		
				boolean match = false;
				if (NumSingleBonds == NumSingleBonds2
						&& NumDoubleBonds == NumDoubleBonds2
						&& NumTripleBonds == NumTripleBonds2
						&& NumAromaticBonds == NumAromaticBonds2
						&& FragmentCharge == a.getFormalCharge()) {

					if (strFragment.equals("SaaN")
							|| strFragment.equals("SaaNH")
							|| strFragment.equals("SsaaN")
							|| strFragment.equals("SsaaNH")
							|| strFragment.equals("SssaaNH")
							|| strFragment.equals("SssaaN")) { // only check
						// total
						// numaromaticbonds
						// for these
						// frags

						if (NumAromaticBondsTotal == NumAromaticBondsTotal2) {
							match = true;
						}
					} else {
						match = true;
					}

				}

				if (match) {
					double value = myField.getDouble(dd);
					myField.setDouble(dd, value + EState[i]);
					Fragment[m.getAtomNumber(a)] = strFragment;

					// need to do this because CDK doesnt do num H's right:
					m.getAtom(i).setImplicitHydrogenCount(NumH);

//					 System.out.println(i+"\t"+strFragment+"\t"+NumH);
					// System.out.println(strFragment+" found for atom"+(i+1));

					int value2 = myField2.getInt(dd);
					myField2.setInt(dd, value2 + 1);

				}

			} // end i atom for loop

			// System.out.println(strFragment+"\t"+myField2.getInt(this)+"\t"+myField.getDouble(this));

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void FindMissingFragments() {
		try {

			for (int i = 0; i <= m.getAtomCount() - 1; i++) {

				if (Fragment[i] instanceof String)
					continue;

				int NumSingleBonds2 = 0;
				int NumDoubleBonds2 = 0;
				int NumTripleBonds2 = 0;
				int NumAromaticBonds2 = 0;
				int NumAromaticBondsTotal2 = 0;
				int HydrogenCount = 0;

				IAtom a = m.getAtom(i);

				String element = a.getSymbol();

				// System.out.println(a.getSymbol()+"\t"+a.getFormalCharge());

				// if (strFragment.equals("SaaNH"))
				// System.out.println((i+1)+"\t"+a.getSymbol());

//				IBond[] bonds = m.getConnectedBonds(a);

				List ca = m.getConnectedAtomsList(a);

				for (int j = 0; j <ca.size(); j++) {

					IAtom caj=(IAtom)ca.get(j);
					IBond b = m.getBond(a, caj);

					if (a.getFlag(CDKConstants.ISAROMATIC)
							&& caj.getFlag(CDKConstants.ISAROMATIC)) {

						boolean SameRing = InSameAromaticRing(m, a,caj, rs);

						if (SameRing) {
							NumAromaticBonds2++;
							if (element.equals("N")) {
								if (b.getOrder().equals(IBond.Order.SINGLE))
									NumAromaticBondsTotal2++;
								if (b.getOrder().equals(IBond.Order.DOUBLE))
									NumAromaticBondsTotal2 = NumAromaticBondsTotal2 + 2;
							}
						} else {
							if (b.getOrder().equals(IBond.Order.SINGLE))
								NumSingleBonds2++;
							if (b.getOrder().equals(IBond.Order.DOUBLE))
								NumDoubleBonds2++;
							if (b.getOrder().equals(IBond.Order.TRIPLE))
								NumTripleBonds2++;
						}

					} else {

						if (b.getOrder().equals(IBond.Order.SINGLE))
							NumSingleBonds2++;
						if (b.getOrder().equals(IBond.Order.DOUBLE))
							NumDoubleBonds2++;
						if (b.getOrder().equals(IBond.Order.TRIPLE))
							NumTripleBonds2++;
					}
				}

				// assign frag here
				Fragment[i] = "S";

				for (int j = 0; j <= NumTripleBonds2 - 1; j++) {
					Fragment[i] += "t";
				}

				for (int j = 0; j <= NumDoubleBonds2 - 1; j++) {
					Fragment[i] += "d";
				}

				for (int j = 0; j <= NumSingleBonds2 - 1; j++) {
					Fragment[i] += "s";
				}

				for (int j = 0; j <= NumAromaticBonds2 - 1; j++) {
					Fragment[i] += "a";
				}

				Fragment[i] += element;
				
				// add number of hydrogens here-
//				ValencyHybridChecker vhc=new ValencyHybridChecker();
//				
//				int numH=vhc.calculateNumberOfImplicitHydrogens(a,m);
				
				int numH=a.getImplicitHydrogenCount();//TODO- do we a class to set this?
				
				
				if (numH>=1) Fragment[i]+="H";
				
				if (numH>1) Fragment[i]+=numH;
				

				if (a.getFormalCharge() == 1) {
					Fragment[i] += "p";
				} else if (a.getFormalCharge() == -1) {
					Fragment[i] += "m";
				}

				Fragment[i] += "*";// flag it as unassigned fragment

				// for now set #H = 0 TODO- add routine to figure out #Hs for
				// possible valence states
				m.getAtom(i).setImplicitHydrogenCount(0);

			} // end i atom for loop
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


		
	public static boolean InSameAromaticRing(IAtomContainer m, IAtom atom1,
			IAtom atom2, IRingSet rs) {
		boolean SameRing = false;

		for (int i = 0; i <= rs.getAtomContainerCount() - 1; i++) {
			Ring r = (Ring) rs.getAtomContainer(i);
			
			if (!r.getFlag(CDKConstants.ISAROMATIC))
				continue;

			if (r.contains(atom1) && r.contains(atom2)) {
				return true;
			}
			

		} // end ring for loop

		return false;
	}

	private void CalcFragmentHEStates(double[] HEState) {

		dd.SHsOH = 0;
		dd.SHother = 0;
		dd.SHdNH = 0;
		dd.SHsSH = 0;
		dd.SHsNH2 = 0;
		dd.SHssNH = 0;
		dd.SHtCH = 0;
		dd.SHCHnX = 0;
		// SHCHnX (a CH or CH2 group w/F or Cl bonded to carbon)

		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			if (!(Fragment[i] instanceof String)) {
				continue;
			}
			if (Fragment[i].equals("SsOH")) {
				dd.SHsOH += HEState[i];
			} else if (Fragment[i].equals("SdNH")) {
				dd.SHdNH += HEState[i];
			} else if (Fragment[i].equals("SsSH")) {
				dd.SHsSH += HEState[i];
			} else if (Fragment[i].equals("SsNH2")) {
				dd.SHsNH2 += HEState[i];
			} else if (Fragment[i].equals("SssNH")
					|| Fragment[i].equals("SaaNH")) {
				dd.SHssNH += HEState[i];
			} else if (Fragment[i].equals("StCH")) {
				dd.SHtCH += HEState[i];
			} else if (Fragment[i].equals("SsCH3")
					|| Fragment[i].equals("SssCH2")
					|| Fragment[i].equals("SsssCH")
					|| Fragment[i].equals("SdsCH")) {

				List ca = m.getConnectedAtomsList(m.getAtom(i));

				boolean b = false;
				for (int j = 0; j <ca.size(); j++) {
					IAtom a = (IAtom)ca.get(j);
					String symbol = a.getSymbol();
					if (symbol.equals("F") || symbol.equals("Cl")) {
						b = true;
						break;
					}
				}
				if (b)
					dd.SHCHnX += HEState[i];
				else
					dd.SHother += HEState[i];

			} else if (Fragment[i].indexOf("H") > -1
					&& Fragment[i].indexOf("Hg") == -1) {
				// System.out.println(Fragment[i]);
				dd.SHother += HEState[i];
			} else {

			}

		}

	}

}
