package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.*;

import org.openscience.cdk.Atom;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;

//import org.openscience.cdk.*;
//import org.openscience.cdk.config.ElementPTFactory;
//import org.openscience.cdk.config.IsotopeFactory;

import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.tools.periodictable.PeriodicTable;


public class EState {
	
	IAtomContainer m;
	double [] D;
	double [] DV;
	double [] DV2;
	int [][] Distance;
	double [] IS;
	double [] EState;
	double [] HEState;
	double [] KHE;
	
	
//	private HybridizationStateATMatcher hsam=null;
	Map <String,Integer>valences;
//	ElementPTFactory eptf = null;
//	IsotopeFactory elfac = null;
	
//	private AtomValenceDescriptor avd = null;
	
	void setupValencesMap() {
		
		//TODO- is this needed? might be already in CDK1.5
		valences = new HashMap();
        valences.put("H", new Integer(1));
        valences.put("Li", new Integer(1));
        valences.put("Be", new Integer(2));
        valences.put("B", new Integer(3));
        valences.put("C", new Integer(4));
        valences.put("N", new Integer(5));
        valences.put("O", new Integer(6));
        valences.put("F", new Integer(7));
        valences.put("Na", new Integer(1));
        valences.put("Mg", new Integer(2));
        valences.put("Al", new Integer(3));
        valences.put("Si", new Integer(4));
        valences.put("P", new Integer(5));
        valences.put("S", new Integer(6));
        valences.put("Cl", new Integer(7));
        valences.put("K", new Integer(1));
        valences.put("Ca", new Integer(2));
        valences.put("Ga", new Integer(3));
        valences.put("Ge", new Integer(4));
        valences.put("As", new Integer(5));
        valences.put("Se", new Integer(6));
        valences.put("Br", new Integer(7));
        valences.put("Rb", new Integer(1));
        valences.put("Sr", new Integer(2));
        valences.put("In", new Integer(3));
        valences.put("Sn", new Integer(4));
        valences.put("Sb", new Integer(5));
        valences.put("Te", new Integer(6));
        valences.put("I", new Integer(7));
        valences.put("Cs", new Integer(1));
        valences.put("Ba", new Integer(2));
        valences.put("Tl", new Integer(3));
        valences.put("Pb", new Integer(4));
        valences.put("Bi", new Integer(5));
        valences.put("Po", new Integer(6));
        valences.put("At", new Integer(7));
        valences.put("Fr", new Integer(1));
        valences.put("Ra", new Integer(2));
        valences.put("Cu", new Integer(2));
        valences.put("Mn", new Integer(2));
        valences.put("Co", new Integer(2));
        
        valences.put("Hg", new Integer(2));//added by TMM
        
        
	}
	
	
	public EState() {
//		avd = new AtomValenceDescriptor();
//		valencesTable = avd.valencesTable;
		
		this.setupValencesMap();
		
//		// add Hg to valences table:
//		valences.put("Hg", new Integer(2));
		
//		try {
//			elfac = IsotopeFactory.getInstance(new Molecule().getBuilder());
//			eptf = ElementPTFactory.getInstance();
//		} catch (Exception exc) {
//			exc.printStackTrace();
//		}
	}
	
	public void Calculate(IAtomContainer m,double [] D, double [] DV,double [] DV2,double [] KHE,int[][] Distance,double[] IS,double[] EState,double [] HEState ) {
		this.m=m;
		this.D=D;
		this.DV=DV;
		this.DV2=DV2;
		
		this.Distance=Distance;
		this.IS=IS;
		this.EState=EState;
		this.HEState=HEState;
		this.KHE=KHE;
		
		this.Initialize();
		
		
		this.CalculateDeltaValues();
		this.CalculateDV();
		
		
		this.CalcIntrinsicState();
		this.CalcKHE_Values();
		
		this.CalculateEState();
		this.CalculateHEState();
		//this.WriteEStates();
		
	}

	private void Initialize() {
		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			D[i] = 0.0;
			DV[i] = 0;
			IS[i] = 0.0;				
			EState[i] = 0;
			HEState[i] = 0;
			KHE[i] = 0;
			
		}

}
	
	
	
	private void CalculateEState() {
//		 System.out.println("");
		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			double sumDeltaIij = 0;

			for (int j = 0; j <= m.getAtomCount() - 1; j++) {
				if (i != j) {
					// System.out.println(i+"\t"+j+"\t"+IS[i]+"\t"+IS[j]+"\t"+Distance[i][j]);
					sumDeltaIij += (IS[i] - IS[j])
							/ Math.pow((double) Distance[i][j] + 1.0, 2.0);
				}
			}
			EState[i] = IS[i] + sumDeltaIij;
			// System.out.println(i+"\t"+IS[i]+"\t"+EState[i]);
		}
	}
	
	private void CalculateHEState() {
		// System.out.println("");
		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			double sumDeltaIij = 0;

			for (int j = 0; j <= m.getAtomCount() - 1; j++) {
				if (i != j) {
					// System.out.println(i+"\t"+j+"\t"+IS[i]+"\t"+IS[j]+"\t"+Distance[i][j]);
					sumDeltaIij += (KHE[j] + 0.2)
							/ Math.pow((double) Distance[i][j] + 1.0, 2.0);
				}
			}
			HEState[i] = KHE[i] + (KHE[i] + 0.2) + sumDeltaIij;
			// System.out.println(i+"\t"+HEState[i]);
		}

	}


	
	private void CalcIntrinsicState() {
		
		double Z=0, Zv=0;

		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			IAtom a = m.getAtom(i);
			
			
//			PeriodicTableElement e = eptf.getElement(a.getSymbol());
//			double N = Double.parseDouble(e.getPeriod());
			
			double N = (double)PeriodicTable.getPeriod(a.getSymbol());;
			
//			
//			try {
//			Zv = ((Integer) valences.get(a.getSymbol())).doubleValue(); // valence
//			Z = e.getAtomicNumber();
//			
//			} catch (Exception ex) {
//				System.out.println("Error getting valence for "+a.getSymbol());
//			}
//			//}
//			
//			double DV2 = DV[i] * (Z - Zv - 1);// need DeltaVprime on page 18
			// of estate book

//			System.out.println(DV2+"\t"+this.DV2[i]);
						
			//above is unnecessary since DV2 = DeltaVprime -TMM 7/15/08
			
			
			IS[i] = (Math.pow(2.0 / N, 2.0) * DV2[i] + 1) / D[i];
//			System.out.println(i+"\t"+a.getSymbol()+"\t"+IS[i]+"\t"+DV[i]);
			// System.out.println(a.getSymbol()+"\t"+IS[i]+"\t"+DV[i]);
		}

	}

	
	
	private void CalculateDeltaValues() {
		// determines D, DV2 for each specific atom

		// .....MODifY DV2 BASED ON ATOM TYPE
		
		//TODO: possibly redo this by using estate fragment strings to determine number of hydrogens
		
		ILoop: 
			for (int I = 0; I <= m.getAtomCount() - 1; I++) {

				DV2[I]=m.getBondOrderSum(m.getAtom(I));
				D[I]=m.getConnectedBondsCount(m.getAtom(I));
							
				String symbol = m.getAtom(I).getSymbol();
				int charge = m.getAtom(I).getFormalCharge();
				
				if (symbol.equals("C") || symbol.equals("Si")
					|| symbol.equals("Pb") || symbol.equals("Sn")) 
				{
					
					// TODO: I am not sure if following block is needed:
					if ((int) (D[I]) != (int) (0.5 + DV2[I])) {
						
						//used for carbons with fractional aromatic bond order
						int ID = (int) D[I];
						
						switch (ID) {
						case 1: { // 41
							
							if ((int) (DV2[I] + 0.5) != 2) {							
								DV2[I] = 3.0;
							}
							break;
						}
						
						case 2: { // 42
							if ((int) (DV2[I] + 0.5) != 3) {							
								DV2[I] = 4.0;
							}
							break;
						}
						
						case 3: { // 43
							if (DV2[I] >= 4) {							
								DV2[I] = 4.0;
							}
							break;
						}
						
						} // end switch ID
					}
					
				} else if (symbol.equals("O")) { // 50 (oxygen)
					
					if (DV2[I] > 1) {
						DV2[I] = 6.0;
					} else {
						DV2[I] = 5.0;
					}
															
				} else if (symbol.equals("N")) { // 60 (nitrogen)
					
					int ID = (int) D[I];
					// System.out.println("DV[3]="+DV[I]);
					
					switch (ID) {
					case 1: // 61-sp3
					{
						DV2[I] = DV2[I] + 2.;
						break;
					}
					
					case 2: // 62-sp2
					{
						if ((int) (DV2[I] + 0.5) == 2)
							DV2[I] = 4.;
						else
							DV2[I] = 5.;
						break;
					}
					
					case 3: // 63-sp
					{
						
						if ((int) (DV2[I] + 0.5) == 3 || (int) (DV2[I] + 0.5) == 4)
							DV2[I] = DV2[I] + 2.;
						else
							DV2[I] = 5.;
						
						break;
					}
					
					} // end switch ID
					
					
					
				} else if (symbol.equals("S")) { // 70 (sulfur)
					// ......SULFUR
					
					//DV2[I] = this.GetDVSulfur(DV2[I], I, D[I]);
										
					double h,Zv;
					
					Zv=6;
					
					if (DV2[I]==1) {
						h=1;
					} else {	
						h=0;
					}
					DV2[I]=Zv-h;
					
					
					
//				} else if (symbol.equals("Hg")) {
//					//TODO- can Hg have hydrogens attached?
//					
//					double Zv = 2;
//					double h = 0; // assume number of hydrogens = 0
//					
//					
//					DV2[I] = (Zv - h);
					
				
				} else if (symbol.equals("P") || symbol.equals("As")) {
					double h,Zv;
					
					Zv=5;
					
					if (DV2[I]==1) { 
						h=2;
					} else if (DV2[I]==2) {	
						h=1;
					} else if (DV2[I]==4) {	
						h=1;
					} else {
						h=0;						
					}
					
					DV2[I]=Zv-h;
					
					
				} else { // F,Cl,Br,I, Hg
					
//					IElement element = null;
					
					try {
//						element = elfac.getElement(symbol);
										
					double h = 0; // assume number of hydrogens = 0
					double Zv = ((Integer) valences.get(symbol)).doubleValue(); // valence
										
					DV2[I] = (Zv - h);

					
					} catch (Exception exc) {
						System.out.println("Error trying to get Zv for "+symbol);
//						exc.printStackTrace();
					}

					
				}
				
				// modify based on charge				
				DV2[I] += -charge; // alternatively could use the charge to modify the number of hydrogens
				
						
			} // end I for loop - 130
	
	
	
	} // end calculate delta values
	
	
	private void misccode() {
	
	}
	
	
	private void CalculateDV() {
		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			String symbol=m.getAtom(i).getSymbol();
			
//			IElement element = null;

			try {
//				element = elfac.getElement(symbol);
								
				double Z=0,Zv=0;
				
				Zv = ((Integer) valences.get(symbol)).doubleValue(); // valence
				Z = PeriodicTable.getAtomicNumber(symbol);
				
				DV[i]=DV2[i]/(Z - Zv - 1);
//				System.out.println(i+"\t"+DV[i]);
				
			} catch (Exception e) {
				System.out.println("Error trying to get valence for "+symbol);
				//e.printStackTrace();
			}
			
		}
	}
	

	private double GetDVSulfur(double DV, int I, double D) {
		double DVnew = 1;

		// right now D=number of neighbors of I, DV=number of bonds attached to
		// I

		double Z = 16; // atomic number

		double Zv = 0; // valence

		double h = 0;

		// DVnew=(Zv-h)/(Z-Zv-1);

		// System.out.println(D);

		switch ((int) D) {

		case 1: {
			if (DV == 1) { // we have SH
				h = 1;
				Zv = 6;
			} else if (DV == 2) { // we have S=X
				h = 0;
				Zv = 6;
			}
			break;
		}

		case 2: {
			// if (DV==2) {
			h = 0;
			Zv = 6;
			// }

			break;
		}

		case 3: {
			// if (DV==3) {
			// Zv=6;
			// h=1;
			// } else if (DV>3) {
			Zv = 6;
			h = 0;
			// }

			break;
		}

		case 4: {
			// if (DV==6) {
			// Zv=12;
			Zv = 6;
			h = 0;
			// }

			break;
		}

		}
		DVnew = (Zv - h);

		// System.out.println(DVnew);

		return DVnew;

	}

	private void CalcKHE_Values() {
		// this method calculates KHE and DV2 (deltaV without the denominator)

//		try {
//			eptf = ElementPTFactory.getInstance();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		double Z, Zv;

		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			IAtom a = m.getAtom(i);
//			PeriodicTableElement e = eptf.getElement(a.getSymbol());
			
			double N = (double)PeriodicTable.getPeriod(a.getSymbol());;
			
//			double N = Double.parseDouble(e.getPeriod());
			
			KHE[i] = (DV2[i] - D[i]) / (N * N);

		}

	}
	private void WriteEStates() {
		try {
			DecimalFormat d4 = new DecimalFormat("0.000000");

			FileWriter fw = new FileWriter("Estates.html");
			fw.write("<html>\n");

			fw.write("<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">\n");

			fw.write("<tr bgcolor=\"#D3D3D3\">\n");
			fw.write("\t<th>Atom</th>\n");
			fw.write("\t<th>Symbol</th>\n");
			fw.write("\t<th>Delta</th>\n");
			fw.write("\t<th>DeltaV</th>\n");
			fw.write("\t<th>Intrinsic State</th>\n");
			fw.write("\t<th>Estate</th>\n");			
			
			fw.write("</tr>\n");

			for (int i = 0; i <= m.getAtomCount() - 1; i++) {
				fw.write("<tr>\n");
				fw.write("\t<td>" + (i + 1) + "</td>\n");
				fw.write("\t<td>" + m.getAtom(i).getSymbol() + "</td>\n");
				fw.write("\t<td>" + d4.format(D[i]) + "</td>\n");
				fw.write("\t<td>" + d4.format(DV[i]) + "</td>\n");
				fw.write("\t<td>" + d4.format(IS[i]) + "</td>\n");
				fw.write("\t<td>" + d4.format(EState[i]) + "</td>\n");				
				fw.write("<tr>\n");
			}

			fw.write("</table>\n");
			fw.write("</html>\n");
			fw.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	public static void main(String[] args) {
		EState e =new EState();
		

	}
	
}
