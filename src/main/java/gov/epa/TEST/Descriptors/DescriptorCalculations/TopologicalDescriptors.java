package gov.epa.TEST.Descriptors.DescriptorCalculations;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

import java.util.ArrayList;
public class TopologicalDescriptors {

	IAtomContainer m;
	DescriptorData dd;
	double []D;
	double []DV;
	double [] EState;
	double [] IS;
	
	int []vdd;
	double A,B,C;
	int[][] Distance;
	int [] Eccentricity;
	int [] nge;
	int Ge;
	
	
	public void Calculate(IAtomContainer m,DescriptorData dd,double []D,double []DV,int [] vdd,int [][]Distance,double [] EState,double []IS) {
		this.m=m;
		this.dd=dd;
		this.D=D;
		this.DV=DV;
		this.vdd=vdd;
		this.Distance=Distance;
		this.EState=EState;
		this.IS=IS;
		
		
		this.A=(double)m.getAtomCount();
		this.B=(double)m.getBondCount();
		this.C=B-A+1;
		
		this.CalculateZM1();
		this.CalculateZM1V();
		this.CalculateZM2();
		this.CalculateZM2V();
		this.CalculateJ();
		this.CalculateJt();
		this.CalculateCentricIndices();
		this.CalculateICR();
		this.CalculateTIE();
		this.CalculateDeltaIDescriptors();
		this.Calculate_W();
		this.Calculate_WA();
	}
	
	private void CalculateDeltaIDescriptors() {
		dd.MAXDN=0;
		dd.MAXDP=0;
		dd.DELS=0;
		
		for (int i=0;i<=A-1;i++) {
			double sumDeltaIij = 0;
			
			for (int j = 0; j <= m.getAtomCount() - 1; j++) {
				if (i != j) {
					sumDeltaIij += (IS[i] - IS[j])
					/ Math.pow((double) Distance[i][j] + 1.0, 2.0);
				}
			}
			dd.DELS+=Math.abs(sumDeltaIij); // need ABS of deltas or DELS would be zero
			
						
			if (sumDeltaIij<0) {
				if (Math.abs(sumDeltaIij)>dd.MAXDN) dd.MAXDN=Math.abs(sumDeltaIij);
			} else if (sumDeltaIij>0) {
				if (Math.abs(sumDeltaIij)>dd.MAXDP) dd.MAXDP=Math.abs(sumDeltaIij);
			}
		}
		
	}
	
	
	
	
	private void CalculateEccentricity() {
		
		Eccentricity=new int[m.getAtomCount()];
		
		for (int i=0;i<=A-1;i++) {
			int max=0;
			for (int j=0;j<=A-1;j++) {
				if (Distance[i][j]>max) max=Distance[i][j];
			}
			Eccentricity[i]=max;
			//System.out.println(Eccentricity[i]);
			
		}
		
		
		
	}
	
	private void Calculate_W() {
		// Wiener Index
		// for CC(C)(C)C(C)CC, value = 63
		// which matches value on pg 129 of Bonchev, 1983
	
		
		dd.W=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			for (int j=0;j<=m.getAtomCount()-1;j++) {
				dd.W+=Distance[i][j];				
			}
		}
		dd.W/=2.0;
		
	
		 		 
	}
	
	private void Calculate_WA() {
		// Mean Wiener Index
		
		double denom=A*(A-1.0);
		
		dd.WA=2.0*dd.W/denom;
	}
	
	private void Calculate_nge() {
		// ng= number with same eccentricity
		Ge=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			if (Eccentricity[i]>Ge) Ge=Eccentricity[i];
		}
		
		nge=new int[Ge+1];
										
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			nge[Eccentricity[i]]++;
		}
		//System.out.println(G);
		
	}
	
	
	private void CalculateICR() {
		this.CalculateEccentricity();
		this.Calculate_nge();
		
		dd.ICR=0;

		double sum=0;
		 
		 for (int g=1;g<=Ge;g++) {
			 if (nge[g]>0) sum+=nge[g]/A*Math.log(nge[g]/A)/Math.log(2);
		 }
		 
		 dd.ICR=-sum;
				
		
	}
	
	private void CalculateZM1() {
		
		dd.ZM1=0;
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			dd.ZM1+=D[i]*D[i];
		}
		//System.out.println(dd.ZM1);
		
	}
	
	private void CalculateZM1V() {
		
		dd.ZM1V=0;
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			dd.ZM1V+=DV[i]*DV[i];
		}
		//System.out.println(dd.ZM1V);
		
	}
	
	private void CalculateZM2() {
		
		dd.ZM2=0;
		for (int b=0;b<m.getBondCount();b++) {
//			IAtom [] atoms=bonds[b].getAtoms();
			int i=m.getAtomNumber(m.getBond(b).getAtom(0));
			int j=m.getAtomNumber(m.getBond(b).getAtom(1));
			dd.ZM2+=D[i]*D[j];
			
		}
		//System.out.println(dd.ZM2);

	}
	
	private void CalculateZM2V() {
		
		dd.ZM2V=0;
		for (int b=0;b<m.getBondCount();b++) {
//			IAtom [] atoms=bonds[b].getAtoms();
			int i=m.getAtomNumber(m.getBond(b).getAtom(0));
			int j=m.getAtomNumber(m.getBond(b).getAtom(1));
			dd.ZM2V+=DV[i]*DV[j];
			
		}
		//System.out.println(dd.ZM2V);
	}

	private void CalculateTIE() {
		
		// See Voelkel, A. (1994). Computers Chem., 18, 1-4
		
		
		dd.TIE=0;
		//System.out.println("B="+B);
		
		double sum=0;
		for (int b=0;b<m.getBondCount();b++) {
			int i=m.getAtomNumber(m.getBond(b).getAtom(0));
			int j=m.getAtomNumber(m.getBond(b).getAtom(1));
			//sum+=Math.pow(EState[i]*EState[j],-0.5);
			// to match dragon:
			sum+=Math.pow(1+Math.exp(EState[i])*Math.exp(EState[j]),-0.5); // need different function since cant take sqrt of negative
			
			
		}
		dd.TIE=B/(C+1)*sum;
		
		//System.out.println(dd.ZM2V);
	}
	
	
	private void CalculateJ() {
		
		dd.J=0;
		for (int b=0;b<m.getBondCount();b++) {
			int i=m.getAtomNumber(m.getBond(b).getAtom(0));
			int j=m.getAtomNumber(m.getBond(b).getAtom(1));
			dd.J+=Math.pow(vdd[i]*vdd[j],-0.5);
			
		}
		dd.J*=B/(C+1);
		//System.out.println(dd.J);
	
	}
	
	private void CalculateJt() {
		
		dd.Jt=0;
		for (int b=0;b<m.getBondCount();b++) {
			int i=m.getAtomNumber(m.getBond(b).getAtom(0));
			int j=m.getAtomNumber(m.getBond(b).getAtom(1));
			dd.Jt+=Math.pow(vdd[i]/D[i]*vdd[j]/D[j],-0.5);
			
		}
		dd.Jt*=B/(C+1);
		//System.out.println(dd.Jt);
	
	}
	
	private void CalculateCentricIndices() {
		dd.BAC=0;
		dd.Lop=0;
		
		
		java.util.ArrayList al=new java.util.ArrayList();
		// compile list connected atoms for each atom
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			ArrayList al2=new ArrayList();
			al.add(al2);
			java.util.List ca=m.getConnectedAtomsList(m.getAtom(i));
			
			for (int j=0;j<ca.size();j++) {
				al2.add(new Integer(m.getAtomNumber(((IAtom)ca.get(j)))));
			}
			
		}
		
		
//		for (int i=0;i<=al.size()-1;i++) {
//			System.out.print(i+":\t");
//			
//			ArrayList al2=(ArrayList)al.get(i);
//			for (int j=0;j<=al2.size()-1;j++) {
//				System.out.print(al2.get(j)+"\t");
//			}
//			System.out.println("");
//		}
		
		int counter=0;
		
		while (true) {
			// figure out which atoms need to be removed:
			counter++;
			ArrayList removelist=new ArrayList();	
			
			int atomsremaining=0;
			
			for (int i=0;i<=al.size()-1;i++) {
				ArrayList al2=(ArrayList)al.get(i);
				
				if (al2.size()==1) {												
					removelist.add(new Integer(i));	
					al2.remove(0);
				}
				if (al2.size()>0) {
					atomsremaining++; //atoms left after current removal
				}
				
			}// end i loop
			
			if (atomsremaining==1) {
				dd.BAC+=1;
				double bob=1.0/(double)m.getAtomCount();
				dd.Lop-=bob*Math.log(bob)/Math.log(2);
			} 
			
			
			if (removelist.size()==0) {				
				
				if (counter==1) break; // to match dragon
				
				if (atomsremaining>0) { // to match dragon
					dd.BAC+=1;
					double bob=1.0/(double)m.getAtomCount();
					dd.Lop-=bob*Math.log(bob)/Math.log(2);
				}
								
				break;
			}
			
			dd.BAC+=Math.pow(removelist.size(),2);
			double bob=(double)removelist.size()/(double)m.getAtomCount();
			dd.Lop-=bob*Math.log(bob)/Math.log(2);			
			
//			System.out.println("removelist:");
//			for (int i=0;i<=removelist.size()-1;i++) {
//				System.out.println(removelist.get(i));
//			}
			
			//System.out.println(atomsremaining+"\t"+removelist.size());
			
			
			// now remove atoms on remove list:
			
			int removecount=0;
			for (int i=0;i<=al.size()-1;i++) {
				ArrayList al2=(ArrayList)al.get(i);
				
				for (int j=0;j<=al2.size()-1;j++) {
					int val=(Integer)al2.get(j);
					
					for (int k=0;k<=removelist.size()-1;k++) {
						int removeval=(Integer)removelist.get(k);
						if (removeval==val){						
							al2.remove(j);
							j=-1;
						}
					}
				}
				
			}// end i loop

			// print out remaining values:			
//			for (int i=0;i<=al.size()-1;i++) {
//				System.out.print(i+":\t");
//				
//				ArrayList al2=(ArrayList)al.get(i);
//				for (int j=0;j<=al2.size()-1;j++) {
//					System.out.print(al2.get(j)+"\t");
//				}
//				System.out.println("");
//			}
			
			int sum=0;
			for (int i=0;i<=al.size()-1;i++) {
						
				ArrayList al2=(ArrayList)al.get(i);
				sum+=al2.size();
			}
			
			if (sum==0) { // no atoms left				
				break;
			}
			
		} // end while true main loop
		
		
//		System.out.println(dd.BAC);
//		System.out.println(dd.Lop);

				
	}
	
//	private boolean CheckAnyBoundToIt(ArrayList al,int num){
//		boolean check=false;
//		
//		ArrayList nal=(ArrayList)al.clone();
//		
//		for (int i=0;i<=nal.size()-1;i++) {
//						
//			ArrayList al2=(ArrayList)nal.get(i);
//			for (int j=0;j<=al2.size()-1;j++) {
//				int val=(Integer)al2.get(j);
//				if (val==num) {
//					return true;
//				}
//			}
//			
//		}
//		
//		
//		
//		return check;
//	}
}
