package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.util.*;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

import org.openscience.cdk.Ring;


//import org.openscience.cdk.*;
//import org.openscience.cdk.interfaces.*;



// book used for verification of values:
// Bonchev 1983, Information Theoretic Indices for characterization of chemical
// structures




public class InformationContentDescriptor {
	
	int [][] Distance;
	IAtomContainer m;
	IRingSet rs;
	DescriptorData d;
	
	double A,B;
	LinkedList[] paths;
	int [] fg;	
	int G=0;// maximum distance value (i.e. the topological diameter)
	
	int [][] EdgeDistance;	
	int [] fgedge;	
	int Gedge=0;// maximum distance value (i.e. the topological diameter)

	int [] ng; // number with each vertex distance degree
	int [] vdd;// vertex distance degree array
	int Gvdd=0; // maximum vertex distance degreee
	
	int Gvertex;
	int [] Fg; // number with each vertex degree
	
	
	double [] D;
	double [] DV;
	double []S;
	double []SV;
	
	// need connectivity index of order 1 to calculate eim
	
	
	
	public  void Calculate(IAtomContainer m,DescriptorData d,int [][] Distance,IRingSet rs,LinkedList[] paths,double [] D,double []DV,double []S,double []SV,int [] vdd) {
		this.Distance=Distance;
		this.m=m;
		this.d=d;
		this.rs=rs;
		this.paths=paths;
		this.D=D;
		this.DV=DV;
		this.S=S;
		this.SV=SV;
		this.vdd=vdd;
		
					
		this.A=(double)m.getAtomCount();
		this.B=(double)m.getBondCount();
		 //System.out.println("A="+A);
		 //System.out.println("B="+B);
				
		
		this.Calculate_fg();
		this.CalculateEdgeDistanceMatrix();
		this.Calculate_fgedge();
		this.Calculate_ng();
		this.Calculate_Fg();
		
		
		Calculate_iadje();
		Calculate_iadjem();
		Calculate_iadjm();
		Calculate_iadjmm();
		Calculate_ieadje();
		Calculate_ieadjem();
		Calculate_ieadjm();
		Calculate_ieadjmm();
		Calculate_ide();
		Calculate_idem();				
		Calculate_idm();
		Calculate_idmm();
		Calculate_iede();
		Calculate_iedem();
		Calculate_We();
		Calculate_iedm();
		Calculate_iedmm();
		Calculate_tvc();
		
		Calculate_iddem();
		Calculate_iddmm();
		Calculate_ivdem();
		Calculate_ivdmm();

		Calculate_cycle_information_indices();
		
		//System.out.println("calculating tti...");
		
		d.tti=Calculate_tti2(1,1,D,S); // Calculate_tti could be modified to calculate tti and ttvi simultaneously to save time		
		d.ttvi=Calculate_tti2(1,1,DV,SV);
		
		ArrayList MasterList=FindSymmetryClasses();
		Calculate_ic(MasterList);
		this.Calculate_si(MasterList);
//		this.Calculate_k0(MasterList);
		this.Calculate_maxic();
		this.Calculate_I();
		this.Calculate_ssi();
		this.Calculate_R();
		this.Calculate_ib();
		this.Calculate_eim();
		
		
		/*
		for (int i=0;i<=d.strIC.length-1;i++) {
	    	try {
	    		Field f=d.getClass().getField(d.strIC[i]);
	    		f.setDouble(d,-9999);
	    		
	    		Method method=this.getClass().getMethod("Calculate_"+d.strIC[i],null);
	    		System.out.println("Calculate_"+d.strIC[i]+"();");
	    		
	    		method.invoke(this,null);
	    		
	    	} catch (Exception e) {
	    		System.out.println(e);
	    	}
	    }
		*/
		
		
				
	}
	
	private void Calculate_ic(ArrayList MasterList) {
		d.ic=0;
		for (int j = 0; j <= MasterList.size() - 1; j++) {
			ArrayList al = (ArrayList) MasterList.get(j);

			double count = Double.parseDouble((String) al.get(1));
//			System.out.println(count);
			d.ic += count * Log(2,count);

			// System.out.println(al.get(0)+"\t"+al.size());
		}
		
	}
	
	private void Calculate_I() {
		// note I=k0 so k0 was omitted
		d.I=d.si*A;
				
	}
	
	private void Calculate_maxic() {
		d.maxic=A*Log(2,A);
	}
	
	private void Calculate_ssi() {
		// note n=sum ng = A
		d.ssi=d.I/d.maxic;
	}
	
	private void Calculate_R() {
		d.R=1-d.ssi;
		
	}
	private void Calculate_ib() {
		d.ib=0;
		
		double [] bondcount=new double [5];
		
		for (int i=1;i<=4;i++) {
			bondcount[i]=0;
		}
		
				
		for (int i=0;i<m.getBondCount();i++) {			
			IAtom a1=m.getBond(i).getAtom(0);
			IAtom a2=m.getBond(i).getAtom(1);
			
			if (a1.getFlag(CDKConstants.ISAROMATIC)
					&& a2.getFlag(CDKConstants.ISAROMATIC)) {

				boolean SameRing = EStateFragmentDescriptor.InSameAromaticRing(m,a1,a2, rs);
				if (SameRing) {
					bondcount[4]++;
				} else {
					bondcount[(int)m.getBond(i).getOrder().numeric()]++;
				}
			} else {
				bondcount[(int)m.getBond(i).getOrder().numeric()]++;
			}
			
		}
		
		double sum=0;
		for (int i=1;i<=4;i++) {
			//System.out.println(i+"\t"+bondcount[i]);
			if (bondcount[i]>0) {
				sum+=bondcount[i]*Log(2,bondcount[i]);
			}
		}
		d.ib=B*Log(2,B)-sum;
		
		
		
	}
	
	private void Calculate_eim() {
		
		
		double sum=0;
		
		for (int i=0;i<m.getBondCount();i++) {			
			IAtom a1=m.getBond(i).getAtom(0);
			IAtom a2=m.getBond(i).getAtom(1);
						
			double d1=D[m.getAtomNumber(a1)];
			double d2=D[m.getAtomNumber(a2)];
			
			double term=Math.pow(d1*d2,-0.5)/d.x1;
			sum+=term*Log(2,term);
						
		}
		d.eim=-sum;
	}
	
	
	private void Calculate_si(ArrayList MasterList) {
		// note si=k0/A
		
		d.si=0;
		double n=0;
		
		for (int j = 0; j <= MasterList.size() - 1; j++) {
			ArrayList al = (ArrayList) MasterList.get(j);

			double count = Double.parseDouble((String) al.get(1));
			//System.out.println(count);
			n+= count;

		}
		
//		System.out.println("n="+n);
//		System.out.println("A="+A);
		
		
		for (int j = 0; j <= MasterList.size() - 1; j++) {
			ArrayList al = (ArrayList) MasterList.get(j);

			double ng = Double.parseDouble((String) al.get(1));
//			System.out.println(ng);
			d.si += ng/n * Log(2,ng/n);
//			d.si += ng/n * Log(10,ng/n);
			

		}
		d.si=-d.si;
		// note si=k0/A
		
		
		
	}
	
	
//	private void Calculate_k0(ArrayList MasterList) {
//		d.k0=0;
//		
//		for (int j = 0; j <= MasterList.size() - 1; j++) {
//			ArrayList al = (ArrayList) MasterList.get(j);
//
//			double count = Double.parseDouble((String) al.get(1));
//			
//			double p = count / (double) m.getAtomCount();
//								
//			d.k0 += p * Log(2,p); // note- to match MDL QSAR need log10 not log2
//			
//			
//			// System.out.println(al.get(0)+"\t"+al.size());
//		}
//		d.k0 *= -(double) m.getAtomCount();
//		
//	}
	
	
	private ArrayList FindSymmetryClasses() {
		double tol=1e-6;
		
		ArrayList MasterList = new ArrayList();
				
		double [] SYM=SV;
				
//		 get things started:
		ArrayList a = new ArrayList();		
		a.add(SYM[0] + "");		
		a.add(1 + "");
		
		MasterList.add(a);
		
		
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			System.out.println(SV[i]);			
//		}
		
		iloop:
		for (int i=1;i<=m.getAtomCount()-1;i++) {
			for (int j = 0; j <= MasterList.size() - 1; j++) {

				ArrayList al = (ArrayList) MasterList.get(j);
				double sv = Double.parseDouble((String) al.get(0));
				
				if (Math.abs(SYM[i] - sv) < tol) {
					
					int count = Integer.parseInt((String) al.get(1));
					count++;
					al.set(1, count + "");
					continue iloop;					
				}
			} // end j for loop

			// no match:
			ArrayList a2 = new ArrayList();			
			a2.add(SYM[i] + "");			
			a2.add(1 + "");
			MasterList.add(a2);
			
		}
					
		return MasterList;
	
		
	}
	
	private void Calculate_iddem() {
		//Mean information content on the distance degree equality (iddem)
		d.iddem=0;
		        
		// value equal to 1.4591 for CC2C1CCC12- matches Dragon 5.0
		
		double sum=0;
		 
		 for (int g=1;g<=Gvdd;g++) {
			 if (ng[g]>0) sum+=ng[g]/A*Log(2,ng[g]/A);
		 }
		 
		 d.iddem=-sum;
				
		
	}
	
	private void Calculate_iddmm() {

		//Mean information content on the distance degree magnitude (iddmm)
		d.iddmm=0;
		
       // 	value equal to 2.5625 for CC2C1CCC12- matches Dragon 5.0
		
		double sum=0;
		 
		for (int i=0;i<=A-1;i++) {
			sum+=vdd[i]/(2*d.W)*this.Log(2,vdd[i]/(2*d.W));
		}
		// alternate way:
		//for (int g=1;g<=Gvdd;g++) {
		//	if (ng[g]>0) sum+=ng[g]*g/(2*d.W)*Log(2,g/(2*d.W));
	    //}
		
		 d.iddmm=-sum;
		
		
	}
	
	private void Calculate_ivdem() {

		//Mean Information Content of the Vertex Degree Equality (ivdem)
		d.ivdem=0;
		// value = 1.4591 for CC2C1CCC12 (molecule g2 in bonchev book)- matches dragon 5.0
		double sum=0;
		 
		 for (int g=1;g<=Gvertex;g++) {
			 if (Fg[g]>0) sum+=Fg[g]/A*Log(2,Fg[g]/A);
		 }
		 
		 d.ivdem=-sum;
				
		
		
	}
	
	private void Calculate_ivdmm() {

		//Mean Information Content of the Vertex Degree magnitude (ivdmm)
		d.ivdmm=-9999;
		// value = 2.5027 for CC2C1CCC12; matches Dragon 5.0
		double sum=0;
		 
		 for (int g=1;g<=Gvertex;g++) {
			 if (Fg[g]>0) sum+=Fg[g]*g/(2*B)*Log(2,g/(2*B));
		 }
		 
		 d.ivdmm=-sum;
		
		
	}


	
	
	
	private void Calculate_fg() {
		// fg = count with distance equal to g
		G=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			for (int j=0;j<=i-1;j++) {// use lower triangular distance submatrix		
				if (Distance[i][j]>G) G=Distance[i][j];
			}
		}
		
		fg=new int[G+1];
										
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			for (int j=0;j<=i-1;j++) {// use lower triangular distance submatrix
				fg[Distance[i][j]]++;		
			}
		}
		//System.out.println(G);
		
	}
	
	private void Calculate_fgedge() {
		// fgedge = count with edge distance equal to g
		Gedge=0;
		
		for (int i=0;i<=EdgeDistance.length-1;i++) {
			for (int j=0;j<=i-1;j++) {// use lower triangular distance submatrix																	
				if (EdgeDistance[i][j]>Gedge) Gedge=EdgeDistance[i][j];
			}
		}
		
		fgedge=new int [Gedge+1];		
				
		for (int i=0;i<=EdgeDistance.length-1;i++) {
			for (int j=0;j<=i-1;j++) {// use lower triangular distance submatrix															
				fgedge[EdgeDistance[i][j]]++;		
			}
		}
		
		
	}
	
	private void Calculate_Fg() {
		// Fg = count with vertex degree equal to g
		Gvertex=0;
		// first determine max vertex degree 
		for (int i=0;i<=m.getAtomCount()-1;i++) {					
			if (D[i]>Gvertex) Gvertex=(int)D[i];			
		}
		
		//now we can dimension the Fg array:
		Fg=new int[Gvertex+1];
						
		for (int i=0;i<=m.getAtomCount()-1;i++) {			
			Fg[(int)D[i]]++;
		}
				
	}
	
	
	
	private void Calculate_ng()
	// ng is the number atoms with vertex distance degree equal to g	
	{
		this.Gvdd=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {		
			if (vdd[i]>Gvdd) Gvdd=vdd[i];			
		}
		
		ng=new int[Gvdd+1];
				
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			ng[vdd[i]]++;					
		}
		
	}
	
	
	private void CalculateEdgeDistanceMatrix()  {
		
		LinkedList ll=new LinkedList();
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			for (int j=0;j<=i;j++) {				
				if (Distance[i][j]==1) {
					LinkedList ll2=new LinkedList();
					ll2.add(i+"");
					ll2.add(j+"");
					ll.add(ll2);
					//System.out.println(ll.size()+"\t"+i+"\t"+j);
				}
			}
		}
		
		EdgeDistance=new int [ll.size()][ll.size()];
		
		
		for (int i=0;i<=ll.size()-1;i++) { // loop over edges
			
			LinkedList ll2=(LinkedList)ll.get(i);
			
			String strone=(String)ll2.get(0);
			String strtwo=(String)ll2.get(1);
			
			int one=Integer.parseInt(strone);
			int two=Integer.parseInt(strtwo);
			
			for (int j=0;j<=ll.size()-1;j++) {
				LinkedList ll3=(LinkedList)ll.get(j);
				String strthree=(String)ll3.get(0);
				String strfour=(String)ll3.get(1);
												
				int three=Integer.parseInt(strthree);
				int four=Integer.parseInt(strfour);
				
				int dist1=Distance[one][three];
				int dist2=Distance[one][four];
				int dist3=Distance[two][three];
				int dist4=Distance[two][four];
				
				
				//if (i==0) System.out.println(dist1+"\t"+dist2+"\t"+dist3+"\t"+dist4);
												
				if (one==three && two==four ) {
					EdgeDistance[i][j]=0;
				} else {
					int mindist=9999;					
					if (dist1<mindist) mindist=dist1;
					if (dist2<mindist) mindist=dist2;
					if (dist3<mindist) mindist=dist3;
					if (dist4<mindist) mindist=dist4;
					EdgeDistance[i][j]=mindist+1; // page 129 todeschini book	
				}
																
			}
			
			
		}
		
		/*for (int i=0;i<=ll.size()-1;i++) { // loop over edges			
			for (int j=0;j<=ll.size()-1;j++) {
				System.out.print(EdgeDistance[i][j]+"\t");			
			}
			System.out.println("");
		}
			*/						
		
	}
	
	

	
	private void Calculate_We() {
		// Edge Wiener Index
		
		//d.We=d.W-A*(A-1.0)/2.0; // wrong for cyclic?
		
		//System.out.println(d.We);
		
		/*
		d.We=0;
		for (int i=0;i<=EdgeDistance.length-1;i++) { // loop over edges			
			for (int j=0;j<=EdgeDistance.length-1;j++) {
				d.We+=EdgeDistance[i][j];				
			}			
		}
		d.We/=2.0;
		
		System.out.println(d.We);
		*/
		
		d.We=0;
		
		for (int g=1;g<=Gedge;g++) {
			 d.We+=fgedge[g]*g;
		 }
		
		//System.out.println(d.We);
		
		
		
	}
	
	private void Calculate_tvc() {
		// Total Vertex Cyclicity
				
		d.tvc=0;
		
		for (int i = 0; i <= rs.getAtomContainerCount() - 1; i++) {
			Ring r = (Ring) rs.getAtomContainer(i);
			d.tvc+=r.getAtomCount();
			
			
		}

	}
	
	private void Calculate_cycle_information_indices() {
		
		// this method calculates the following:
		//Total Information on the Vertex Cycle Matrix Equality (icyce)
		//Mean Information on the Vertex Cycle Matrix Equality (icycem)
		//Total Information on the Vertex Cycle Matrix Magnitude (icycm)
		//Mean Information on the Vertex Cycle Matrix Magnitude (icycmm)
				
		// for CC2C1CCC12, icycem = 0.9183 
		// which matches value on pg 126 of Bonchev, 1983
		
		// for CC2C1CCC12, icycmm = 3.5850 
		// which matches value on pg 126 of Bonchev, 1983
		
		
		int C=rs.getAtomContainerCount();
		
		//System.out.println("rs.size="+rs.size());
		
		
		if (C>0) {
		
			int n0=0;		
			int n1=0;
			
			for (int i = 0; i <= rs.getAtomContainerCount() - 1; i++) {
				Ring r = (Ring) rs.getAtomContainer(i);
				n1+=r.getAtomCount();
				
//				System.out.println(r.getAtomCount());
				
				n0+=A-r.getAtomCount();			
			}
			
			//System.out.println("n0="+n0);
			//System.out.println("n1="+n1);
			
			if (n0>0) {
				d.icyce=A*C*Log(2,A*C)-n1*Log(2,n1)-n0*Log(2,n0);
			} else {
				d.icyce=A*C*Log(2,A*C)-n1*Log(2,n1);
			}
			
			d.icycem=d.icyce/(A*C);
			d.icycm=n1*Log(2,n1);
			d.icycmm=d.icycm/n1;
			
		} else {
			d.icyce=0;
			d.icycem=0;
			d.icycm=0;
			d.icycmm=0;
		}
		
		
		
	}
	
	
	
	
	
	
	private void Calculate_iadje() {
		// Total Information Content On The Adjacency Equality
		d.iadje=0;
		
		double A2=A*A;
		
		d.iadje=A2*Log(2,A2)-2*B*Log(2,2*B)-(A2-2*B)*Log(2,A2-2*B);
						
		
	}
	
	private void Calculate_iadjm() {
		//Total Information Content On The Adjacency Magnitude
		d.iadjm=2*B*Log(2,2*B);
	}
	
	private void Calculate_iadjmm() {
		//Mean Information Content on the Adjacency Magnitude 
		// for CC(C)(C)C(C)CC, value = 3.8074 
		// which matches value on pg 104 of Bonchev, 1983
	
		
		
		d.iadjmm=1.0+Log(2,B);
	}
	
	private void Calculate_iadjem() {
		// Mean Information Content On The Adjacency Equality
		// for CC(C)(C)C(C)CC, value = 0.7579 
		// which matches value on pg 104 of Bonchev, 1983
		
		d.iadjem=0;
		
		double A2=A*A;
		
		double b=2.0*B/A2;
				
		d.iadjem=-b*Log(2,b)-(1-b)*Log(2,1-b);
		
		//System.out.println(d.iadjem);
		//System.out.println(d.iadje/A2); // alternate way of calculating iadjem
		
	}
	
	

	
	
	 private void Calculate_ieadje() {
		 //Total Information Content on the Edge Adjacency Equality
		 
		 double N2=paths[2].size();
		 
		 double B2=B*B;
		 if (N2>0) {
			 //Error fixed in formula 10/24/08
			 d.ieadje=B2*Log(2,B2)-2*N2*Log(2,2*N2)-(B2-2*N2)*Log(2,B2-2*N2);
		 } else {
			 d.ieadje=0;
		 }
		 
	 }
	 
	 private void Calculate_ieadjem() {
		// Total Information Content on the Edge Adjacency Equality
		// for CC(C)(C)C(C)CC, value = 0.9755
		// which matches value on pg 110 of Bonchev, 1983

		d.ieadjem = 0;
		double N2 = paths[2].size();

		double B2 = B * B;
		if (N2>0) {
		d.ieadjem = -2.0 * N2 / B2 * Log(2, 2.0 * N2 / B2)
				- (1 - 2.0 * N2 / B2) * Log(2, 1 - 2.0 * N2 / B2);
		} else {
			d.ieadjem=0;
		}

	}
	 
	 
	 private void Calculate_ieadjm() {
		// Total Information Content on the Edge Adjacency Magnitude (ieadjmm)
		
		 d.ieadjm=0;
		 
		 double N2=paths[2].size();
		 double Ae=2*N2;
		 if (N2>0) {
			 d.ieadjm=Ae*Log(2,Ae);
		 } else {
			 d.ieadjm=0;	 
		 }
		 
		 
	 }
	 
	 private void Calculate_ieadjmm() {		 
		 //Mean Information Content on the Edge Adjacency Magnitude (ieadjmm)

		 // for CC(C)(C)C(C)CC, value = 4.3219
		 // which matches value on pg 110 of Bonchev, 1983

		 
		 d.ieadjmm=0;
		 
		 double N2=paths[2].size();
		 if (N2>0) {
			 d.ieadjmm=1.0+Log(2,N2);
		 } else {
			 d.ieadjmm=0;	 
		 }
		 
		 
	 }
	 
	 
	 private void Calculate_idm() {
		 //Total Information Content on the Distance Magnitude (idm)
		 
		 // note: confirmed by Dragon 5.0 (idm=IDMT in Dragon)
		 // for test molecule 119-64-2
		 
		 d.idm=0;
		 
		 double sum=0;
		 
		 for (int g=1;g<=G;g++) {
			 sum+=fg[g]*g*Log(2,g);
		 }
		 d.idm=d.W*Log(2,d.W)-sum;
		 
		 
	 }
	 
	 private void Calculate_idmm() {
		 //Mean Information Content on the Distance Magnitude  (idmm)
		 // note: confirmed by Dragon 5.0 (idmm=IDM in Dragon)
		 // for test molecule 119-64-2
		 
		// for CC(C)(C)C(C)CC, value = 4.6751
		// which matches value on pg 130 of Bonchev, 1983
		 
		 
		 d.idmm=0;
		 
		 double sum=0;
		 
		 for (int g=1;g<=G;g++) {
			 sum+=fg[g]*g/d.W*Log(2,g/d.W);
		 }
		 
		 d.idmm=-sum;
		 
	 }
	 
	 
	 
	 private void Calculate_ide() {
		 //Total Information Content on the Distance Equality (ide)
//		 note: confirmed by Dragon 5.0 (ide=IDET in Dragon)
		 // for test molecule 119-64-2
		 
		 
		 d.ide=0;
		 
		 double num=A*(A-1)/2.0;
		 
		 double sum=0;
		 for (int g=1;g<=G;g++) {
			 sum+=fg[g]*Log(2,fg[g]);
		 }
		 
		 d.ide=num*Log(2,num)-sum;
		 
	 }

	 
	 private void Calculate_idem() {
		 //Mean Information Content on the Distance Equality (idem)
		 // note: confirmed by Dragon 5.0 (idem=IDE in dragon) for test molecule 119-64-2
		 
		// for CC(C)(C)C(C)CC, value = 1.8922
		// which matches value on pg 130 of Bonchev, 1983

		 d.idem=0;
		 
		 double sum=0;
		 
		 double denom=A*(A-1)/2.0;
		 		 
		 for (int g=1;g<=G;g++) {
			 sum+=fg[g]/denom*Log(2,fg[g]/denom);
		 }
		 d.idem=-sum;	
		 
		 //System.out.println(d.idem);
		 //System.out.println(d.ide*2/(A*(A-1)));
		 
		 
	 }
	 
	 
	 private void Calculate_iede() {
		 //Total Information Content on the Edge Distance Equality (iede)

		 
		 d.iede=0;
		 
		 double num=B*(B-1)/2.0;
		 
		 double sum=0;
		 
		 for (int g=1;g<=Gedge;g++) {
			 sum+=fgedge[g]*Log(2,fgedge[g]);
		 }
		 if (num>0) {
		 d.iede=num*Log(2,num)-sum;
		 } else {
			 d.iede=0;	 
		 }
		 
	 }
	 
	 private void Calculate_iedem() {
		 //Mean Information Content on the Edge Distance Equality (iedem)
		 
		 d.iedem=0;
		 //System.out.println(B);
		 
		 
		 double num=B*(B-1)/2.0;
		 
		 double sum=0;
		 for (int g=1;g<=Gedge;g++) {
			 sum+=fgedge[g]/num*Log(2,fgedge[g]/num);
		 }
		 
		 d.iedem=-sum;
		 
		 //System.out.println("iedem="+d.iedem);
		 //System.out.println("iedem2="+d.iede/num);
	 
		 
		 
		 
		 
	 }
	 
	 private void Calculate_iedm() {
		 //Total Information Content on the Edge Distance Magnitude (iedm)
		 
		 d.iedm=0;
		 
		 double sum=0;
		 for (int g=1;g<=Gedge;g++) {
			 sum+=fgedge[g]*g*Log(2,g);
		 }
		 if (d.We>0) {
			 d.iedm=d.We*Log(2,d.We)-sum;
		 } else {
			 d.iedm=0;	 
		 }
		 
		 
	 }
	 
	 private void Calculate_iedmm() {
		 //Mean Information Content on the Edge Distance Magnitude (iedmm)
		 d.iedmm=0;
		 
		 double sum=0;
		 for (int g=1;g<=Gedge;g++) {
			 sum+=fgedge[g]*g/d.We*Log(2,g/d.We);
			 			 
		 }
		 d.iedmm=-sum;
		 
		 //System.out.println("iedmm="+d.iedmm);
		 //System.out.println("iedmm2="+d.iedm/d.We);
	 	 		 
		 //reference 
		 
		 
	 }
	
	 
	 private int [][][] convertPaths() {
		 
		 int maxj=0;
		 int maxk=0;
		 
		 for (int i=1;i<=paths.length-1;i++) {
			 LinkedList al=(LinkedList)paths[i];
			 										 
			 ListIterator li=al.listIterator();
			
			 int j=0;
			 
			 while (li.hasNext()) {

				 LinkedList <Integer>al2=(LinkedList)li.next();
				 
				 ListIterator li2=al2.listIterator();
				 
				 int k=0;
				 while (li2.hasNext()) {
//					 System.out.println(i+"\t"+j+"\t"+k);
					 
					 int num=(Integer)(li2.next());
					 
					 k++;
					 if (k>maxk) maxk=k;
					 
				 }
				 j++;
				 
				 if (j>maxj) maxj=j;
				 
			 }
		 }
		 
		 int [][][] bob=new int[paths.length][maxj][maxk];
		 
		 
		 for (int i=1;i<=paths.length-1;i++) {
			 LinkedList al=(LinkedList)paths[i];
			 										 
			 ListIterator li=al.listIterator();
			
			 int j=0;
			 
			 while (li.hasNext()) {

				 LinkedList <Integer>al2=(LinkedList)li.next();
				 
				 ListIterator li2=al2.listIterator();
				 
				 int k=0;
				 while (li2.hasNext()) {
//					 System.out.println(i+"\t"+j+"\t"+k);
					 int num=(Integer)(li2.next());
					 
					 bob[i][j][k]=num;
					 k++;
				 }
				 j++;
			 }
		 }
		 
//		 System.out.println("maxj="+maxj);
//		 System.out.println("maxk="+maxk);
//		 System.out.println("atomCount="+m.getAtomCount());
		 
		 return bob;
	 }
	 
//	 private double calctti2(int b,int c,double [] Delta,double []S) {
//		 
//		 int maxj=0;
//		 int maxk=0;
//		 
//		 for (int i=1;i<=paths.length-1;i++) {
//			 LinkedList al=(LinkedList)paths[i];
//			 										 
//			 ListIterator li=al.listIterator();
//			
//			 int j=0;
//			 
//			 while (li.hasNext()) {
//
//				 LinkedList <Integer>al2=(LinkedList)li.next();
//				 
//				 ListIterator li2=al2.listIterator();
//				 
//				 int k=0;
//				 while (li2.hasNext()) {
////					 System.out.println(i+"\t"+j+"\t"+k);
//					 
//					 int num=(Integer)(li2.next());
//					 
//					 k++;
//					 if (k>maxk) maxk=k;
//					 
//				 }
//				 j++;
//				 
//				 if (j>maxj) maxj=j;
//				 
//			 }
//		 }
//		 
//		 int [][][] paths=new int[paths.length][maxj][maxk];
//		 
//		 
//		 for (int i=1;i<=paths.length-1;i++) {
//			 LinkedList al=(LinkedList)paths[i];
//			 										 
//			 ListIterator li=al.listIterator();
//			
//			 int j=0;
//			 
//			 while (li.hasNext()) {
//
//				 LinkedList <Integer>al2=(LinkedList)li.next();
//				 
//				 ListIterator li2=al2.listIterator();
//				 
//				 int k=0;
//				 while (li2.hasNext()) {
////					 System.out.println(i+"\t"+j+"\t"+k);
//					 int num=(Integer)(li2.next());
//					 
//					 paths[i][j][k]=num;
//					 k++;
//				 }
//				 j++;
//			 }
//		 }
//
//		 
//		 double [][] t=new double [(int)A][(int)A];
//		 
//		 try {
//		 
//		int totalnumpaths=0;	 
//		
//			 
//		 for (int i=0;i<=A-1;i++) {
//			 
//			 for (int j=0;j<=i;j++) {
//				 t[i][j]=0;
//				 
//				 
//				 for (int ii=0;ii<paths.length;ii++) {
//					 for (int jj=0;jj<maxj;jj++) {
//						 for (int kk=0;kk<maxk;kk++) {
//							 
//							 if ((first==i && last==j) || (first==j && last==i) ) {
//								 
//								 totalnumpaths++;
//								 							 
//								 double product=1;
//								 
//								 int nij=al2.size();
//								 
//								 ListIterator li2=al2.listIterator();
//								 
//								 while (li2.hasNext()) {							 
////									 String strnum=(String)li2.next();
////									 int num=Integer.parseInt(strnum);
//									 int num=(Integer)(li2.next());								 
//									 product*=Delta[num];								 								 
//								 }
//								 
//								 double GM=Math.pow(product,1.0/(double)nij);
//								 double fn=1.0/(double)nij;
//								 
//								 														
//								 t[i][j]+=Math.pow(GM,b)*Math.pow(fn,c);
//								 
//								 
//							 } // end if 
//						 }
//					 }
//				 }
//				 
//				 if (i!=j) t[j][i]=t[i][j];
//				 if (i==j) t[i][j]=Delta[i];
//				 
//				 
//				 
//			 }// end j loop
//			 
//		 } // end i loop
//		 
//		 
//		 
//		 for (int i=0;i<=m.getAtomCount()-1;i++) {
//			 S[i]=0;
//			 for (int j=0;j<=m.getAtomCount()-1;j++) {			 
//				 S[i]+=t[i][j];			 
//			 }
//			 //System.out.println("S["+i+"]\t"+S[i]+"\t"+Delta[i]);
//		 }
//		 
//		 // calculate sum for diagonal:
//		 double sumtii=0;
//		 for (int i=0;i<=A-1;i++) {			 
//			 sumtii+=t[i][i];
//		 }
//		 
//		 // calculate sum for upper triangular portion of matrix:
//		 double sumtij=0;
//		 for (int i=0;i<=A-2;i++) {			
//			 for (int j=i+1;j<=A-1;j++) {
//				 sumtij+=t[i][j];
//			 }
//		 }
//			 
//		 double tti=sumtii+sumtij;  // this calculation could be done in one loop if desired
//		 
//		 //System.out.println("tti="+tti);
//		 
//		 return tti;
//		 		 
//		 } catch (Exception e) {
//			 e.printStackTrace();
//			 return -9999;
//		 }
//		 
//	 }
	 
	 @Deprecated // use Calculate_tti2 instead
	 private double Calculate_tti(int b,int c,double [] Delta,double []S) {
		 
		 
		 
		 //also calculates S array (S or SV depending on what is sent to it)
		 
		 double [][] t=new double [(int)A][(int)A];
		 
		 try {
		 
		int totalnumpaths=0;	 
		
			 
		 for (int i=0;i<=A-1;i++) {
			 
			 for (int j=0;j<=i;j++) {
				 t[i][j]=0;
				 
				 for (int k=1;k<=paths.length-1;k++) {
					 LinkedList al=(LinkedList)paths[k];
					 										 
					 ListIterator li=al.listIterator();
					 
					 while (li.hasNext()) {
					 //for (int l=0;l<=al.size()-1;l++) {
						 
						 LinkedList <Integer>al2=(LinkedList)li.next();
						 				
//						 String strfirst=(String)al2.get(0);
//						 int first=Integer.parseInt(strfirst);
//						 
//						 String strlast=(String)al2.get(al2.size()-1);
//						 int last=Integer.parseInt(strlast);
												 
						 int first=al2.get(0);						 
						 int last=al2.get(al2.size()-1);
						 
						 						 
						 if ((first==i && last==j) || (first==j && last==i) ) {
							 
							 totalnumpaths++;
							 
							 							 
							 double product=1;
							 
							 int nij=al2.size();
							 
							 ListIterator li2=al2.listIterator();
							 
							 while (li2.hasNext()) {							 
//								 String strnum=(String)li2.next();
//								 int num=Integer.parseInt(strnum);
								 int num=(Integer)(li2.next());								 
								 product*=Delta[num];								 								 
							 }
							 
							 double GM=Math.pow(product,1.0/(double)nij);
							 double fn=1.0/(double)nij;
							 
							 														
							 t[i][j]+=Math.pow(GM,b)*Math.pow(fn,c);
							 
							 
						 } // end if
						 
					 } // end loop over paths of length k				
				 } // end k path length loop
				 				 				 
				 if (i!=j) t[j][i]=t[i][j];
				 if (i==j) t[i][j]=Delta[i];
				 
				 
				 
			 }// end j loop
			 
		 } // end i loop
		 
		 
		 
		 for (int i=0;i<=m.getAtomCount()-1;i++) {
			 S[i]=0;
			 for (int j=0;j<=m.getAtomCount()-1;j++) {			 
				 S[i]+=t[i][j];			 
			 }
			 //System.out.println("S["+i+"]\t"+S[i]+"\t"+Delta[i]);
		 }
		 //System.out.println("");
		 //System.out.println("totalnumpaths="+totalnumpaths);
		 
		 /*
		 for (int i=0;i<=m.getAtomCount()-1;i++) {
			 
			 for (int j=0;j<=m.getAtomCount()-1;j++) {			 
				 System.out.print(t[i][j]+"\t");
			 }
			 System.out.println("");
		 }
		*/
		 
		 // calculate sum for diagonal:
		 double sumtii=0;
		 for (int i=0;i<=A-1;i++) {			 
			 sumtii+=t[i][i];
		 }
		 
		 // calculate sum for upper triangular portion of matrix:
		 double sumtij=0;
		 for (int i=0;i<=A-2;i++) {			
			 for (int j=i+1;j<=A-1;j++) {
				 sumtij+=t[i][j];
			 }
		 }
			 
		 double tti=sumtii+sumtij;  // this calculation could be done in one loop if desired
		 
		 //System.out.println("tti="+tti);
		 
		 return tti;
		 		 
		 } catch (Exception e) {
			 e.printStackTrace();
			 return -9999;
		 }
		 
		 
	 }
	
	 
	 
     
    private double Calculate_tti2(int b, int c, double[] Delta, double[] S) {
        // also calculates S array (S or SV depending on what is sent to it)
        final int ac = m.getAtomCount();
        double[][] t = new double[ac][ac];
        try {
            for (LinkedList<LinkedList<Integer>> path : paths) {
                if (path == null) continue;
                
                for (LinkedList<Integer> al2 : path) {
                    int first = Math.min(al2.get(0), al2.get(al2.size() - 1));
                    int last = Math.max(al2.get(0), al2.get(al2.size() - 1));
                    int nij = al2.size();
                    
                    if (first != last && first < ac && last < ac) {
                        double product = 1.0d;
                        for (Integer idx : al2) {
                            product *= Delta[idx];
                        }
                        
                        double GM = Math.pow(product, 1.0d / nij);
                        double fn = 1.0d / (double) nij;
                        t[first][last] += Math.pow(GM, b) * Math.pow(fn, c);
                    }
                }
            }

            double tti = 0d;
            for (int i = 0; i <= A - 1; i++) {
                for (int j = i; j <= A-1; j++) {
                    if (i != j)
                        t[j][i] = t[i][j];
                    if (i == j)
                        t[i][j] = Delta[i];
                } // end j loop
            } // end i loop

            for (int i = 0; i <= m.getAtomCount() - 1; i++) {
                S[i] = 0;
                for (int j = 0; j <= m.getAtomCount() - 1; j++) {
                    S[i] += t[i][j];
                }
            }

            // calculate sum for diagonal:
            double sumtii = 0;
            for (int i = 0; i <= A - 1; i++) {
                sumtii += t[i][i];
            }

            // calculate sum for upper triangular portion of matrix:
            double sumtij = 0;
            for (int i = 0; i <= A - 2; i++) {
                for (int j = i + 1; j <= A - 1; j++) {
                    sumtij += t[i][j];
                }
            }

            tti = sumtii + sumtij; // this calculation could be done in one loop
                                   // if desired
                
            return tti;

        } catch (Exception e) {
            e.printStackTrace();
            return -9999;
        }
    }
	 
	 private static double Log(int base,double x) {
		
		double Logbx=0;
		
		Logbx=Math.log10(x)/Math.log10((double)base);
		
		return Logbx;
		
		
	}
	
}
