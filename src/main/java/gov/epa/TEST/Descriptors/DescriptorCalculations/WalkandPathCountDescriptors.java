package gov.epa.TEST.Descriptors.DescriptorCalculations;
import java.lang.reflect.Field;
import java.util.LinkedList;
import java.util.Iterator;
import java.util.ListIterator;

import Jama.*;
//import ToxPredictor.Utilities.Matrix2;
import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

public class WalkandPathCountDescriptors {

	/**
	 * @param args
	 */
	
	
	IAtomContainer m;
	DescriptorData dd;
	LinkedList[] paths;
	double [] D;
	int []vdd;
	IRingSet rs;
	
	public void Calculate(IAtomContainer m,DescriptorData dd,LinkedList[] paths,double [] D,int [] vdd,IRingSet rs) {
		this.m=m;
		this.dd=dd;
		this.paths=paths;
		this.D=D;
		this.vdd=vdd;
		this.rs=rs;
		
		this.CalculateMolecularWalkCounts();
		this.CalculateSelfReturningWalkCounts();
		this.CalculatePathCounts();
		this.CalculateCID();
		this.CalculateBID();
		this.CalculateMultiplePathCountDescriptors();
				
		//benzimidazlone: c1ccc2ncnc2(c1)
	}
	
	private void CalculateCID() {
				
		double sum=0;
		for (int k=1;k<=paths.length-1;k++){
			LinkedList ll=paths[k];
			
			for (int i=0;i<=ll.size()-1;i++) {
				LinkedList <Integer>ll2=(LinkedList)ll.get(i);
		
				double wij=1;
				for (int j=0;j<=ll2.size()-2;j++) {
//					String strOne=(String)ll2.get(j);
//					String strTwo=(String)ll2.get(j+1);
//					int One=Integer.parseInt(strOne);
//					int Two=Integer.parseInt(strTwo);
					
					int One=ll2.get(j);
					int Two=ll2.get(j+1);
					
					double d1=D[One];
					double d2=D[Two];
					
					//System.out.println("k="+k+"\ti="+i+"\t"+d1+"\t"+d2);
					wij*=Math.pow(d1*d2,-0.5);
				}
				sum+=wij;
			}
		}
		dd.CID=m.getAtomCount()+sum;
		dd.CID2=dd.CID/m.getAtomCount();
	}
	
	private void CalculateBID() {
		
		double sum=0;
		for (int k=1;k<=paths.length-1;k++){
			LinkedList ll=paths[k];
			
			for (int i=0;i<=ll.size()-1;i++) {
				LinkedList <Integer>ll2=(LinkedList)ll.get(i);
		
				double wij=1;
				for (int j=0;j<=ll2.size()-2;j++) {
//					String strOne=(String)ll2.get(j);
//					String strTwo=(String)ll2.get(j+1);
//					int One=Integer.parseInt(strOne);
//					int Two=Integer.parseInt(strTwo);
					
					int One=(Integer)ll2.get(j);
					int Two=(Integer)ll2.get(j+1);

					double d1=vdd[One];
					double d2=vdd[Two];
					
					//System.out.println("k="+k+"\ti="+i+"\t"+d1+"\t"+d2);
					wij*=Math.pow(d1*d2,-0.5);
				}
				sum+=wij;
			}
		}
		dd.BID=m.getAtomCount()+sum;
		
	}
	
	
	private void CalculatePathCounts() {
		String s;
		dd.TPC=0;
		
		double TPCnew=0;
		
		try {
				
			for (int k=1;k<=paths.length-1;k++) {
				
				
				if (k<=9) s="MPC0"+k;
				else s="MPC10";
				
				Field myField = dd.getClass().getField(s);
				
				
				if (k<=10) {
					myField.setDouble(dd, paths[k].size());// init to zero				
				}
				
				TPCnew+=paths[k].size();
				
				dd.TPC+=Math.log(1+paths[k].size());
//				dd.TPC+=Math.log(1+paths[k].size()+dd.nAT); // to match dragon
				
			}
			dd.TPC+=m.getAtomCount();
//			dd.TPC+=dd.nAT; // to match dragon
			TPCnew+=m.getAtomCount();
			TPCnew=Math.log(1+TPCnew);
			
//			System.out.println("TPCnew="+TPCnew);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private void CalculateMolecularWalkCounts() {

		try {
		// calculate adjacency matrix:
		double [][] am = this.getAdjacencyMatrix(m);
		
		Matrix matrix=new Matrix(am);
		Matrix newmatrix=new Matrix(am);
		
		String s=null;
		
		dd.TWC=0;
		
		double TWCnew=0;////in version 6.0 of dragon MWC1 is calc using ln(1+x) transformation at the very end
		double MWC1new=0;//in version 6.0 of dragon MWC1 is calc using ln(1+x) transformation
		
		
		Field myField=null;
		
		int maxk=0;
		
		if (m.getAtomCount()<10) {
			maxk=10;
		} else {
			maxk=m.getAtomCount();
		}
		
		for (int k=1;k<=maxk;k++) {
			
			if (k<=10) {
				if (k<=9) s="MWC0"+k;
				else s="MWC10";					
				myField = dd.getClass().getField(s);	
				myField.setDouble(dd,0.0);
			}
		
			if (k>1) {
				newmatrix=newmatrix.times(matrix);
			}
											

			double MWCk; 
			double MWCk2;
			
			if (k==1) {
				MWCk=m.getBondCount();
				MWCk2=MWCk;
				
				MWC1new=Math.log(MWCk+1);
				
				TWCnew+=MWCk;
			} else {
				MWCk=this.sumMatrix(newmatrix);
				// why dont divide by 2?- see rucker and rucker, 2000
				// doesnt agree with pg 481 of Todeschini book
				
				MWCk2=MWCk;
				
				if (k<=10) {
					TWCnew+=MWCk;
				}
				
				MWCk=Math.log(MWCk+1); // do log transform to match dragon
			}
			
			
			//System.out.println(k+"\t"+MWCk2);
			
			if (k<=10) {
				myField.setDouble(dd, MWCk);
			}
			
//			TODO: Apparently Dragon adds 'em up all the way to 10 instead 
			// of stopping at A-1.  Mistake??? - see pg 482 todeschini book
//			if (k<=m.getAtomCount()-1) {
			if (k<=10) { // to match dragon 
				dd.TWC+=MWCk; 
			}
			
			//System.out.println(s+"\t"+myField.getDouble(dd));
			
		}// end k for loop

		dd.TWC+=m.getAtomCount();
		
		TWCnew+=m.getAtomCount();
		TWCnew=Math.log(TWCnew+1);
//		System.out.println(TWCnew);
		
//		System.out.println("MWC1new="+MWC1new);
		
		//System.out.println("TWC="+dd.TWC);
		
		
					
		
		
		} catch (Exception e) {
			e.printStackTrace();
		}
	
	}
	
	

private void CalculateSelfReturningWalkCounts() {
		
		try {
		// calculate adjacency matrix:
		double [][] am = this.getAdjacencyMatrix(m);
		
		Matrix matrix=new Matrix(am);
		
//		matrix.print(6, 4);
		
		Matrix newmatrix=new Matrix(am);
		
		String s=null;
		
		
		for (int k=1;k<=10;k++) {
			if (k<=9) s="SRW0"+k;
			else s="SRW10";					
			
			Field myField = dd.getClass().getField(s);
			myField.setDouble(dd, 0); //init to zero
			
			if (k==1) {
				myField.setDouble(dd, m.getAtomCount());
			} else {
				newmatrix=newmatrix.times(matrix);
				myField.setDouble(dd, newmatrix.trace());				
			}
							
			//System.out.println(s+"\t"+myField.getDouble(dd));
											
		}
				
		
		} catch (Exception e) {
			e.printStackTrace();
		}
	
	}

	
	public static double[][] getAdjacencyMatrix(IAtomContainer container) {
		//made a new version of this so can have double [][]
		int [][] matrix=org.openscience.cdk.graph.matrix.AdjacencyMatrix.getMatrix(container);
		double[][] conMat = new double[container.getAtomCount()][container.getAtomCount()];
		for (int i=0;i<container.getAtomCount();i++) {
			for (int j=0;j<container.getAtomCount();j++) {
				conMat[i][j]=(double)matrix[i][j];
			}
		}
		return conMat;
	}
		
	private void CalculateMultiplePathCountDescriptors() {
		//benzimidazlone: c1ccc2ncnc2(c1)
		
		double piID=0;
		double piIDNew=0;//attempt to match piID in Dragon v. 6.0
		
		Field myField=null;
		
		for (int k=1;k<=paths.length-1;k++) {
			 
			if (k<=10) {
				String s;
				if (k<=9) s="piPC0"+k;
				else s="piPC10";					
				try {
					myField = dd.getClass().getField(s);
					myField.setDouble(dd,0.0);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		
			
			
			LinkedList al=(LinkedList)paths[k];
			 //System.out.println(k);
			 
			 ListIterator li=al.listIterator();
			 
			 double piPCk=0;
			 
			 while (li.hasNext()) {
			 //for (int l=0;l<=al.size()-1;l++) {
				 
				 LinkedList <Integer>al2=(LinkedList)li.next();
				 
				 double wij=1;
				 for (int i=0;i<=al2.size()-2;i++) {
//					 String strfirst=(String)al2.get(i);
//					 int first=Integer.parseInt(strfirst);
//					 
//					 String strsecond=(String)al2.get(i+1);
//					 int second=Integer.parseInt(strsecond);
					 
					 
					 int first=(al2.get(i));
					 int second=(al2.get(i+1));

					 IAtom a1=m.getAtom(first);
					 IAtom a2=m.getAtom(second);
					 
					 
					 
					 double bondorder;

					 // make sure aromatic bond orders are 1.5:
					 if (a1.getFlag(CDKConstants.ISAROMATIC)
								&& a2.getFlag(CDKConstants.ISAROMATIC)) {
						 
						 // check if in same ring:
						 boolean SameRing = EStateFragmentDescriptor.InSameAromaticRing(m,a1,
									a2, rs);
						 if (SameRing) {
							 bondorder=1.5;
						 } else {
							 IBond b=m.getBond(a1,a2);
							 bondorder=(double)b.getOrder().numeric();
						 }
						 						 						 
					 } else {
						 IBond b=m.getBond(a1,a2);
						 bondorder=(double)b.getOrder().numeric();
					
					 }
					
					 
					 //System.out.println(first+"\t"+second+"\t"+bondorder);
					 
//					 int d1=m.getBondCount(m.getAtom(i));
//					 int d2=m.getBondCount(m.getAtom(i+1));
//					 
//					 wij*=Math.pow((double)(d1*d2),-0.5);
					 
					 
					
					 
					 //wij*=Math.pow(bondorder,-1); // use fractional bond order
					 wij*=bondorder;					 
					 
					 //System.out.println(first+"\t"+second+"\t"+wij);
					 
				 }
				//System.out.print("\n");
				// System.out.println(wij);
				 
				 piPCk+=wij;
				 						 
				 
			 } // end loop over paths of length k

			 piIDNew+=piPCk;
			 
			 //System.out.println(k+"\t"+piIDk+"\t"+Math.log(1+piIDk));
			 piPCk=Math.log(1+piPCk);
			 
			 
			 if (k<=10) {
				 try {
					 myField.setDouble(dd,piPCk);
				 } catch (Exception e) {
					 e.printStackTrace();
				 }
			 }
			 
			 piID+=piPCk;
			 
			 
			 
		} // end k path length loop
		
		piIDNew+=m.getAtomCount();	
		piIDNew=Math.log(1+piIDNew);
		
//		System.out.println(piIDNew);
		
		piID+=m.getAtomCount();
		dd.piID=piID;
		
		//System.out.println("piID="+piID);
	}

	public double sumMatrix (Matrix matrix) {
		double sum = 0;
		
		int m = matrix.getRowDimension();
		int n = matrix.getColumnDimension();
		
		double [][]A=matrix.getArray();
		
		for (int i=0;i<=m-1;i++) {
			for (int j=0;j<=n-1;j++) {
				sum+=A[i][j];
			}
			
		}
		return sum;
	}
	
	
	
	
}
