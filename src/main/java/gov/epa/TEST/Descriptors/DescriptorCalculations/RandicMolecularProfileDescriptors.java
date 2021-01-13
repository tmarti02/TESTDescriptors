package gov.epa.TEST.Descriptors.DescriptorCalculations;
//TODO finish these descriptors!

import java.lang.reflect.Field;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import Jama.*;
//import ToxPredictor.Utilities.CDKUtilities;
import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

public class RandicMolecularProfileDescriptors {
	
	//Molecule m;
	DescriptorData dd;
	
	public void Calculate(IAtomContainer molecule,DescriptorData dd) {
		//this.m=m;
		this.dd=dd;
				
		IAtomContainer m2=null;
		
		try {
		 m2=(IAtomContainer)molecule.clone();
		
		 //		m2=3d hydrogen suppressed graph
		 m2=(IAtomContainer)AtomContainerManipulator.removeHydrogens(m2); 		 
//		 CDKUtilities.RemoveHydrogens(m2);
		 
		} catch( Exception e) {
			e.printStackTrace();
		}
		
		Matrix Gm=this.CalculateGeometryMatrix(m2);		
		//Gm.print(6,4);
		
//		for (int i=0;i<=m2.getAtomCount()-1;i++) {
//			System.out.println(m.getAtom(i).getSymbol());
//		}
		
		double fact=1;
		
		for (int k=1;k<=20;k++) {
			fact*=k;
			double DSum=0;
			double SSum=0;
			double SCount=0;
			
			for (int i=0;i<=Gm.getRowDimension()-1;i++) {
				
				IAtom atomi=m2.getAtom(i);
				
				if (m2.getConnectedBondsCount(atomi)<3) {
					SCount++;
				}
				
				for (int j=0;j<=Gm.getColumnDimension()-1;j++) {
					IAtom atomj=m2.getAtom(j);
					DSum+=Math.pow(Gm.get(i,j),k);
					
					if (m2.getConnectedBondsCount(atomi)<3 && m2.getConnectedBondsCount(atomj)<3) {
//						if (k==1) System.out.println(atomi.getSymbol()+"\t"+i+"\t"+j);
						SSum+=Math.pow(Gm.get(i,j),k);
					}
										
				}
			}
			 
			//if (k==1) System.out.println(m2.getAtomCount()+"\t"+SCount);
			double DPk=DSum/(m2.getAtomCount()*fact);
			DPk=Math.log(DPk+1);
			
			double SPk=SSum/(SCount*fact);
			SPk=Math.log(SPk+1);
			
			String strNum="";
			
			if (k<=9) strNum="0"+k;
			else strNum=k+"";
						
			try {
				Field myField = dd.getClass().getField("DP"+strNum);
				myField.setDouble(dd, DPk);
				
				myField = dd.getClass().getField("SP"+strNum);
				myField.setDouble(dd, SPk);				
			
			} catch (Exception e) {
				e.printStackTrace();
			}

			
			//System.out.println(k+"\t"+DPk+"\t"+SPk);
		}
		
		//The average shape profile index of order 2 (SHP2)
		//is the shape profile index of order 2 (SP02) divided 
		// by the number of atoms with H-depleted connectivity 
		// equal to 1 or 2. 

		int count=0;
		for (int i=0;i<=m2.getAtomCount()-1;i++) {
						
			if (m2.getConnectedBondsCount(m2.getAtom(i))<=2) {
				count++;
			}
		}
		
		dd.SHP2=dd.SP02/(double)count;
		
	}

	private static Matrix CalculateGeometryMatrix(IAtomContainer m) {
		double [][] G=new double [m.getAtomCount()][m.getAtomCount()];
		
		double Gmax=0;		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			IAtom ai=m.getAtom(i);
			
			G[i][i]=0;
			for (int j=i+1;j<=m.getAtomCount()-1;j++) {
				IAtom aj=m.getAtom(j);
				G[i][j]=Distance(ai,aj);
				//G[i][j]/=/1.39;
				
				
				if (G[i][j]>Gmax) Gmax=G[i][j];
				G[j][i]=G[i][j];
								
			}
		}
		
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			for (int j=0;j<=m.getAtomCount()-1;j++) {
//				G[i][j]/=Gmax;
//			}
//		}

		//System.out.println("Gmax="+Gmax);
		Matrix Gm=new Matrix(G);
		return Gm;
	}

	private static double  Distance(IAtom ai,IAtom aj) {
		double distance=0;
		
		distance+=Math.pow(ai.getPoint3d().x-aj.getPoint3d().x,2);
		distance+=Math.pow(ai.getPoint3d().y-aj.getPoint3d().y,2);
		distance+=Math.pow(ai.getPoint3d().z-aj.getPoint3d().z,2);
		distance=Math.sqrt(distance);		
		
		return distance;
		
		
		
		
	}

}




