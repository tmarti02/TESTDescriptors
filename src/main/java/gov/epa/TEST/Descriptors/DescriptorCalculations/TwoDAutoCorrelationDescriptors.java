package gov.epa.TEST.Descriptors.DescriptorCalculations;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.AtomicProperties;

//import ToxPredictor.Utilities.CDKUtilities;

import java.lang.reflect.Field;

public class TwoDAutoCorrelationDescriptors {

	IAtomContainer m;
	DescriptorData dd;
	int Distance [][];
	double [] wm;
	double [] wv;
	double [] we;
	double [] wp;

	
	
	public void Calculate(IAtomContainer m, DescriptorData dd, int [][]Distance) {
		this.m=m;
		
		this.dd=dd;
		this.Distance=Distance;
		
//		try {
//			m2=(Molecule)m.clone();
//			CDKUtilities.RemoveHydrogens(m2); //m2=3d hydrogen suppressed graph
//		} catch( Exception e) {
//			e.printStackTrace();
//		}
//		
//		
//		ToxPredictor.Utilities.CDKUtilities.AddHydrogens(m2);
//		Distance=new int[m2.getAtomCount()][m2.getAtomCount()];
//		PathFinder.CalculateDistanceMatrix(m2,Distance);
		
		this.GetWeights();
		
		this.Calculate2(wm,"m");
		this.Calculate2(wv,"v");
		this.Calculate2(we,"e");
		this.Calculate2(wp,"p");
	}
	
	
	/**
	 * Calculate the descriptors for a given set of weights
	 * @param w
	 * @param strw
	 */
	private void Calculate2(double []w,String strw) {

		// E-mail from Todeschini regarding ATS descriptors:
		//
		// The calculations are performed on the whole molecule (H included),
		// using the true weights, but the final value is calculated as ln(1+w),
		// avoiding too big numbers for big molecules
		// However, in the Dragon version 5.4 (that will be released in the next weeks)
		// we made some changes.
		// We came back to the H-depleted calculations, but using the scaled
		// weights (on the Carbon atom). This modification avoid the strong
		// correlation of ATS descriptors when too different weights are used
		// (see articles of Hollas).
		// It will remain the logarithm transformation.
		
		
		double Deltaij,SumDeltaij,SumGeary,wmean,Sumw,cd;
		double SumMoran,Id;

		double SumMB, ATSd;

		
		
		wmean=0;
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			wmean+=w[i];
		}
		wmean/=m.getAtomCount();
		
		//System.out.println("wmean="+wmean);
		
		Sumw=0;
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			Sumw+=Math.pow(w[i]-wmean,2);			
		}
		
		
		for (int d=1;d<=8;d++) {
			SumDeltaij=0;
			SumGeary=0;
			SumMoran=0;
			SumMB=0;
			
			for (int i=0;i<=m.getAtomCount()-1;i++) {								
				for (int j=0;j<=m.getAtomCount()-1;j++) {
					//System.out.println(i+"\t"+j+"\t"+Distance[i][j]);					
					if(Distance[i][j]==d) {
						Deltaij=1;
					} else {
						Deltaij=0;
					}
					
					if (j>i) 
						SumMB+=Deltaij*w[i]*w[j];
					
					SumMoran+=Deltaij*(w[i]-wmean)*(w[j]-wmean);
					SumGeary+=Deltaij*Math.pow(w[i]-w[j],2.0);										
					
					SumDeltaij+=Deltaij;
				}
			}
			
//			System.out.println(SumDeltaij);
//			System.out.println(Sum);
			
			if (SumDeltaij>0) {				
				ATSd=Math.log(1+SumMB);//to match Dragon
				//System.out.println(strw+d+"\t"+ATSd);
				
				if (Sumw>0) {
					Id=SumMoran*m.getAtomCount()/SumDeltaij/Sumw;
					cd=SumGeary*(m.getAtomCount()-1)/(2*SumDeltaij)/Sumw;
				} else {
					Id=1;
					cd=0;
				}
				
//				System.out.println(strw+"\t"+SumMoran+"\t"+SumDeltaij+"\t"+Sumw);
				
			} else {
				cd=0;				
				Id=0;
				ATSd=0;
			}
									
			try {
			
				Field myField;
				
				myField = dd.getClass().getField("ATS"+d+strw);
				myField.setDouble(dd,ATSd);
								
				myField = dd.getClass().getField("GATS"+d+strw);
				myField.setDouble(dd,cd);
				
				myField = dd.getClass().getField("MATS"+d+strw);
				myField.setDouble(dd,Id);
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			//System.out.println("c"+d+"\t"+cd);
		} // end d for loop
		
	}

	private void GetWeights() {
		
		
		wm=new double [m.getAtomCount()] ;
		wv=new double [m.getAtomCount()] ;
		we=new double [m.getAtomCount()] ;
		wp=new double [m.getAtomCount()] ;
		
		
		AtomicProperties ap;
		
		try {
			ap=AtomicProperties.getInstance();
			
			for (int i=0;i<=m.getAtomCount()-1;i++) {
				String s=m.getAtom(i).getSymbol();
				
				wm[i]=ap.GetNormalizedMass(s);
				wv[i]=ap.GetNormalizedVdWVolume(s);
				we[i]=ap.GetNormalizedElectronegativity(s);
				wp[i]=ap.GetNormalizedPolarizability(s);
//				System.out.println(s+"\t"+wm[i]+"\t"+wv[i]+"\t"+we[i]+"\t"+wp[i]);
			}
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
}
