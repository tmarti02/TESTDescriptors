package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.lang.reflect.Field;
import java.util.*;

import org.openscience.cdk.interfaces.*;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
//import ToxPredictor.Utilities.CDKUtilities;
import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.AtomicProperties;

//TODO: figure out why eigenvalues for EState weights doesnt match DRAGON

public class WHIMDescriptors {
	IAtomContainer m;
	
	DescriptorData dd;

	double [] EState;
	
	double [] wu;
	double [] wm;
	double [] wv;
	double [] we;
	double [] wp;
	double [] ws;
	//int [][] DistanceWithH;
	
	int HCount=0;
	
	private static double Log(int base,double x) {
		
		double Logbx=0;
		
		Logbx=Math.log10(x)/Math.log10((double)base);
		
		return Logbx;
		
		
	}	
	
	public  void Calculate(IAtomContainer m,DescriptorData dd,double [] EState) {
		
		this.dd=dd;
		this.EState=EState;
		this.m=m;
		
		HCount=0;
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			if (m.getAtom(i).getSymbol().equals("H")) {
				HCount++;
			}
		}
		
		
		this.GetWeights(m);
		
		Matrix M=this.CalculateCoordinateMatrix(m);
		//Matrix M2=this.CalculateCoordinateMatrix(m2);  // hydrogen suppressed coord matrix
		
		PCA(M,wu,"u");
		PCA(M,wm,"m");
		PCA(M,wv,"v");
		PCA(M,we,"e");
		PCA(M,wp,"p");
		PCA(M,ws,"s");
		
	}
	
	
	
	private void PCA(Matrix M,double [] w,String strw) {

		
		//System.out.println(strw+":");
		Field myField;
		
		
		double covxx=0,covxy=0,covxz=0;
		double covyx=0,covyy=0,covyz=0;
		double covzx=0,covzy=0,covzz=0;
				
		double sumw=0;
		for (int i=0;i<=M.getRowDimension()-1;i++) {
			
			covxx+=w[i]*M.get(i,0)*M.get(i,0);
			covyy+=w[i]*M.get(i,1)*M.get(i,1);
			covzz+=w[i]*M.get(i,2)*M.get(i,2);
			
			covxy+=w[i]*M.get(i,0)*M.get(i,1);
			covxz+=w[i]*M.get(i,0)*M.get(i,2);
			covyz+=w[i]*M.get(i,1)*M.get(i,2);
			
			sumw+=w[i];			
		}
			
		covxx/=sumw;
		covyy/=sumw;
		covzz/=sumw;
		
		covxy/=sumw;
		covxz/=sumw;
		covyz/=sumw;
				
		covyx=covxy;
		covzx=covxz;
		covzy=covyz;
		
		
		double[][] covmatrix = { { covxx, covxy, covxz },
				{ covyx, covyy, covyz }, { covzx, covzy, covzz } };

		
		Matrix C=new Matrix(covmatrix);
		
		//C.print(6,4);
		
		EigenvalueDecomposition eid=C.eig();
		
		this.SortEID(eid);
			

		double [] eigenvalues=eid.getRealEigenvalues();		
		Matrix eigenvectors=eid.getV();

		for (int i=0;i<=eigenvalues.length-1;i++) {
			if (eigenvalues[i]<1e-6) eigenvalues[i]=0;
			//System.out.println(eigenvalues[i]);
		}

		
		
		//eigenvectors.print(6,4);
		for (int i=0;i<=eigenvalues.length-1;i++) {
			if (eigenvalues[i]<1e-6) eigenvalues[i]=0;
		//	System.out.println(eigenvalues[i]);
		}
	
		//M.print(6,4);
		Matrix T=M.times(eigenvectors);
//		System.out.println(strw);
//		T.print(14,12);
		
		
		
		//Matrix T=eigenvectors.transpose().times(M.transpose()).transpose();
//		T.print(14,12);
		
		
		try {
			
			double sum=0;
			
			for (int i=1;i<=3;i++) {
				myField = dd.getClass().getField("L"+i+strw);
				myField.setDouble(dd, eigenvalues[i-1]); // eigenvalues are in ascending order- need in descending
				//System.out.println("L"+i+strw+"\t"+eigenvalues[i-1]);
				sum+=eigenvalues[3-i];
			}
			
			myField = dd.getClass().getField("T"+strw);
			myField.setDouble(dd, sum); // 
			
			double A=eigenvalues[0]*eigenvalues[1];
			A+=eigenvalues[1]*eigenvalues[2];
			A+=eigenvalues[0]*eigenvalues[2];
			
			myField = dd.getClass().getField("A"+strw);
			myField.setDouble(dd, A); // 
			
			double V=sum+A+eigenvalues[0]*eigenvalues[1]*eigenvalues[2];
			myField = dd.getClass().getField("V"+strw);
			myField.setDouble(dd, V); // 
			
			// ***************************************************
			// calculate K
			double K=0;
			
			for (int i=1;i<=3;i++) {
				K+=Math.abs(eigenvalues[3-i]/sum-1.0/3.0);
			}			
			K*=3.0/4.0;
			
			myField = dd.getClass().getField("K"+strw);
			myField.setDouble(dd, K); 
						
			// ***************************************************
						
			for (int i=1;i<=2;i++) {
				myField = dd.getClass().getField("P"+i+strw);
				myField.setDouble(dd, eigenvalues[i-1]/sum); // eigenvalues are in ascending order- need in descending
				//System.out.println("P"+i+strw+"\t"+eigenvalues[3-i]/sum);			
			}
			
			double D=0;
			
			for (int mm=0;mm<=2;mm++) {
				
				//http://en.wikipedia.org/wiki/Moment_%28mathematics%29
					
				double sum3=0;
				
				for (int i=0;i<=M.getRowDimension()-1;i++) {
					sum3+=Math.pow(T.get(i,mm),4);
				}
				
				//System.out.println(sum2);
				double etam;
				if (sum3==0) {
					etam=0;
				} else {
					etam=(Math.pow(eigenvalues[mm],2)*M.getRowDimension())/sum3;
				}
				
				D+=etam;
								
				myField = dd.getClass().getField("E"+(mm+1)+strw);
				myField.setDouble(dd, etam); // eigenvalues are in ascending order- need in descending
										
			}

			D/=3.0;
						
			myField = dd.getClass().getField("D"+strw);
			myField.setDouble(dd, D); // eigenvalues are in ascending order- need in descending

			
			
			//System.out.println("");
			
			
			this.Symmetry(T,strw);
			//this.Symmetry(M,eigenvectors);
			
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}
	
	private void SortEID(EigenvalueDecomposition eid) {
		
		double [] eigenvalues=(double [])eid.getRealEigenvalues().clone();
		Matrix ev=(Matrix)eid.getV().clone();
	
		double max=-1,min=99999999,middle=0;
		int imax=-1,imin=-1,imid=-1;
		
		for (int i=0;i<=2;i++) {
			if (eigenvalues[i]>max) {
				max=eigenvalues[i];
				imax=i;				
			}
			if (eigenvalues[i]<min) {
				min=eigenvalues[i];
				imin=i;
			}
		}
		for (int i=0;i<=2;i++) {
			if (i!=imin && i!=imax) {
				imid=i;
			}
		}
		
		//System.out.println(imax+"\t"+imid+"\t"+imin);
		
		eid.getRealEigenvalues()[0]=eigenvalues[imax];
		eid.getRealEigenvalues()[1]=eigenvalues[imid];
		eid.getRealEigenvalues()[2]=eigenvalues[imin];
		
		for (int i=0;i<=2;i++) {
			eid.getV().set(i,0,ev.get(i,imax));
			eid.getV().set(i,1,ev.get(i,imid));
			eid.getV().set(i,2,ev.get(i,imin));
			//eid.get
		}
		
		
		
		
	}
	
	private void Symmetry(Matrix T,String strw) {
		// for atoms to be symmetric wrt to an axis the atoms must have 
		// opposite coordinates along the considered axis or must lie on the axis
		double xi,yi,zi,xj,yj,zj;
		
		double tol=0.005;//FIXME maybe can reduce tol if our coordinates are more accurate
		
		double coordi,coordj;
		
		double [] ns=new double [3];
		Field myField;
		
		ArrayList [] al=new ArrayList[3];
		
		for (int i=0;i<=2;i++) {
			al[i]=new ArrayList();
			ns[i]=0;
		}
		
//		T.print(6,4);
		
		for (int i=0;i<=T.getRowDimension()-1;i++) {
			
			if (strw.equals("s")) {
				if (this.m.getAtom(i).getSymbol().equals("H")) {
					continue;
				}
			}
			
			for (int k=0;k<=2;k++) {
			
				coordi=T.get(i,k);
				
				if (Math.abs(coordi)<tol) {
					al[k].add(new Integer(i));
				}
			
				for (int j=i+1;j<=T.getRowDimension()-1;j++) {
				
					if (strw.equals("s")) {
						if (this.m.getAtom(j).getSymbol().equals("H")) {
							continue;
						}
					}
					
					
					coordj=T.get(j,k);
					
					
					
					if (Math.abs(coordi+coordj)<tol) {
						al[k].add(new Integer(i));
						al[k].add(new Integer(j));
					}
																	
				}//end j for loop
			} // end k loop
			
		}// end i for loop
	
		// remove duplicates:
		for (int k=0;k<=2;k++) {
			Collections.sort(al[k]);
			
			for (int i=0;i<=al[k].size()-2;i++) {
				if (al[k].get(i).equals(al[k].get(i+1))) {
					al[k].remove(i);
					i=-1;
				}
			}			
			ns[k]=al[k].size();
		}
		
		

		
		double n;
		
		if (strw.equals("s")) {
			n=(double)T.getRowDimension()-HCount;
		} else {
			n=(double)T.getRowDimension();
		}
		
		
		double G=1;
		for (int k=0;k<=2;k++) {
						
			double na=n-ns[k];
			
//			System.out.println(strw+"\t"+k);
//			System.out.println("ns"+k+"\t"+ns[k]);
//			System.out.println("na"+k+"\t"+na);
//			System.out.println("n"+k+"\t"+n+"\n");
			
			double gammap;
			if (ns[k]>0) {
				gammap=-(ns[k]/n*Log(2,ns[k]/n)+na*(1/n*Log(2,1/n)));
			} else {
				gammap=-(Log(2,1/n));
			}
			double gamma=1/(1+gammap);
			
			G*=gamma;
									
			//System.out.println(k+"\t"+gamma);
			try {
				myField = dd.getClass().getField("G"+(k+1)+strw);
				myField.setDouble(dd, gamma); // eigenvalues are in ascending order- need in descending
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			
		}// end k for loop
		
		G=Math.pow(G,1.0/3.0);
		
		try {
			if (strw.equals("u") || strw.equals("m") || strw.equals("s")) {
				myField = dd.getClass().getField("G"+strw);
				myField.setDouble(dd, G); // 
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
			
	}
	
	
	
	

	
	
	private Matrix CalculateCoordinateMatrix(IAtomContainer m) {
	    double [][] Coords=new double [m.getAtomCount()][3];
	    
	    double Xavg=0,Yavg=0,Zavg=0;
	    
	    // first make coords relative to the geometric center:
	    
	    for (int i=0;i<=m.getAtomCount()-1;i++){			
	    	IAtom a=m.getAtom(i);	    	
	    	Xavg+=a.getPoint3d().x;
	    	//System.out.println(a.getX3d());
	    	Yavg+=a.getPoint3d().y;
	    	Zavg+=a.getPoint3d().z;		
	    }
		
	    Xavg/=m.getAtomCount();
	    Yavg/=m.getAtomCount();
	    Zavg/=m.getAtomCount();
	    
	   // System.out.println(Xavg+"\t"+Yavg+"\t"+Zavg);
	    
	    for (int i=0;i<=m.getAtomCount()-1;i++){
			IAtom a=m.getAtom(i);			
			
	    	Coords[i][0]=a.getPoint3d().x-Xavg;
			Coords[i][1]=a.getPoint3d().y-Yavg;
			Coords[i][2]=a.getPoint3d().z-Zavg;		

			
	    }
	    Matrix M=new Matrix(Coords);
	    return M;
	    //M.print(6,4);
	}
	
	
	
	private void GetWeights(IAtomContainer m) {
		
		wu=new double [m.getAtomCount()] ;
		wm=new double [m.getAtomCount()] ;
		wv=new double [m.getAtomCount()] ;
		we=new double [m.getAtomCount()] ;
		wp=new double [m.getAtomCount()] ;
		ws=new double [m.getAtomCount()] ;
		
		AtomicProperties ap;
		
		try {
			ap=AtomicProperties.getInstance();
			
			for (int i=0;i<=m.getAtomCount()-1;i++) {
				String s=m.getAtom(i).getSymbol();
				wu[i]=1;
				wm[i]=ap.GetNormalizedMass(s);
				wv[i]=ap.GetNormalizedVdWVolume(s);
				we[i]=ap.GetNormalizedElectronegativity(s);
				wp[i]=ap.GetNormalizedPolarizability(s);
				if (i<EState.length) ws[i]=EState[i]+7;
				else ws[i]=0;
								
			}
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		}
	
}

