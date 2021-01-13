package gov.epa.TEST.Descriptors.DescriptorCalculations;

import org.openscience.cdk.interfaces.*;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.AtomicProperties;

public class PrincipalMomentOfInertiaDescriptors {
	
	double CMX,CMY,CMZ,SUMWT;
	IAtomContainer m;
	DescriptorData dd;
	
	double  [] wm;
	
	public  void Calculate(IAtomContainer m,DescriptorData dd) {
		this.m=m;
		this.dd=dd;
		
		Matrix M=this.CalculateCoordinateMatrix(m);
		this.GetWeights(m);
		this.CalculateCenterOfMass(M,wm);
		this.CalcMoments(M,wm);
		
	}
	
	
	
	private void CalcMoments(Matrix M,double [] wm) {
		// *** CALCULATE THE MOMENTS AND PRODUCTS OF INERTIA ABOUT THE C OF MASS
	      double RGX=0.0;
	      double RGY=0.0;
	      double RGZ=0.0;
	      double RGXY=0.0;
	      double RGXZ=0.0;
	      double RGYZ=0.0;
	      
	      for (int i=0;i<=m.getAtomCount()-1;i++) {	
	         double XDIF = M.get(i,0) - CMX;
	         double YDIF = M.get(i,1) - CMY;
	         double ZDIF = M.get(i,2) - CMZ;
	         	         
	         double XSQ = XDIF*XDIF;
	         double YSQ = YDIF*YDIF;
	         double ZSQ = ZDIF*ZDIF;
	         
	         RGX += wm[i]*(YSQ + ZSQ);
	         RGY += wm[i]*(XSQ + ZSQ);
	         RGZ += wm[i]*(XSQ + YSQ);
	         
	         RGXY += -wm[i]*XDIF*YDIF;
	         RGXZ += -wm[i]*XDIF*ZDIF;
	         RGYZ += -wm[i]*YDIF*ZDIF;
	      }
	
	      double [][]T=	{ { RGX, RGXY, RGXZ },
					{ RGXY, RGY, RGYZ }, { RGXZ, RGYZ, RGZ } };

	      Matrix C=new Matrix(T);
						
	      EigenvalueDecomposition eid=C.eig();
	      //this.SortEID(eid);
	      
	      double [] eigenvalues=eid.getRealEigenvalues();
	      
//	      for (int i=0;i<=eigenvalues.length-1;i++) {
//	    	  System.out.println(eigenvalues[i]);
//	      }
	      	 
	      dd.MOMI1=eigenvalues[2];
	      dd.MOMI2=eigenvalues[1];
	      dd.MOMI3=eigenvalues[0];
	      
	      if (dd.MOMI2>1e-3) {
	    	  dd.MOMI4=dd.MOMI1/dd.MOMI2;
	      } else {
	    	  dd.MOMI4=-1;
	      }
	      
	      if (dd.MOMI3>1e-3) {
	    	  dd.MOMI5=dd.MOMI1/dd.MOMI3;
	      } else {
	    	  dd.MOMI5=-1;
	      }

	      if (dd.MOMI3>1e-3) {
	    	  dd.MOMI6=dd.MOMI2/dd.MOMI3;
	      } else {
	    	  dd.MOMI6=-1;
	      }

	      
	      
	      double RG=0;
	      
	      for (int i=0;i<=m.getAtomCount()-1;i++) {	
	    	  double XDIF = M.get(i,0) - CMX;
	    	  double YDIF = M.get(i,1) - CMY;
	    	  double ZDIF = M.get(i,2) - CMZ;
	    	  
	    	  double XSQ = XDIF*XDIF;
	    	  double YSQ = YDIF*YDIF;
	    	  double ZSQ = ZDIF*ZDIF;
	    	  
	    	  double RSQ=XSQ+YSQ+ZSQ;
	    	  RG+=wm[i]*RSQ;	    	  	    	  
	      }
	      
	      RG/=SUMWT;
	      RG=Math.sqrt(RG);
	      
	      dd.MOMI7=RG;
	      
	      
	      //System.out.println("RG="+RG);
	      
	      
//	      double RG2;
//	      
//	      if (eigenvalues[0]<1e-5) {
//	    	  RG2=Math.pow(eigenvalues[1]*eigenvalues[2],1.0/2.0);  
//	      } else {
//	    	  RG2=Math.pow(eigenvalues[0]*eigenvalues[1]*eigenvalues[2],1.0/3.0);
//	      }
//	      	      
//	      RG2/=SUMWT;
//	      RG2=Math.sqrt(RG2);
//	      
//	      System.out.println("RG2="+RG2);
	      
	      //Matrix Mnew=M.times(eid.getV());
	      //Mnew.print(6,4);
	      
	      
	      /*
	      double IA=0,IB=0,IC=0;
	      RGXY=0;
	      RGXZ=0;
	      RGYZ=0;
	      
	      for (int i=0;i<=m.getAtomCount()-1;i++) {	
	    	  
	    	  double XDIF = Mnew.get(i,0);
	    	  double YDIF = Mnew.get(i,1);
	    	  double ZDIF = Mnew.get(i,2);
	    	  
	    	  double XSQ = XDIF*XDIF;
	    	  double YSQ = YDIF*YDIF;
	    	  double ZSQ = ZDIF*ZDIF;
	    	  
	    	  
	    	  IA += wm[i]*(YSQ + ZSQ);
	    	  IB += wm[i]*(XSQ + ZSQ);
	    	  IC += wm[i]*(XSQ + YSQ);
	    	  
	          RGXY += wm[i]*XDIF*YDIF;
	          RGXZ += wm[i]*XDIF*ZDIF;
	          RGYZ += wm[i]*YDIF*ZDIF;
	    	  
	      }
	      System.out.println(IA+"\t"+IB+"\t"+IC);
	      System.out.println(RGXY+"\t"+RGXZ+"\t"+RGYZ);
	      */
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
	
	private void CalculateCenterOfMass(Matrix M,double [] wm) {
		// *** CALCULATE CENTER OF MASS
	      double XSUM=0.0;
	      double YSUM=0.0;
	      double ZSUM=0.0;
	      SUMWT=0.0;
	      
	      for (int i=0;i<=m.getAtomCount()-1;i++) {	     
	         XSUM+=M.get(i,0)*wm[i]; 
	         YSUM+=M.get(i,1)*wm[i];
	         ZSUM+=M.get(i,2)*wm[i];
	         SUMWT+=wm[i];
	      }
	      
	      CMX=XSUM/SUMWT;
	      CMY=YSUM/SUMWT;
	      CMZ=ZSUM/SUMWT;
	     // System.out.println(SUMWT);
	      
	     //System.out.println(CMX+"\t"+CMY+"\t"+CMZ);
	      
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
		
	
		wm=new double [m.getAtomCount()] ;
		
		AtomicProperties ap;
		
		try {
			ap=AtomicProperties.getInstance();
			
			for (int i=0;i<=m.getAtomCount()-1;i++) {
				String s=m.getAtom(i).getSymbol();
				wm[i]=ap.GetMass(s);
//				System.out.println(wm[i]);
			}
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		}
	
}
