package gov.epa.TEST.Descriptors.DescriptorCalculations;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;
import Jama.*;
import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.AtomicProperties;
import gov.epa.TEST.Descriptors.DescriptorUtilities.PathFinder;
import gov.epa.TEST.Descriptors.DescriptorUtilities.Utilities;

import java.util.*;
import java.lang.reflect.Field;
import java.text.DecimalFormat;

public class GetawayDescriptors {
	IAtomContainer m;
	DescriptorData dd;
	
	//Matrix M;
	//Matrix M_wo_Hydrogens;
	
	
	Matrix H;
	Matrix H_wo_Hydrogens;
	
	double [] wu;
	double [] wm;
	double [] wv;
	double [] we;
	double [] wp;
	int [][] DistanceWithH;
	int TD; // topological diameter
	double [][] R;
	
	
	public  void Calculate(IAtomContainer m,DescriptorData dd) {
		this.m=m;
		this.dd=dd;
		DistanceWithH=new int [m.getAtomCount()][m.getAtomCount()];
		
		PathFinder.CalculateDistanceMatrix(m,DistanceWithH);

		this.CalculateTopologicalDiameter();
		
		this.GetWeights();
		
		H=this.CalculateMolecularInfluenceMatrix(m);
		//H.print(6,4);
		
//		IMolecule m2=null;
//		try {
//		
//			m2=(Molecule)m.clone();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//		
//		CDKUtilities.RemoveHydrogens(m2);		
				
		this.CalculateITH2();
		this.CalculateHIC();
		
		this.CalculateHGM();
		this.CalculateHkw();
		this.CalculateHATSkw();
		
		this.CalculateInfluenceDistanceMatrix();
		this.CalculateRARS();
		this.CalculateRCON();
		this.CalculateREIG();
		this.CalculateRkw();
		this.CalculateRpluskw();
	}
	
	private void CalculateTopologicalDiameter() {
		this.TD=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			for (int j=i+1;j<=m.getAtomCount()-1;j++) {
				
//				System.out.println(i+"\t"+j+"\t"+DistanceWithH[i][j]);
				
				if (DistanceWithH[i][j]>this.TD) this.TD=DistanceWithH[i][j];
			}
		}
		
		
	}
	
	
	private void GetWeights() {
		wu=new double [m.getAtomCount()] ;
		wm=new double [m.getAtomCount()] ;
		wv=new double [m.getAtomCount()] ;
		we=new double [m.getAtomCount()] ;
		wp=new double [m.getAtomCount()] ;
		
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
				
			}
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		}
	
	private Matrix CalculateMolecularInfluenceMatrix(IAtomContainer m) {
	    double [][] Coords=new double [m.getAtomCount()][3];
		Matrix M=new Matrix(Coords);
		
		
		
	    double Xavg=0,Yavg=0,Zavg=0;
	    
	    // first make coords relative to the geometric center:
	    
	    for (int i=0;i<=m.getAtomCount()-1;i++){			
	    	IAtom a=m.getAtom(i);	    	
	    	Xavg+=a.getPoint3d().x;
	    	Yavg+=a.getPoint3d().y;
	    	Zavg+=a.getPoint3d().z;;		
	    }
		
	    Xavg/=m.getAtomCount();
	    Yavg/=m.getAtomCount();
	    Zavg/=m.getAtomCount();
	    
	   // System.out.println(Xavg+"\t"+Yavg+"\t"+Zavg);
	    
	    for (int i=0;i<=m.getAtomCount()-1;i++){
			IAtom a=m.getAtom(i);			
//	    	Coords[i][0]=a.getX3d();
//			Coords[i][1]=a.getY3d();
//			Coords[i][2]=a.getZ3d();
			
	    	Coords[i][0]=a.getPoint3d().x-Xavg;
			Coords[i][1]=a.getPoint3d().y-Yavg;
			Coords[i][2]=a.getPoint3d().z-Zavg;		

			
	    }
	    
	    		
		//M.print(6,4);
		
		double covxx=0,covxy=0,covxz=0;
		double covyx=0,covyy=0,covyz=0;
		double covzx=0,covzy=0,covzz=0;
				
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			covxx+=Coords[i][0]*Coords[i][0]; // Xavg=0 since used centered coords
			covyy+=Coords[i][1]*Coords[i][1];
			covzz+=Coords[i][2]*Coords[i][2];
			
			covxy+=Coords[i][0]*Coords[i][1];
			covxz+=Coords[i][0]*Coords[i][2];
			covyz+=Coords[i][1]*Coords[i][2];
			
		}
//		covxx/=(m.getAtomCount()-1);
//		covyy/=(m.getAtomCount()-1);
//		covzz/=(m.getAtomCount()-1);
//		
//		covxy/=(m.getAtomCount()-1);
//		covxz/=(m.getAtomCount()-1);
//		covyz/=(m.getAtomCount()-1);
		
		covxx/=(m.getAtomCount());
		covyy/=(m.getAtomCount());
		covzz/=(m.getAtomCount());
		
		covxy/=(m.getAtomCount());
		covxz/=(m.getAtomCount());
		covyz/=(m.getAtomCount());
		
		
		covyx=covxy;
		covzx=covxz;
		covzy=covyz;
		
		
		double[][] covmatrix = { { covxx, covxy, covxz },
				{ covyx, covyy, covyz }, { covzx, covzy, covzz } };
		
		Matrix C=new Matrix(covmatrix);
		
//		C.print(6,4);
		
		EigenvalueDecomposition eid=C.eig();
		
		Matrix eigenvectors=eid.getV();
		double [] eigenvalues=eid.getRealEigenvalues();
		
		
//		M.print(6,4);
		
		M=M.times(eigenvectors);// use PCA adjusted coordinates
		
		Coords=M.getArray(); //revise coordinate array with PCA adjusted array
		
		boolean [] delete=new boolean [3];
		
		int colcount=0;
		
		
		
		for (int i=0;i<=eigenvalues.length-1;i++) {
//			System.out.println(i+"\t"+eigenvalues[i]);
			if (eigenvalues[i]<1e-6) {
				delete[i]=true;
			} else {
				delete[i]=false;	
				colcount++;
			}						
		}
		

		
		double [][] newCoords=new double [m.getAtomCount()][colcount];
		// delete coordinate axes that have small eigenvalues 
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			
			if (!delete[0] && !delete[1] && !delete[2]) {
				newCoords[i][0]=Coords[i][0];
				newCoords[i][1]=Coords[i][1];
				newCoords[i][2]=Coords[i][2];
			} else if (!delete[0] && !delete[1] && delete[2]) {
				newCoords[i][0]=Coords[i][0];
				newCoords[i][1]=Coords[i][1];				
			} else if (!delete[0] && delete[1] && !delete[2]) {
				newCoords[i][0]=Coords[i][0];
				newCoords[i][1]=Coords[i][2];								
			} else if (delete[0] && !delete[1] && !delete[2]) {
				newCoords[i][0]=Coords[i][1];
				newCoords[i][1]=Coords[i][2];
			} else if (!delete[0] && delete[1] && delete[2]) {
				newCoords[i][0]=Coords[i][0];				
			} else if (delete[0] && !delete[1] && delete[2]) {
				newCoords[i][0]=Coords[i][1];				
			} else if (delete[0] && delete[1] && !delete[2]) {
				newCoords[i][0]=Coords[i][2];				
			}
		}
		
		
		//M.print(6,4);
		
		M=new Matrix(newCoords); // use new matrix with deleted coordinate axes
//		M=new Matrix(Coords);
		
		
//		M.print(6,4);
						
		
		//Matrix M2=eigenvectors.transpose().times(M.transpose());
		//Matrix eigenvectors=M.eig().getV();
//		System.out.println("eigenvectors:");
//		eigenvectors.print(6,4);
		
//		M=M.times(eigenvectors);		
		//M.print(6,4);
		
//		double [] eigenvalues=M.eig().getRealEigenvalues();
//		for (int i=0;i<=eigenvalues.length-1;i++) {
//			System.out.println(eigenvalues[i]);
//			
//		}
		
		
		
		H=M.times(M.transpose().times(M).inverse()).times(M.transpose());
//		System.out.println("H=");
//		H.print(6,4);
		
//		Matrix MT=M.transpose();
//		Matrix MTM=MT.times(M);
//		Matrix invMTM=MTM.inverse();
		
//		invMTM.print(6,4);
				
//		Matrix MinvMTM=M.times(invMTM);
//		MinvMTM.print(6,4);
//		
//		
//		Matrix H2=MinvMTM.times(MT);
//		System.out.println("H2=");
//		H2.print(6,4);
//		
//		Matrix invMTM_MT=invMTM.times(MT);
//		Matrix H3=M.times(invMTM_MT);
//		System.out.println("H3=");
//		H3.print(6,4);
		
		
		return H;
	}

	
	
	private void CalculateHGM() {
		
		dd.HGM=1;
		double count=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++){
			// Atoms with a leverage value less than or equal to 10-6 are
			// excluded from the calculation of HGM.
			if (H.get(i,i)>1e-6) { 
				dd.HGM*=H.get(i,i);
				//System.out.println(H.get(i,i));
				count++; 
			}						
		}
		
		dd.HGM=Math.pow(dd.HGM,1.0/count);
		dd.HGM*=100.0;
		//System.out.println("HGM="+dd.HGM);
		
		
		
	}
	
	
//	private void CalculateITH() {
//		
//		//TODO : figure out how to consistently match values in Consonni et al., 2002
//		// possible soln = redo H matrix without H's attached
//		
//		
//		
////		Molecule m2=(Molecule)m.clone();		
////		CDKUtilities.RemoveHydrogens(m2);
////		int A0=m2.getAtomCount();
////		
////		System.out.println("A0="+A0);
//		
//		int A0=0;
////		ArrayList al=new ArrayList();
//		
//		Hashtable ht=new Hashtable();
//				
//		DecimalFormat df=new DecimalFormat("0.0000"); //allows 4-digit equivalence of leverage values
//		
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			String s=m.getAtom(i).getSymbol();
//			if (!(s.equals("H"))) {
//				A0++;
//				String strnum=df.format(H.get(i,i));				
//				ht.put(strnum,new Integer(0));								
//			}
//		}
//		
//										
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			String s=m.getAtom(i).getSymbol();
//			if (!(s.equals("H"))) {				
//				String strnum=df.format(H.get(i,i));				
//				int count=(Integer)ht.get(strnum);
//				
//				count++;
//				ht.put(strnum,count);
//			}
//		}
//		
//		
//		double SumNlogN=0.0;
//		for (int i=0;i<=m.getAtomCount()-1;i++) {
//			String s=m.getAtom(i).getSymbol();
//			if (!(s.equals("H"))) {				
//				String strnum=df.format(H.get(i,i));				
//				int count=(Integer)ht.get(strnum);
//				
//				//System.out.println(strnum+"\t"+count);
//				SumNlogN+=count*Log(2,count);
//			}
//			
//		}
//		
//		double ITH=A0*Log(2,A0)-SumNlogN;
//		double ISH=ITH/(A0*Log(2,A0));
//		
//		System.out.println("ITH="+ITH);
//		System.out.println("ISH="+ISH);
//		
//		
//	}
	
	private void CalculateITH2() {
		
		
		int A0=0;
		ArrayList al=new ArrayList();
		
		double tol=0.0001;
										
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			String s=m.getAtom(i).getSymbol();
			
//			System.out.println(s+"\t"+H.get(i,i));
			
			if (!(s.equals("H"))) {				
				A0++;

				if (al.size()==0) {
					al.add(H.get(i,i)+"\t"+1);
				} else {
					this.CheckForMatches(al,H.get(i,i),tol);
				}

			}									
			
			
		}
		
		
		double SumNlogN=0.0;
		
		for (int i=0;i<=al.size()-1;i++) {
			String strLine=(String)al.get(i);
			List l=Utilities.Parse(strLine,"\t");
			String strval=(String)l.get(0);
			String strcount=(String)l.get(1);
			double count=Double.parseDouble(strcount);
			
//			System.out.println(count+"\t"+strval);
			SumNlogN+=count*Log(2,count);
			
			
		}
		
		
		
		
//		System.out.println("AO="+A0);
//		System.out.println("Number of atoms="+dd.nAT);
//		System.out.println("A0logA0="+A0*Log(2,A0));
//		System.out.println("SumNlogN="+SumNlogN);
		
		dd.ITH=A0*Log(2,A0)-SumNlogN;
		dd.ISH=dd.ITH/(A0*Log(2,A0));
		
		//System.out.println("ITH="+dd.ITH);
		//System.out.println("ISH="+dd.ISH);
		
		
	}
	
	private void CheckForMatches(ArrayList al,double val,double tol) {
		
		for (int i=0;i<=al.size()-1;i++) {
			String strLine=(String)al.get(i);
			List l=Utilities.Parse(strLine,"\t");
			String strval=(String)l.get(0);
			String strcount=(String)l.get(1);
			
			double value=Double.parseDouble(strval);
			int count=Integer.parseInt(strcount);
									
			if (Math.abs(val-value)<=tol) {
				count++;
				al.set(i,val+"\t"+count);
				return;
			}
			
		}
		al.add(val+"\t"+1);
		
	}

	private void CalculateHIC() {
		double D=H.trace();
		
		//System.out.println("D="+D);
		
		double sum=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			sum+=H.get(i,i)/D*Log(2,H.get(i,i)/D);
		}
		dd.HIC=-sum;
		
		//System.out.println("HIC="+dd.HIC);
	}
	
	
//	private void CheckForMatches(ArrayList al,double val) {
//		for (int i=0;i<=al.size()-1;i++) {
//		}
//	}
	
	
	private void CalculateHkw() {
		
		
		double [] Hu=this.CalculateHvalues(wu);
		double [] Hm=this.CalculateHvalues(wm);
		double [] Hv=this.CalculateHvalues(wv);
		double [] He=this.CalculateHvalues(we);
		double [] Hp=this.CalculateHvalues(wp);
		
		//TODO: assign these to dd class
		try {
		for (int i=0;i<=8;i++) {
	
				Field myField = dd.getClass().getField("H"+i+"u");
				myField.setDouble(dd, Hu[i]);			
				
				myField = dd.getClass().getField("H"+i+"m");
				myField.setDouble(dd, Hm[i]);			
				
				myField = dd.getClass().getField("H"+i+"v");
				myField.setDouble(dd, Hv[i]);			
				
				myField = dd.getClass().getField("H"+i+"p");
				myField.setDouble(dd, Hp[i]);			
				
				myField = dd.getClass().getField("H"+i+"e");
				myField.setDouble(dd, He[i]);			
				
			}
	
		
		dd.HTu=Hu[9];
		dd.HTm=Hm[9];
		dd.HTv=Hv[9];
		dd.HTp=Hp[9];
		dd.HTe=He[9];
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		

				
	}
	
	private void CalculateHATSkw() {
		
		
		double [] HATSu=this.CalculateHATSvalues(wu);
		double [] HATSm=this.CalculateHATSvalues(wm);
		double [] HATSv=this.CalculateHATSvalues(wv);
		double [] HATSp=this.CalculateHATSvalues(wp);
		double [] HATSe=this.CalculateHATSvalues(we);
		
		//TODO: assign these to dd class
		try {
		
									
			for (int i=0;i<=8;i++) {
								
				Field myField = dd.getClass().getField("HATS"+i+"u");
				myField.setDouble(dd, HATSu[i]);			
				
				myField = dd.getClass().getField("HATS"+i+"m");
				myField.setDouble(dd, HATSm[i]);			
				
				myField = dd.getClass().getField("HATS"+i+"v");
				myField.setDouble(dd, HATSv[i]);			
				
				myField = dd.getClass().getField("HATS"+i+"p");
				myField.setDouble(dd, HATSp[i]);			
				
				myField = dd.getClass().getField("HATS"+i+"e");
				myField.setDouble(dd, HATSe[i]);			
				
			}
			
			dd.HATSu=HATSu[9];
			dd.HATSm=HATSm[9];
			dd.HATSv=HATSv[9];
			dd.HATSp=HATSp[9];
			dd.HATSe=HATSe[9];
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		

				
	}
	
	
	private double [] CalculateHvalues(double [] w) {
		
		double [] Hw=new double [TD+1];
		double [] Hw2=new double [10];
		
		Hw[0]=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++){			
			Hw[0]+=H.get(i,i)*w[i]*w[i];
		}
		
		double sumH=0;
		for (int k=1;k<=TD;k++) {
			Hw[k]=0;
			
		//	System.out.println(H.getRowDimension());
//			System.out.println(H.getColumnDimension());
	//		System.out.println(m.getAtomCount());
			
			for (int i=0;i<=m.getAtomCount()-2;i++) {
				for (int j=i+1;j<=m.getAtomCount()-1;j++) {
					if (H.get(i,j) > 0 && DistanceWithH[i][j]==k) {
						Hw[k]+=H.get(i,j)*w[i]*w[j];
					}
				}
			}
			
			sumH+=Hw[k];
		}
		
		
		if (TD>=8) {
			for (int k=0;k<=8;k++) {
				Hw2[k]=Hw[k];
			}
		} else {
			for (int k=0;k<=TD;k++) {
				Hw2[k]=Hw[k];
			}
			for (int k=TD+1;k<=8;k++) {
				Hw2[k]=0;
			}
		}
		
		Hw2[9]=Hw[0]+2*sumH;
		
		//System.out.println(Hw2[9]);
		return Hw2;
		
	}
	
private double [] CalculateHATSvalues(double [] w) {
		
		double [] HATSw=new double [TD+1];
		double [] HATSw2=new double [10];
		
		
		HATSw[0]=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++){			
			HATSw[0]+=Math.pow(H.get(i,i)*w[i],2);
		}
		
		double sumHATS=0;
		

		
		for (int k=1;k<=TD;k++) {
			HATSw[k]=0;
			
		//	System.out.println(H.getRowDimension());
//			System.out.println(H.getColumnDimension());
	//		System.out.println(m.getAtomCount());
			
			for (int i=0;i<=m.getAtomCount()-2;i++) {
				for (int j=i+1;j<=m.getAtomCount()-1;j++) {
					if (DistanceWithH[i][j]==k) {						
						HATSw[k]+=w[i]*H.get(i,i)*w[j]*H.get(j,j);
					}
				}
			}
			sumHATS+=HATSw[k];
		}
		
		if (TD>=8) {
			for (int k=0;k<=8;k++) {
				HATSw2[k]=HATSw[k];
			}
		} else {
			for (int k=0;k<=TD;k++) {
				HATSw2[k]=HATSw[k];
			}
			for (int k=TD+1;k<=8;k++) {
				HATSw2[k]=0;
			}
		}
		
		HATSw2[9]=HATSw[0]+2*sumHATS;
		
		return HATSw2;
		
	}
	


	private void CalculateInfluenceDistanceMatrix() {
		R=new double[m.getAtomCount()][m.getAtomCount()];
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			R[i][i]=0;
			for (int j=i+1;j<=m.getAtomCount()-1;j++) {
				R[i][j]=Math.sqrt(H.get(i,i)*H.get(j,j))/Calculaterij(i,j);								
				R[j][i]=R[i][j];				
			}
		}
//		Matrix R2=new Matrix(R);
//		R2.print(6,4);
		
		
	}

	private double Calculaterij(int i,int j) {
		double rij=0;
		
		IAtom atomi=m.getAtom(i);
		IAtom atomj=m.getAtom(j);
		
		double term1=atomi.getPoint3d().x-atomj.getPoint3d().x;
		double term2=atomi.getPoint3d().y-atomj.getPoint3d().y;
		double term3=atomi.getPoint3d().z-atomj.getPoint3d().z;
		
		
		rij=Math.sqrt(term1*term1+term2*term2+term3*term3);
		
		
		return rij;
		
	}
	
	private void CalculateRARS() {
		dd.RARS=0;
		for (int i=0;i<=m.getAtomCount()-1;i++) {		
			for (int j=0;j<=m.getAtomCount()-1;j++) {
				dd.RARS+=R[i][j];
			}
		}
		
		dd.RARS/=m.getAtomCount();
	}
	
	private void CalculateRCON() {
		dd.RCON=0;		
		
		
		for (int b=0;b<m.getBondCount();b++) {
			int i=m.getAtomNumber(m.getBond(b).getAtom(0));
			int j=m.getAtomNumber(m.getBond(b).getAtom(1));			
			
			double RSi=0;
			for (int jj=0;jj<=m.getAtomCount()-1;jj++) {
				RSi+=R[i][jj];
			}
			
			double RSj=0;
			for (int jj=0;jj<=m.getAtomCount()-1;jj++) {
				RSj+=R[j][jj];
			}
			
			dd.RCON+=Math.sqrt(RSi*RSj);
						
		}
	}
	
	private void CalculateREIG() {
		Matrix Rmat=new Matrix(R);
		double [] eigenvalues=Rmat.eig().getRealEigenvalues();
		
		dd.REIG=eigenvalues[eigenvalues.length-1];
		
//		for (int i=0;i<=eigenvalues.length-1;i++) {
//			System.out.println(eigenvalues[i]);
//		}
		
		
	}
	
	
	private void CalculateRkw() {
		
		
		double [] Ru=this.CalculateRvalues(wu);
		double [] Rm=this.CalculateRvalues(wm);
		double [] Rv=this.CalculateRvalues(wv);
		double [] Rp=this.CalculateRvalues(wp);
		double [] Re=this.CalculateRvalues(we);
		
		//TODO: assign these to dd class
		try {
		
									
			for (int i=1;i<=8;i++) {
								
				Field myField = dd.getClass().getField("R"+i+"u");
				myField.setDouble(dd, Ru[i]);			
				
				myField = dd.getClass().getField("R"+i+"m");
				myField.setDouble(dd, Rm[i]);			
				
				myField = dd.getClass().getField("R"+i+"v");
				myField.setDouble(dd, Rv[i]);			
				
				myField = dd.getClass().getField("R"+i+"p");
				myField.setDouble(dd, Rp[i]);			
				
				myField = dd.getClass().getField("R"+i+"e");
				myField.setDouble(dd, Re[i]);			
				
			}
			
			dd.RTu=Ru[9];
			dd.RTm=Rm[9];
			dd.RTv=Rv[9];
			dd.RTp=Rp[9];
			dd.RTe=Re[9];
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
private void CalculateRpluskw() {
		
		
		double [] Ruplus=this.CalculateRplusvalues(wu);
		double [] Rmplus=this.CalculateRplusvalues(wm);
		double [] Rvplus=this.CalculateRplusvalues(wv);
		double [] Rpplus=this.CalculateRplusvalues(wp);
		double [] Replus=this.CalculateRplusvalues(we);
		
		//TODO: assign these to dd class
		try {
		
									
			for (int i=1;i<=8;i++) {
								
				Field myField = dd.getClass().getField("R"+i+"uplus");
				myField.setDouble(dd, Ruplus[i]);			
				
				myField = dd.getClass().getField("R"+i+"mplus");
				myField.setDouble(dd, Rmplus[i]);			
				
				myField = dd.getClass().getField("R"+i+"vplus");
				myField.setDouble(dd, Rvplus[i]);			
				
				myField = dd.getClass().getField("R"+i+"pplus");
				myField.setDouble(dd, Rpplus[i]);			
				
				myField = dd.getClass().getField("R"+i+"eplus");
				myField.setDouble(dd, Replus[i]);			
				
			}
			
			dd.RTuplus=Ruplus[9];
			dd.RTmplus=Rmplus[9];
			dd.RTvplus=Rvplus[9];
			dd.RTpplus=Rpplus[9];
			dd.RTeplus=Replus[9];
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private double [] CalculateRvalues(double [] w) {
		
		double [] Rw=new double [TD+1];
		double [] Rw2=new double [10];
		
		
		Rw[0]=0;
		
		
		double sumR=0;
		
		for (int k=1;k<=TD;k++) {
			Rw[k]=0;
						
			for (int i=0;i<=m.getAtomCount()-2;i++) {
				for (int j=i+1;j<=m.getAtomCount()-1;j++) {
					if (DistanceWithH[i][j]==k) {
						Rw[k]+=R[i][j]*w[i]*w[j];
					}
				}
			}
			sumR+=Rw[k];
		}
		
		if (TD>=8) {
			for (int k=0;k<=8;k++) {
				Rw2[k]=Rw[k];
			}
		} else {
			for (int k=0;k<=TD;k++) {
				Rw2[k]=Rw[k];
			}
			for (int k=TD+1;k<=8;k++) {
				Rw2[k]=0;
			}
		}
		
		Rw2[9]=2*sumR;
		
		return Rw2;
		
	}
	
private double [] CalculateRplusvalues(double [] w) {
		
		double [] Rwplus=new double [TD+1];
		double [] Rwplus2=new double [10];
		
		
		Rwplus[0]=0;
		
		
		double overallMax=0;
		
		for (int k=1;k<=TD;k++) {
			Rwplus[k]=0;
						
			for (int i=0;i<=m.getAtomCount()-2;i++) {
				for (int j=i+1;j<=m.getAtomCount()-1;j++) {
					if (DistanceWithH[i][j]==k) {
						double val=R[i][j]*w[i]*w[j];
						if (val>Rwplus[k]) Rwplus[k]=val;
					}
				}
			}			
			
			if (Rwplus[k]>overallMax) overallMax=Rwplus[k];
		}
		
		if (TD>=8) {
			for (int k=0;k<=8;k++) {
				Rwplus2[k]=Rwplus[k];
			}
		} else {
			for (int k=0;k<=TD;k++) {
				Rwplus2[k]=Rwplus[k];
			}
			for (int k=TD+1;k<=8;k++) {
				Rwplus2[k]=0;
			}
		}
		
		Rwplus2[9]=overallMax;
		
		return Rwplus2;
		
	}
	
	
	
	private static double Log(int base,double x) {
		
		double Logbx=0;
		
		Logbx=Math.log10(x)/Math.log10((double)base);
		
		return Logbx;
		
		
	}	
}
