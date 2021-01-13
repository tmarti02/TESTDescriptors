package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.lang.reflect.Field;
import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.PathFinder;

import java.util.LinkedList;

public class ChiDescriptor {

	DescriptorData dd;

	private IAtomContainer m;
	private double[] DV;
	private double[] D;
	private LinkedList[] paths;
	private IRingSet rs;

	public void Calculate(IAtomContainer m, DescriptorData dd, double[] D,
			double[] DV, LinkedList[] paths, IRingSet rs) {

		this.dd = dd;
		this.m = m;
		this.DV = DV;
		this.D = D;
		this.paths = paths;
		this.rs = rs;

		Initialize();
		this.Calculatex0();
		this.CalculatePathDescriptors();
		this.CalculateClusterDescriptors();
		this.CalculatePathClusterDescriptors();
		this.CalculateChainDescriptors();
		
		dd.knotp=dd.xc3-dd.xpc4;
		dd.knotpv=dd.xvc3-dd.xvpc4;
		//xc3 - xpc4
		
		
		/*
		if (debug)
			System.out.print("done\n");

		if (debug)
			System.out.print("finding clusters...");
					if (debug)
			System.out.print("done\n");

		if (debug)
			System.out.print("finding path clusters...");
		
		if (debug)
			System.out.print("done\n");
			*/
		
	}

	private void Initialize() {
		//		 initialize chi arrays:

		for (int i = 0; i <= 10; i++) {
			try {
				Field myField = dd.getClass().getField(dd.strxp[i]);
				myField.setDouble(dd, 0.0);

				myField = dd.getClass().getField(dd.strxvp[i]);
				myField.setDouble(dd, 0.0);

				if (i >= 3) {
					myField = dd.getClass().getField("xch" + i);
					myField.setDouble(dd, 0.0);

					myField = dd.getClass().getField("xvch" + i);
					myField.setDouble(dd, 0.0);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

		}

		dd.xc3 = 0;
		dd.xc4 = 0;
		dd.xpc4 = 0;

		dd.xvc3 = 0;
		dd.xvc4 = 0;
		dd.xvpc4 = 0;

	}

	private void Calculatex0() {
		//	 ***********************************************************************
		// ......CALCULATE CHI0
		for (int M = 0; M <= m.getAtomCount() - 1; M++) {
			dd.xv0 += 1 / Math.sqrt(Math.abs(DV[M]));
			dd.x0 += 1 / Math.sqrt(Math.abs(D[M]));
		}
	}

	private void CalculatePathDescriptors() {

		try {
			for (int M = 1; M <= 10; M++) {

				Field myField = dd.getClass().getField(dd.strxp[M]);
				myField.setDouble(dd, CalcChi(paths[M], D));
				
				//System.out.println(M+"\t"+paths[M].size());
				

				myField = dd.getClass().getField(dd.strxvp[M]);
				myField.setDouble(dd, CalcChi(paths[M], DV));

			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private double CalcChi(LinkedList MasterList, double[] D) {
		// in this method Chi is calculated for an LinkedList
		// which contains one or more array lists with the atom numbers

		double CHI = 0;

		//System.out.println(MasterList.size());
		
		for (int i = 0; i <= MasterList.size() - 1; i++) {

			LinkedList <Integer>al = (LinkedList) MasterList.get(i);

			double product = 1;
			
			for (int j = 0; j <= al.size() - 1; j++) {			
				int J = al.get(j);
				product *= D[J];
			}
			CHI += Math.pow(product, -0.5);

		}
		return CHI;
	}

	private double CalcChi2(LinkedList <Integer> al, double[] D) {
		// in this method chi is calculated for an LinkedList containing the atom
		// numbers
		double CHI = 0;

		double product = 1;

		for (int j = 0; j <= al.size() - 1; j++) {
//			int J = Integer.parseInt((String) al.get(j));
			int J = al.get(j);
			product *= D[J];
		}
		CHI += Math.pow(product, -0.5);

		return CHI;
	}

	private void CalculateClusterDescriptors() {

		for (int M = 3; M <= 4; M++) {

			LinkedList clusters = PathFinder.FindClusters(m, M);

			try {
				Field myField = dd.getClass().getField("xc" + M);
				myField.setDouble(dd, CalcChi(clusters, D));

				myField = dd.getClass().getField("xvc" + M);
				myField.setDouble(dd, CalcChi(clusters, DV));
			} catch (Exception e) {
				e.printStackTrace();
			}

		}
	}

	private void CalculatePathClusterDescriptors() {
		LinkedList clusters = PathFinder.FindPathClusters(m);
		dd.xpc4 = this.CalcChi(clusters, D);
		dd.xvpc4 = this.CalcChi(clusters, DV);

	}

	private void CalculateChainDescriptors() {

		for (int i = 0; i <= rs.getAtomContainerCount() - 1; i++) {
			IRing r = (IRing) rs.getAtomContainer(i);

			LinkedList al = new LinkedList();

			for (int j = 0; j <= r.getAtomCount() - 1; j++) {
				IAtom a = r.getAtom(j);				
				al.add(new Integer(m.getAtomNumber(a)));
			}

			try {

				if (r.getAtomCount() <= 10) {
					Field myField = dd.getClass().getField(
							"xch" + r.getAtomCount());
					double value = myField.getDouble(dd);
					value += this.CalcChi2(al, D);
					myField.setDouble(dd, value);

					myField = dd.getClass().getField(
							"xvch" + r.getAtomCount());
					value = myField.getDouble(dd);
					value += this.CalcChi2(al, DV);
					myField.setDouble(dd, value);
				}

			} catch (Exception e) {
				e.printStackTrace();
			}

		}

	}

	
	/*
	 * private void WriteChiValues() { try { fw.write("\n\nSimple\n"); // 9659
	 * 
	 * 
	 * for (int i=0;i<=10;i++) { fw.write("\tCHI"+i); }
	 * 
	 * fw.write("\n");
	 * 
	 * fw.write("P\t"); fw.write(myF.format(x0)+"\t");
	 * 
	 * for (int I=1;I<=10;I++) { Field
	 * myField=this.getClass().getField(strXP[I]);
	 * fw.write(myF.format(myField.getDouble(this))+"\t"); }
	 * 
	 * fw.write("\n");
	 * 
	 * fw.write("C\t\t\t\t"); // 4030 for (int I=3;I<=10;I++) {
	 * fw.write(myF.format(XC[I])+"\t"); } fw.write("\n");
	 * 
	 * fw.write("PC\t\t\t\t\t");
	 * 
	 * for (int I=4;I<=10;I++) { fw.write(myF.format(XPC[I])+"\t"); }
	 * fw.write("\n");
	 * 
	 * fw.write("CH\t\t\t\t");
	 * 
	 * for (int I=3;I<=10;I++) { fw.write(myF.format(XCH[I])+"\t"); }
	 * fw.write("\n");
	 * 
	 * 
	 * fw.write("\n\nValence\n"); // 9659
	 * 
	 * fw.write("P\t"); fw.write(myF.format(XPV[0])+"\t"); for (int I=1;I<=10;I++) {
	 * fw.write(myF.format(XPV[I])+"\t"); } fw.write("\n");
	 * 
	 * 
	 * 
	 * fw.write("C\t\t\t\t");
	 * 
	 * for (int I=3;I<=10;I++) { fw.write(myF.format(XCV[I])+"\t"); }
	 * fw.write("\n");
	 * 
	 * fw.write("PC\t\t\t\t\t");
	 * 
	 * for (int I=4;I<=10;I++) { fw.write(myF.format(XPCV[I])+"\t"); }
	 * fw.write("\n");
	 * 
	 * fw.write("CH\t\t\t\t");
	 * 
	 * for (int I=3;I<=10;I++) { fw.write(myF.format(XCHV[I])+"\t"); }
	 * fw.write("\n"); // goto 4500 } catch (Exception e) { } }
	 */

	

}//end class
