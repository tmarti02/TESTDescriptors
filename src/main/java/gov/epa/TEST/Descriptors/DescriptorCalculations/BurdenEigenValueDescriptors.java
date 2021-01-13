package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.lang.reflect.Field;
import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.ringsearch.AllRingsFinder;

import Jama.*;
//import ToxPredictor.Utilities.CDKUtilities;
//import ToxPredictor.Utilities.HueckelAromaticityDetector;
import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.AtomicProperties;

//import ToxPredictor.Utilities.CDKUtilities;

//They are derived for a hydrogen-included molecular graph 
// and calculated from the Burden matrix, which is a modified 
// adjacency matrix defined as the following: the diagonal elements
// are atomic properties; the off-diagonal elements corresponding 
// to pairs of bonded atoms are the square roots of conventional 
// bond order; all other matrix elements are set at 0.001.

public class BurdenEigenValueDescriptors {

	/**
	 * @param args
	 */

	// IsotopeFactory elfac = null;
	public void Calculate(IAtomContainer m, IRingSet rs,DescriptorData dd) {
		
		double[][] Bm = new double[m.getAtomCount()][m.getAtomCount()];
		double[][] Bv = new double[m.getAtomCount()][m.getAtomCount()];
		double[][] Be = new double[m.getAtomCount()][m.getAtomCount()];
		double[][] Bp = new double[m.getAtomCount()][m.getAtomCount()];

		AtomicProperties ap = null;

		try {
			ap = AtomicProperties.getInstance();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

		try {
			
			
			for (int i = 0; i < m.getAtomCount(); i++) {
				for (int j = 0; j < m.getAtomCount(); j++) {
					if (i == j) {
						String symbol = m.getAtom(i).getSymbol();
						Bm[i][i] = ap.GetNormalizedMass(symbol);
						Bv[i][i] = ap.GetNormalizedVdWVolume(symbol);
						Be[i][i] = ap.GetNormalizedElectronegativity(symbol);
						Bp[i][i] = ap.GetNormalizedPolarizability(symbol);
						
//						System.out.println(Bm[i][i]);

					} else {

						IAtom atomi = m.getAtom(i);
						IAtom atomj = m.getAtom(j);

						IBond b = m.getBond(atomi, atomj);
						
						if (b instanceof IBond) {//need IBond 5/4/17

							double val = 0;

							if (EStateFragmentDescriptor.InSameAromaticRing(m,
									atomi, atomj, rs))

								val = Math.sqrt(1.5);
							else
								val = Math.sqrt(b.getOrder().numeric());

							Bm[i][j] = val;// see pg 310
							// Burden, 1997
							Bv[i][j] = val;
							Be[i][j] = val;
							Bp[i][j] = val;
							
							// List cai = m.getConnectedAtomsList(atomi);
							// List caj = m.getConnectedAtomsList(atomj);

							// Elements of B corresponding to bonds to
							// terminal atoms (i.e. atoms with one connection
							// only) are augmented by 0.01

							// if (cai.length == 1 || caj.length == 1) {
							// Bm[i][j] += 0.01; //add 0.1 to match Dragon 6.0?
							// }

						} else {
							Bm[i][j] = 0.001;
							Bv[i][j] = 0.001;
							Be[i][j] = 0.001;
							Bp[i][j] = 0.001;

						}
					}// end i!=j else

				}// end j for loop
			}// end i for loop

		} catch (Exception exc) {
			exc.printStackTrace();
		}

		Matrix matrix;
		EigenvalueDecomposition ed;
		double[] eigenvalues;

		matrix = new Matrix(Bm);
		ed = new EigenvalueDecomposition(matrix);
		eigenvalues = ed.getRealEigenvalues();
		AssignDescriptors("m", eigenvalues, dd);

		// for (int i=0;i<=eigenvalues.length-1;i++) {
		// System.out.println(eigenvalues[i]);
		// }

		matrix = new Matrix(Bv);
		ed = new EigenvalueDecomposition(matrix);
		eigenvalues = ed.getRealEigenvalues();
		AssignDescriptors("v", eigenvalues, dd);

		matrix = new Matrix(Be);
		ed = new EigenvalueDecomposition(matrix);
		eigenvalues = ed.getRealEigenvalues();
		AssignDescriptors("e", eigenvalues, dd);

		matrix = new Matrix(Bp);
		ed = new EigenvalueDecomposition(matrix);
		eigenvalues = ed.getRealEigenvalues();
		AssignDescriptors("p", eigenvalues, dd);

		// Matrix bob =ed.getD();
		// bob.print(new java.text.DecimalFormat("0.0000"),10);

		// CDKUtilities.RemoveHydrogens(m);

	}

	private static void AssignDescriptors(String strvar, double[] eigenvalues,
			DescriptorData dd) {
		int poscount = 0, negcount = 0;
		for (int i = 0; i <= eigenvalues.length - 1; i++) {
			if (eigenvalues[i] < 0)
				negcount++;
		}
		poscount = eigenvalues.length - negcount;

		// System.out.println("negcount="+negcount);

		try {

			if (negcount <= 8) {
				for (int i = 0; i <= negcount - 1; i++) {
					String fieldname = "BEL" + strvar + (i + 1);
					Field myField = dd.getClass().getField(fieldname);
					myField.setDouble(dd, Math.abs(eigenvalues[i]));
					// System.out.println(fieldname+"\t"+myField.getDouble(dd));
				}
				for (int i = negcount; i <= 8 - 1; i++) {
					String fieldname = "BEL" + strvar + (i + 1);
					Field myField = dd.getClass().getField(fieldname);
					myField.setDouble(dd, 0.0);
					// System.out.println(fieldname+"\t"+myField.getDouble(dd));

				}
			} else {
				for (int i = 0; i <= 8 - 1; i++) {
					String fieldname = "BEL" + strvar + (i + 1);
					Field myField = dd.getClass().getField(fieldname);
					myField.setDouble(dd, Math.abs(eigenvalues[i]));
					// System.out.println(fieldname+"\t"+myField.getDouble(dd));
				}
			}

			if (poscount <= 8) {
				for (int i = 1; i <= poscount; i++) {
					String fieldname = "BEH" + strvar + (i);
					Field myField = dd.getClass().getField(fieldname);
					myField.setDouble(dd, eigenvalues[eigenvalues.length - i]);
				}

				for (int i = poscount + 1; i <= 8; i++) {
					String fieldname = "BEH" + strvar + (i);
					Field myField = dd.getClass().getField(fieldname);
					myField.setDouble(dd, 0.0);
					// System.out.println(fieldname+"\t"+myField.getDouble(dd));

				}

			} else {
				for (int i = 1; i <= 8; i++) {
					String fieldname = "BEH" + strvar + (i);
					Field myField = dd.getClass().getField(fieldname);
					myField.setDouble(dd, eigenvalues[eigenvalues.length - i]);
				}
			}

		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}



	private IRingSet FindRings(IAtomContainer m) {
		IRingSet rs = null;

		try {
			AllRingsFinder arf = new AllRingsFinder();
			arf.setTimeout(100000);
			rs = arf.findAllRings(m);

		} catch (Exception e) {
			e.printStackTrace();
		}

		return rs;
	}

}
