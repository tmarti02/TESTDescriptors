package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.util.LinkedList;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

public class KappaDescriptor {
	/**
	 * 
	 * @param m - hydrogen suppressed molecule
	 * @param paths - paths in hydrogen suppressed molecule
	 * @param dd - descriptor object to store results
	 */
	public static void Calculate(IAtomContainer m,LinkedList[] paths,DescriptorData dd) {

		
		double A = m.getAtomCount();

		// kappa[0]=c*Math.log10(c);

		if (A == 1) {
			dd.k1 = 0;
			dd.k2 = 0;
			dd.k3 = 0;
		} else {

			double NumSinglePaths = paths[1].size();
			dd.k1 = (((A) * ((A - 1) * (A - 1))) / (NumSinglePaths * NumSinglePaths));

			if (A == 2) {
				dd.k2 = 0;
				dd.k3 = 0;
			} else {

				double NumDoublePaths = paths[2].size();

				dd.k2 = (((A - 1) * ((A - 2) * (A - 2))) / (NumDoublePaths * NumDoublePaths));
				if (A == 3) {
					dd.k3 = 0;
				} else {

					double NumTriplePaths = paths[3].size();

					if (A % 2 != 0) {
						dd.k3 = (((A - 1) * ((A - 3) * (A - 3))) / (NumTriplePaths * NumTriplePaths));
					} else {
						dd.k3 = (((A - 3) * ((A - 2) * (A - 2))) / (NumTriplePaths * NumTriplePaths));
					}
					if (NumTriplePaths == 0)
						dd.k3 = 0;

				}
			}
		}
		// System.out.println("kappa0 ="+kappa[0]);
		// System.out.println("kappa1 ="+kappa[1]);
		// System.out.println("kappa2 ="+kappa[2]);
		// System.out.println("kappa3 ="+kappa[3]);
		//ka3
	}

}
