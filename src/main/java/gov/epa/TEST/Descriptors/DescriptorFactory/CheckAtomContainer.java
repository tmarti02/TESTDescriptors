package gov.epa.TEST.Descriptors.DescriptorFactory;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;

public class CheckAtomContainer {

	public static final String ERROR_CODE_STRUCTURE_ERROR = "SE";//TODO
	static final String ERROR_CODE_APPLICABILITY_DOMAIN_ERROR = "AD";
	static final String ERROR_CODE_DESCRIPTOR_CALCULATION_ERROR = "DE";
	private static final String ERROR_CODE_APPLICATION_ERROR = "AE";

	
	static void checkAtomContainer(AtomContainer m) {

		if (HaveBadElement(m)) {
			m.setProperty("Error", "Molecule contains unsupported element");
			m.setProperty("ErrorCode", ERROR_CODE_APPLICABILITY_DOMAIN_ERROR);
		} else if (m.getAtomCount() == 1) {
			m.setProperty("Error", "Only one nonhydrogen atom");
			m.setProperty("ErrorCode", ERROR_CODE_APPLICABILITY_DOMAIN_ERROR);
		} else if (m.getAtomCount() == 0) {
			m.setProperty("Error", "Number of atoms equals zero");
			m.setProperty("ErrorCode", ERROR_CODE_APPLICABILITY_DOMAIN_ERROR);
		} else if (!HaveCarbon(m)) {
			m.setProperty("Error", "Molecule does not contain carbon");
			m.setProperty("ErrorCode", ERROR_CODE_APPLICABILITY_DOMAIN_ERROR);
		}

		AtomContainerSet  moleculeSet = (AtomContainerSet)ConnectivityChecker.partitionIntoMolecules(m);
		if (moleculeSet.getAtomContainerCount() > 1) {
			//			m.setProperty("Error","Multiple molecules, largest fragment retained");
			m.setProperty("Error","Multiple molecules");
			m.setProperty("ErrorCode", ERROR_CODE_APPLICABILITY_DOMAIN_ERROR);
		}

		if (m.getProperty("Error")==null) m.setProperty("Error", "");

	}
	
static boolean HaveBadElement(IAtomContainer mol) {
		
		try {
						
			for (int i=0; i<mol.getAtomCount();i++) {
	
				String var = mol.getAtom(i).getSymbol();
	
				// OK: C, H, O, N, F, Cl, Br, I, S, P, Si, As, Hg, Sn
	
				if (!var.equals("C") && !var.equals("H") && !var.equals("O")
						&& !var.equals("N") && !var.equals("F")
						&& !var.equals("Cl") && !var.equals("Br")
						&& !var.equals("I") && !var.equals("S")
						&& !var.equals("P") && !var.equals("Si")
						&& !var.equals("As") && !var.equals("Hg")
						&& !var.equals("Sn")) {
	
					return true;
			
			
				}
			}
		
			
			return false;
	
		} catch (Exception e) {
			return true;
		}
		
		
	}

	static boolean HaveCarbon(IAtomContainer mol) {
		
		try {
			
		for (int i=0; i<mol.getAtomCount();i++) {
	
			String var = mol.getAtom(i).getSymbol();
	
			// OK: C, H, O, N, F, Cl, Br, I, S, P, Si, As, Hg, Sn
	
			if (var.equals("C")) {
				return true;
			}
		}
	
		return false;
	
	} catch (Exception e) {
		return true;
	}
		
	}
}
