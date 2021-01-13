package gov.epa.TEST.Descriptors.DescriptorUtilities;

import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.formula.MolecularFormula;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import com.epam.indigo.Indigo;
import com.epam.indigo.IndigoException;
import com.epam.indigo.IndigoInchi;
import com.epam.indigo.IndigoObject;


import net.sf.jniinchi.INCHI_OPTION;
import net.sf.jniinchi.INCHI_RET;

public class CDKUtilities {

	private static final Logger logger = LogManager.getLogger(CDKUtilities.class);

	public static String generateSmiles(IAtomContainer ac) {
		return generateSmiles(ac, SmiFlavor.Unique);
	}

	public static String generateSmiles(IAtomContainer ac, int flavor) {
		try {
			SmilesGenerator sg = new SmilesGenerator(flavor);
			String smiles = sg.create(ac);
			return smiles;

		} catch (Exception ex) {
			logger.catching(ex);
		}
		return null;
	}

	
	
	
	
	public static Inchi generateInChiKey(IAtomContainer ac) {
		String warning = "";

		try {

			// Generate factory - throws CDKException if native code does not load
//			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
//			// Get InChIGenerator
//			InChIGenerator gen = factory.getInChIGenerator(ac);

			List<INCHI_OPTION> options = new ArrayList<INCHI_OPTION>();
			// FIXME: uncomment after updating DB
//			options.add(INCHI_OPTION.FixedH);//makes sure  tautomers come out different! 

			// TODO- do we need any of these options???
			// https://github.com/cdk/cdk/issues/253
//			options.add(INCHI_OPTION.SAbs);
//			options.add(INCHI_OPTION.SAsXYZ);
//			options.add(INCHI_OPTION.SPXYZ);
//			options.add(INCHI_OPTION.FixSp3Bug);
//			options.add(INCHI_OPTION.AuxNone);
			

			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			InChIGenerator gen = factory.getInChIGenerator(ac, options);

			INCHI_RET ret = gen.getReturnStatus();

			if (ret == INCHI_RET.WARNING) {
				// InChI generated, but with warning message
//				System.out.println();

				warning = "InChI warning: " + gen.getMessage();

			} else if (ret != INCHI_RET.OKAY) {
				// InChI generation failed
				warning = "InChI failed: " + ret.toString() + " [" + gen.getMessage() + "]";
			}

			Inchi inchi=new Inchi();
			inchi.inchi = gen.getInchi();
			inchi.inchiKey = gen.getInchiKey();
			inchi.warning=warning;
			

			return inchi;

			// TODO: distinguish between singlet and undefined spin multiplicity
			// TODO: double bond and allene parities
			// TODO: problem recognising bond stereochemistry

		} catch (CDKException | IllegalArgumentException ex) {
			Inchi inchi=new Inchi();
			inchi.warning=ex.getMessage();					
			return inchi;
		}
	}

	public static AtomContainer getAtomContainer(String inchi) {
		// Generate factory - throws CDKException if native code does not load
		InChIGeneratorFactory factory;
		try {
			factory = InChIGeneratorFactory.getInstance();
			
			// Get InChIToStructure
			InChIToStructure intostruct = factory.getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());

			INCHI_RET ret = intostruct.getReturnStatus();
			if (ret == INCHI_RET.WARNING) {
			// Structure generated, but with warning message
			System.out.println("InChI warning: " + intostruct.getMessage());
			} else if (ret != INCHI_RET.OKAY) {
			// Structure generation failed
			throw new CDKException("Structure generation failed failed: " + ret.toString()
			+ " [" + intostruct.getMessage() + "]");
			}

			AtomContainer container = (AtomContainer) intostruct.getAtomContainer();
			return container;

		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}
	
	
	public static String FixSmiles(String Smiles) {
		Smiles = Smiles.replace("CL", "Cl");
		Smiles = Smiles.replace("(H)", "([H])");
		Smiles = Smiles.replace("N(=O)(=O)", "[N+](=O)([O-])");

//		Smiles=Smiles.replace("O=n","[O-][n+]");		
//		Smiles=Smiles.replace("n(=O)","[n+][O-]");
//		Smiles=Smiles.replace("n=O","[n+][O-]");		
//		Smiles=Smiles.replace("n1=O","[n+][O-]");
//		Smiles=Smiles.replace("n1(=O)","[n+][O-]");
		//
		Smiles = Smiles.replace("N=N#N", "N=[N+]=[N-]");
//		Smiles=Smiles.replace("#N=O","#[N+][O-]");
		Smiles = Smiles.trim();

		return Smiles;

	}

	//
	/**
	 * Add explicit hydrogens using CDK's hydrogen adder
	 * 
	 * @param atomContainer
	 * @return
	 */
	public static AtomContainer addHydrogens(IAtomContainer atomContainer) {
		AtomContainer ac;
		try {
			ac = (AtomContainer) atomContainer.clone();
			AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(ac);
			CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(ac.getBuilder());
			hAdder.addImplicitHydrogens(ac);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(ac);

//            for (int i=0;i<ac.getAtomCount();i++) {
//            	System.out.println(i+"\t"+ac.getAtom(i).getSymbol());
//            }

			return ac;
		} catch (Exception e) {
			return null;
		}
	}

	/**
	 * Convert charged nitrogens to uncharged version using brute force rather than
	 * normalizer
	 * 
	 * @param m
	 * @return
	 */
	public static boolean FixNitroGroups2(IAtomContainer m) {
		// changes nitros given by [N+](=O)[O-] to N(=O)(=O)
		boolean changed = false;
		try {
			for (int i = 0; i <= m.getAtomCount() - 1; i++) {
				IAtom a = m.getAtom(i);
				if (a.getSymbol().equals("N")) {
					List<IAtom> ca = m.getConnectedAtomsList(a);

					if (ca.size() == 3) {

						IAtom[] cao = new Atom[2];

						int count = 0;

						for (int j = 0; j <= 2; j++) {
							IAtom caj = (IAtom) ca.get(j);
							if (caj.getSymbol().equals("O")) {
								count++;
							}
						}

						if (count > 1) {

							count = 0;
							for (int j = 0; j <= 2; j++) {
								IAtom caj = (IAtom) ca.get(j);
								if (caj.getSymbol().equals("O")) {
									if (m.getConnectedAtomsCount(caj) == 1) {// account for possibility of ONO2
										cao[count] = caj;
										count++;
									}
								}
							}

							IBond.Order order1 = m.getBond(a, cao[0]).getOrder();
							IBond.Order order2 = m.getBond(a, cao[1]).getOrder();

							// int totalobonds=0;
							// totalobonds+=m.getBond(a,cao[0]).getOrder();
//						totalobonds+=m.getBond(a,cao[1]).getOrder();

							// if (totalobonds==4) { // need to fix

							boolean cond1 = order1.equals(IBond.Order.SINGLE) && order2.equals(IBond.Order.DOUBLE);
							boolean cond2 = order2.equals(IBond.Order.SINGLE) && order1.equals(IBond.Order.DOUBLE);

							if (cond1 || cond2) {
								a.setFormalCharge(0);
								cao[0].setFormalCharge(0); // pick first O arbitrarily
								cao[1].setFormalCharge(0); // pick first O arbitrarily
								m.getBond(a, cao[0]).setOrder(IBond.Order.DOUBLE);
								m.getBond(a, cao[1]).setOrder(IBond.Order.DOUBLE);
								changed = true;
							}
						} // end if count>1

					} // end ca==3 if

				} // end symbol == N

			}

			return changed;
		} catch (Exception e) {
			return changed;
		}
	}

	public static String generateFormula(AtomContainer m) {

		MolecularFormula mf = (MolecularFormula) MolecularFormulaManipulator.getMolecularFormula(m);
		String formula = MolecularFormulaManipulator.getString(mf);
		return formula;

	}

	/**
	 * Generate html formula from AtomContainer
	 * 
	 * @param m
	 * @return
	 */
	public static String generateHTMLFormula(AtomContainer m) {

		MolecularFormula mf = (MolecularFormula) MolecularFormulaManipulator.getMolecularFormula(m);
		String formula = MolecularFormulaManipulator.getHTML(mf);
		return formula;
	}

	/**
	 * Puts formula in correct order
	 * 
	 * @param formula
	 * @return
	 */
	public static String fixFormula(String formula) {
		MolecularFormula mf = (MolecularFormula) MolecularFormulaManipulator.getMolecularFormula(formula,
				DefaultChemObjectBuilder.getInstance());
		String formulaNew = MolecularFormulaManipulator.getString(mf);
		return formulaNew;
	}

	/**
	 * Generate formula from AtomContainer
	 * 
	 * @param m
	 * @return formula as string
	 */
	public static String calculateFormula(AtomContainer m) {
		MolecularFormula mf = (MolecularFormula) MolecularFormulaManipulator.getMolecularFormula(m);
		String formula = MolecularFormulaManipulator.getString(mf);
		
		return formula;
	}

	/**
	 * Calculate molecular weight from formula
	 * 
	 * @param formula
	 * @return molecular weight
	 */
	public static double calculateMolecularWeight(String formula) {
		MolecularFormula mf = (MolecularFormula) MolecularFormulaManipulator.getMolecularFormula(formula,
				DefaultChemObjectBuilder.getInstance());
		double MW = MolecularFormulaManipulator.getTotalExactMass(mf);
		return MW;
	}

	public static double calculateMolecularWeight(IAtomContainer m) {
		MolecularFormula mf = (MolecularFormula) MolecularFormulaManipulator.getMolecularFormula(m);
		double MW = MolecularFormulaManipulator.getTotalExactMass(mf);
		return MW;
	}

	static void testINCHI() {
		try {
			String smiles = "CC(=O)OC1=CC=CC=C1C(=O)O";// aspirin
			SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			AtomContainer molecule = (AtomContainer) sp.parseSmiles(smiles);
			Inchi inchi = CDKUtilities.generateInChiKey(molecule);
			
			String InChi = inchi.inchi;
			String InChiKey = inchi.inchiKey;

			System.out.println("SMILES:\t" + smiles);
			System.out.println("InChI:\t" + InChi);
			System.out.println("InChIKey:\t" + InChiKey);

			// Output:
			// BSYNRYMUTXBXSQ-WXRBYKJCNA-N
			// InChI=1/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)/f/h11H

			// Pubchem:
			// BSYNRYMUTXBXSQ-UHFFFAOYSA-N
			// InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)

		} catch (InvalidSmilesException ex) {
			ex.printStackTrace();
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
//		String formula="H6C6";
//		System.out.println(CDKUtilities.fixFormula(formula));
		CDKUtilities.testINCHI();

	}
	
	

}
