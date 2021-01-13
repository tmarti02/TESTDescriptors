package gov.epa.TEST.Descriptors.DescriptorFactory;
import gov.epa.TEST.Descriptors.DescriptorCalculations.*;
import gov.epa.TEST.Descriptors.DescriptorUtilities.CDKUtilities;
import gov.epa.TEST.Descriptors.DescriptorUtilities.HueckelAromaticityDetector;
import gov.epa.TEST.Descriptors.DescriptorUtilities.PathFinder;

//java imports:
import java.io.*;
import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import java.lang.reflect.Field;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor;

import org.openscience.cdk.qsar.result.*;
import org.openscience.cdk.qsar.*;
import org.openscience.cdk.qsar.result.DoubleResult;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.config.Isotopes;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.formula.MolecularFormula;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.normalize.Normalizer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.w3c.dom.Document;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.Ring;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

 

public class DescriptorFactory {

	private static final Logger logger = LogManager.getLogger(DescriptorFactory.class);

	/**
	 * @param args
	 */

	// private AdaptMolecule m;
	// FileWriter fw;
	DecimalFormat myF = new DecimalFormat("0.000000");

	public static boolean debug = false;
	public boolean done;

	public String errorMsg = "";

	// words of caution- need to make sure when fortran changes a scalar value
	// in subroutine we
	// make the change in the place where the routine was called from

	// TMM changes:
	// 1. Changed how DV for certain atoms is calculated

	public boolean UnassignedAtoms; // flag for case where there is atoms with
									// an

	// undefined estate fragment

	// IAtomContainer m;

	DescriptorData dd;

	ChiDescriptor cd = new ChiDescriptor();

	EStateFragmentDescriptor esfd = new EStateFragmentDescriptor();

	public MoleculeFragmenter4 moleculefragmenter = new MoleculeFragmenter4();

	TwoDMolecularDescriptors twod = new TwoDMolecularDescriptors();

	TopologicalDescriptors td = new TopologicalDescriptors();

	InformationContentDescriptor icd = new InformationContentDescriptor();

	WalkandPathCountDescriptors wapcd = new WalkandPathCountDescriptors();

	GetawayDescriptors gd = new GetawayDescriptors();

	WHIMDescriptors wd = new WHIMDescriptors();

	TwoDAutoCorrelationDescriptors tdacd = new TwoDAutoCorrelationDescriptors();

	ConstitutionalDescriptors Con = new ConstitutionalDescriptors();

	RandicMolecularProfileDescriptors rmpd = new RandicMolecularProfileDescriptors();

	EState estate = new EState();

	PrincipalMomentOfInertiaDescriptors pmoid = new PrincipalMomentOfInertiaDescriptors();

	MolecularPropertyDescriptors mpd = new MolecularPropertyDescriptors();

	BurdenEigenValueDescriptors bevd = new BurdenEigenValueDescriptors();

	MomentOfInertiaDescriptor mid = new MomentOfInertiaDescriptor();

	public ALOGP alogp = new ALOGP();

	public MLOGP mlogp = new MLOGP();
	// KowWin2 kw=new KowWin2();

	KowDLL kd;

	XLogPDescriptor xlpd = new XLogPDescriptor();

	// needed intermediate values:
	// public IRingSet rs = null;
	// public IRingSet rs2 = null;//ringset with hydrogens in molecule

	IRingSet srs = null;

	double[] D = null;

	double[] DV = null;

	double[] DV2 = null;

	double[] S = null;

	double[] SV = null;

	double[] IS = null; // Intrinsic state

	int[][] Distance = null; // topological distance matrix

	int[] vdd = null; // vertex distance degrees - see pg 113 todeschini book

	double[] EState = null; // E-state for each atom

	double[] HEState = null; // H E-state for each atom

	public String[] Fragment = null; // E-state fragment assigned to each atom

	double[] KHE = null; // kier hall

	LinkedList[] paths; // path list - array of lists

	DecimalFormat df = new DecimalFormat("0.000");

	boolean isStandAlone = false;
	public boolean Calculate3DDescriptors = true;

	public IAtomContainer processedMolecule;

	public DescriptorFactory(boolean isStandAlone) {

		this.isStandAlone = isStandAlone;

		if (isStandAlone) {
			kd = new KowDLL();
		}

	}



	private IAtomContainer Create_molecule_with_Hydrogens(IAtomContainer m) {
		//

		IAtomContainer m2 = null;

		try {
			m2 = (IAtomContainer) m.clone();

			Isotopes isf = Isotopes.getInstance();

			for (int i = 0; i < m.getAtomCount(); i++) {
				IIsotope isotope = isf.getMajorIsotope("H");
				IAtom atom = m2.getAtom(i);

				for (int j = 1; j <= m.getAtom(i).getImplicitHydrogenCount(); j++) {
					Atom hydrogen = new Atom("H");
					isf.configure(hydrogen, isotope);
					m2.addAtom(hydrogen);

					// System.out.println(m2.getAtomNumber(hydrogen));

					Bond newBond = new Bond(atom, hydrogen, IBond.Order.SINGLE);
					m2.addBond(newBond);
				}

			}

			// for (int i=0;i<m2.getAtomCount();i++) {
			// System.out.println(i+"\t"+m2.getAtom(i).getSymbol());
			// }

		} catch (Exception e) {
			logger.catching(e);
		}
		return m2;
	}

	/**
	 * 
	 * @param m
	 * @param dd
	 * @param WriteResultsToFile
	 * @param FindFragments
	 * @param OutputFolder
	 * @return
	 * 
	 * 		TODO: add different error codes for each method so can better
	 *         output error to output file
	 * 
	 */
	public int CalculateDescriptors(IAtomContainer m, DescriptorData dd, boolean FindFragments) {
		// this.m = m;
		this.dd = dd;
		this.errorMsg = "";
		this.done = false;

		IAtomContainer m3d = null;

		// System.out.println("enter calculate descriptors");

		// TODO: add code to catch errors when calculating descriptors
		// int IERRCODE = 0;

		// Normalize before structure image gen:
		this.Normalize(m);

		// System.out.println("here");

		double time1 = System.currentTimeMillis() / 1000.0;


		// clone it:
		try {
			m3d = m.clone();
		} catch (Exception e) {
			logger.catching(e);
			errorMsg = "clone";
		}

		// remove hydrogens from molecule:
		m = (IAtomContainer) AtomContainerManipulator.removeHydrogens(m);// remove
																			// all
																			// hydrogens
																			// first
																			// so
																			// they
																			// all
																			// get
																			// added
																			// to
																			// end

		// IAtomContainer m2=this.Create_molecule_with_Hydrogens(m);
		IAtomContainer m2 = CDKUtilities.addHydrogens(m);

		if (this.done || !errorMsg.equals(""))
			return -1;

		// Dimension all arrays:
		this.DimensionArrays(m);
		if (this.done || !errorMsg.equals(""))
			return -1;

		//TMM-1/10/2018
		for (int i=0;i<m.getAtomCount();i++) {
	        if (m.getAtom(i).getImplicitHydrogenCount()==null) {
	            throw new RuntimeException("Implicit hydrogen count is null");
	        }
		}
		
		// ----------------------------------------------------------------
		// Preliminary Calculations
		// find all the rings (well actually cycles)
		IRingSet rs = FindRings(m);
	
		
		if (this.done || !errorMsg.equals(""))
			return -1;

		IRingSet rs2 = FindRings(m2);

		if (this.done || !errorMsg.equals(""))
			return -1;

		// System.out.print("Finding paths...");

		double timeaa = System.currentTimeMillis() / 1000.0;
		// find all paths in the molecule
		FindPaths(m);
		if (this.done || !errorMsg.equals(""))
			return -1;
		double timebb = System.currentTimeMillis() / 1000.0;
		// System.out.print("done\t"+(timebb-timeaa)+"\r\n");

		// double timecc = System.currentTimeMillis() / 1000.0;
		// System.out.print("Calculating new chi descriptors...");
		// org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor
		// ccd=new
		// org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor();
		// ccd.calculate(m);
		// double timedd = System.currentTimeMillis() / 1000.0;
		// System.out.print("done\t"+(timedd-timecc)+"\r\n");

		// // calculate distance matrix
		PathFinder.CalculateDistanceMatrix(m, Distance);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate vertex distance degrees
		PathFinder.CalculateVertexDistanceDegrees(m, Distance, vdd);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // detect aromaticity of rings and atoms:
		DetectAromaticityOld(m, rs);
		DetectAromaticityOld(m2, rs2);

		// DetectAromaticity(m,rs);
		// DetectAromaticity(m2,rs2);

		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate Estates
		CalculateEStates(m);
		if (this.done || !errorMsg.equals(""))
			return -1;
		// // ----------------------------------------------------------------
		// // begin descriptor calcs
		//
		// // calculate kier hall chi descriptors:
		CalculateChiDescriptors(m, rs);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//

		// // assign estate fragments to each atom:
		CalculateEStateFragments(m, rs);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate Kappa and kappa alpha descriptors
		CalculateKappaDescriptors(m);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // determine our fragments:
		if (FindFragments)
			CalculateFragments(m, m2, rs, rs2);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate 2 dimensional descriptors:
		Calculate2DDescriptors(m);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate topological descriptors:
		CalculateTopologicalDescriptors(m);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate information content descriptors:
		CalculateInformationContentDescriptors(m, rs);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate constitutional descriptors:
		CalculateConstitutionalDescriptors(m, rs);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate walk and path count descriptors:
		CalculateWalkAndPathCountDescriptors(m, rs);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // ADAPT's Molecular Distance Edge Descriptors();
		CalculateMolecularDistanceEdgeDescriptors(m);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate 2d autocorrelation descriptors:
		Calculate2DAutocorrelationDescriptors(m);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		//
		// // calculate Burden eigenvalue descriptors:
		CalculateBurdenEigenvalueDescriptors(m2, rs2);
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		// // calculate ALOGP:
		alogp.Calculate(m, dd, Fragment, rs);
		if (this.done || !errorMsg.equals(""))
			return -1;

		// // calculate MLOGP:
		// mlogp.Calculate(m,dd,Fragment,rs);
		// if (this.done || !errorMsg.equals("")) return -1;
		//
		//
		// if (this.done) return -1;
		//
		// calculate molecular property descriptors:
		mpd.Calculate(m, dd, Fragment);
		if (this.done || !errorMsg.equals(""))
			return -1;

		
//		if (isStandAlone) {
//			this.CalculateKOWWIN_From_DLL(m);
//		}
		
		
		if (this.done || !errorMsg.equals(""))
			return -1;
		//
		//
		this.CalculateXLogP(m2, rs2);// new code adds explicit hs

		return 0;

	}
	
	

	private boolean HasHydrogensInCa(IAtomContainer mol) {

		for (int i = 0; i < mol.getAtomCount(); i++) {
			List ca = mol.getConnectedAtomsList(mol.getAtom(i));
			for (int j = 0; j < ca.size(); j++) {
				IAtom caj = (IAtom) ca.get(j);
				if (caj.getSymbol().equals("H")) {
					System.out.println("--Found hydrogens in ca");
					return true;
				}
			}

		}

		return false;

	}

	public void Normalize(IAtomContainer m) {
		try {
			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = dbf.newDocumentBuilder();

			String DataFile = "fixsmiles.xml";
			InputStream ins = this.getClass().getClassLoader().getResourceAsStream(DataFile);

			Document d = db.parse(ins);

			Normalizer.normalize(m, d);
		} catch (Exception e) {
			logger.catching(e);
			this.errorMsg = "normalize";
		}
	}

	private void CalculateMomentOfInertia(IAtomContainer m3d, DescriptorData dd) {
		try {

			DescriptorValue dv = mid.calculate(m3d); //

			DoubleArrayResult dar = (DoubleArrayResult) dv.getValue();

			for (int i = 0; i <= dar.length() - 1; i++) {
				Field myField = dd.getClass().getField("MOMI" + (i + 1));
				myField.setDouble(dd, dar.get(i));
				// System.out.println("MOMI"+(i+1)+"\t"+dar.get(i));
			}

			// mid.calculate(m3d, dd);
		} catch (Exception e) {
			logger.catching(e);
		}
	}

	private void CalculateXLogP(IAtomContainer m2, IRingSet rs2) {
		try {
			// convert nitro groups to uncharged versions for XLOGP to work
			// right:
			// this.FixNitroGroups2(m2);//TODO- doesnt seem affect new xlogp in
			// CDK1.5
			xlpd.setCheckAromaticity(false);
			xlpd.setSalicylFlag(false);

			// DescriptorValue dv = this.xlpd.calculate(m2);
			DescriptorValue dv = this.xlpd.calculate(m2, rs2);
			DoubleResult dr = (DoubleResult) dv.getValue();
			dd.XLOGP = dr.doubleValue();
			dd.XLOGP2 = dd.XLOGP * dd.XLOGP;

		} catch (Exception e) {
			logger.catching(e);
		}
	}

	private void CalculateConstitutionalDescriptors(IAtomContainer m, IRingSet rs) {
		if (debug)
			System.out.print("Calculating Constitutional Descriptors...");
		double time3 = System.currentTimeMillis() / 1000.0;

		Con.Calculate(m, dd, EState, rs, srs);

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	private void CalculateBurdenEigenvalueDescriptors(IAtomContainer m2, IRingSet rs2) {
		if (debug)
			System.out.print("Calculating BurdenEigenvalueDescriptors...");
		double time3 = System.currentTimeMillis() / 1000.0;
		bevd.Calculate(m2, rs2, dd);
		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");

		// System.out.println("done in DF = "+this.done);
	}

//	private void CalculateKOWWIN_From_DLL(IAtomContainer m) {
//		SmilesGenerator sg = new SmilesGenerator();
//		Lookup lookup = new Lookup();
//
//		if (debug)
//			System.out.print("Calculating KOWWIN from DLL...");
//		double time3 = System.currentTimeMillis() / 1000.0;
//
//		String smiles = "";
//
//		if (dd.ID instanceof String) {
//			// first try to look up the smiles from the database:
//			smiles = lookup.LookUpSmiles(dd.ID);
//
//			// System.out.println(smiles);
//		}
//
//		// if dont have in the database generate the smiles using CDK:
//
//		// smiles="missing";
//
//		if (smiles.equals("missing")) {
//			smiles = sg.createSMILES(m);
//			// System.out.println("missing so created for "+dd.CAS);
//		}
//		// System.out.println(dd.CAS+"\t"+smiles);
//
//		dd.KLOGP = kd.Calculate(smiles);
//		dd.KLOGP2 = dd.KLOGP * dd.KLOGP;
//
//		// System.out.println("dd.KLOGP="+dd.KLOGP);
//		// System.out.println("LogPexp="+kd.KowExp);
//
//		double time4 = System.currentTimeMillis() / 1000.0;
//
//		if (debug)
//			System.out.println("done-" + df.format(time4 - time3) + " secs");
//
//	}

	private void CalculateKOWWIN_From_DLL(String smiles) {

		if (debug)
			System.out.print("Calculating KOWWIN from DLL...");
		double time3 = System.currentTimeMillis() / 1000.0;

		dd.KLOGP = kd.Calculate(smiles);
		dd.KLOGP2 = dd.KLOGP * dd.KLOGP;

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");

	}

	/**
	 * @param m
	 * @param dd
	 */
	private void Calculate2DAutocorrelationDescriptors(IAtomContainer m) {
		if (debug)
			System.out.print("Calculating 2DAutocorrelationDescriptors...");
		double time3 = System.currentTimeMillis() / 1000.0;
		tdacd.Calculate(m, dd, Distance);
		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");

	}

	private void CalculateFragments(IAtomContainer m, IAtomContainer m2, IRingSet rs, IRingSet rs2) {

		if (debug)
			System.out.print("Calculating fragments...");

		double time3 = System.currentTimeMillis() / 1000.0;

		// determine molecular fragments:

		dd.FragmentList = moleculefragmenter.Calculate(m, m2, rs, rs2, dd, Fragment);
		// dd.FragmentList = moleculefragmenter.fragmentMolecule(m);

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");

		/*
		 * String fragnames=(String)dd.FragmentList.get(0); List
		 * namelist=Utilities.Parse(fragnames,"|");
		 * 
		 * for (int i=0;i<=namelist.size()-1;i++) { Integer
		 * fragval=(Integer)dd.FragmentList.get(i+1); int
		 * fragval2=fragval.intValue(); if (fragval2>0) {
		 * System.out.println(namelist.get(i)+"\t"+fragval2); } }
		 */

	}

	// private void CalculateKappa0() {
	// if (debug)
	// System.out.print("Calculating Kappa0...");
	// double time3 = System.currentTimeMillis() / 1000.0;
	//
	// // calculate Kappa0 using SV topological states of each atom:
	// Kappa0Descriptor2.Calculate(m, SV, dd);
	//
	// double time4 = System.currentTimeMillis() / 1000.0;
	//
	// if (debug)
	// System.out.println("done-" + df.format(time4 - time3) + " secs");
	// }

	/**
	 * @param m
	 * @param dd
	 */
	private void CalculateMolecularDistanceEdgeDescriptors(IAtomContainer m) {
		// -----------------------------------------------------------------
		// calculate Molecular distance edge descriptors:
		// MolecularDistanceEdgeDescriptors.Calculate(m,dd,Fragment,Distance);

		if (debug)
			System.out.print("Calculating MolecularDistanceEdgeDescriptors...");
		double time3 = System.currentTimeMillis() / 1000.0;

		MolecularDistanceEdgeDescriptors.Calculate2(m, dd, Distance);
		// MolecularDistanceEdgeDescriptors.Calculate(m, dd,Fragment,Distance);

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	private void CalculateWalkAndPathCountDescriptors(IAtomContainer m, IRingSet rs) {
		if (debug)
			System.out.print("Calculating walk and path count descriptors...");
		double time3 = System.currentTimeMillis() / 1000.0;
		wapcd.Calculate(m, dd, paths, D, vdd, rs);
		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	private void CalculateInformationContentDescriptors(IAtomContainer m, IRingSet rs) {
		if (debug)
			System.out.print("Calculating information content descriptors...");

		double time3 = System.currentTimeMillis() / 1000.0;

		icd.Calculate(m, dd, Distance, rs, paths, D, DV, S, SV, vdd);

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	private void CalculateTopologicalDescriptors(IAtomContainer m) {

		if (debug)
			System.out.print("Calculating topological descriptors...");
		double time3 = System.currentTimeMillis() / 1000.0;

		// calculate topological descriptors:
		td.Calculate(m, dd, D, DV, vdd, Distance, EState, IS);
		double time4 = System.currentTimeMillis() / 1000.0;
		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");

	}

	private void Calculate2DDescriptors(IAtomContainer m) {
		if (debug)
			System.out.print("Calculating 2D descriptors...");

		double time3 = System.currentTimeMillis() / 1000.0;

		twod.calculate(m, dd, IS, EState, HEState, Fragment);

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	private void CalculateKappaDescriptors(IAtomContainer m) {
		if (debug)
			System.out.print("Calculating kappa descriptors...");

		double time3 = System.currentTimeMillis() / 1000.0;

		// calculate Kappa descriptors:
		KappaDescriptor.Calculate(m, paths, dd);

		// calculate Kappa alpha descripors:
		KappaAlphaDescriptor.Calculate(m, paths, dd);

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	/**
	 * @param m
	 * @param dd
	 */
	private void CalculateEStateFragments(IAtomContainer m, IRingSet rs) {

		if (debug)
			System.out.print("Calculating fragment descriptors...");

		double time3 = System.currentTimeMillis() / 1000.0;

		// calculates fragment estates and assigns fragments to each atom
		esfd.Calculate(m, dd, EState, HEState, rs, Fragment);

		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			if (!(Fragment[i] instanceof String)) {
				System.out.print("Atom " + (i + 1) + " has no fragment assigned\n");
				UnassignedAtoms = true;
			}
		}

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	/**
	 * @param m
	 * @param dd
	 */
	private void CalculateChiDescriptors(IAtomContainer m, IRingSet rs) {
		if (debug)
			System.out.print("Calculating chi descriptors...");

		double time3 = System.currentTimeMillis() / 1000.0;

		// calculate chi descriptors:
		cd.Calculate(m, dd, D, DV, paths, rs);

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	// private void CalculateKowWin() {
	// if (debug)
	// System.out.print("Calculating Kowwin...");
	//
	// double time3 = System.currentTimeMillis() / 1000.0;
	//
	// this.kw.Calculate(m,dd,Fragment,rs);
	//
	// double time4 = System.currentTimeMillis() / 1000.0;
	//
	// if (debug)
	// System.out.println("done-" + df.format(time4 - time3) + " secs");
	// }

	/**
	 * @param m
	 */
	private void CalculateEStates(IAtomContainer m) {

		UnassignedAtoms = false;

		if (debug)
			System.out.print("Calculating estates...");

		double time3 = System.currentTimeMillis() / 1000.0;

		// calculate D,DV,Distance matrix, Intrinsic states, Estates:
		estate.Calculate(m, D, DV, DV2, KHE, Distance, IS, EState, HEState);

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");
	}

	private void DimensionArrays(IAtomContainer m) {
		D = new double[m.getAtomCount()];
		DV = new double[m.getAtomCount()];
		DV2 = new double[m.getAtomCount()];
		S = new double[m.getAtomCount()];
		SV = new double[m.getAtomCount()];
		IS = new double[m.getAtomCount()];
		Distance = new int[m.getAtomCount()][m.getAtomCount()];
		vdd = new int[m.getAtomCount()];
		EState = new double[m.getAtomCount()];
		HEState = new double[m.getAtomCount()];
		Fragment = new String[m.getAtomCount()];
		KHE = new double[m.getAtomCount()];
		dd.FragmentList = null;
	}

	/**
	 * @param m
	 * @param rs
	 */
	public void DetectAromaticityOld(IAtomContainer m, IRingSet rs) {
		double time3 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.print("Finding aromaticity...");

		try {

			// for (int i=0;i<=m.getAtomCount()-1;i++) {
			// m.getAtom(i).setFlag(CDKConstants.ISAROMATIC,false);
			// }
			//
			// for (int i=0;i<=rs.getAtomContainerCount()-1;i++) {
			// Ring r=(Ring)rs.getAtomContainer(i);
			// r.setFlag(CDKConstants.ISAROMATIC,false);
			// }

			// do it multiple times to make sure it gets it right:
			for (int i = 0; i <= rs.getAtomContainerCount() - 1; i++) {
				HueckelAromaticityDetector.detectAromaticity(m, rs, false);
			}

			// for (int i=0;i<=rs.getAtomContainerCount()-1;i++) {
			// Ring r=(Ring)rs.getAtomContainer(i);
			// System.out.println("Ring size = "+r.getAtomCount()+", aromaticity
			// = "+r.getFlag(CDKConstants.ISAROMATIC));
			// }

			// for (int i=0;i<=m.getAtomCount()-1;i++) {
			// System.out.println(m.getAtom(i).getFlag(CDKConstants.ISAROMATIC));
			// }

			// figure out which atoms are in aromatic rings:
			// HueckelAromaticityDetector.detectAromaticity(m,true);
			// figure out which rings are aromatic:
			// HueckelAromaticityDetector.detectAromaticity(m, rs,true);

		} catch (Exception e) {
			logger.catching(e);
			errorMsg = "DetectAromaticity";
		}

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");

	}

	/**
	 * @param m
	 * @param rs
	 * 
	 */

	// just in case set all flags to false:
	// for (int i=0;i<m.getAtomCount();i++) {
	// m.getAtom(i).setFlag(CDKConstants.ISAROMATIC,false);
	// }

	// for (int i=0;i<rs.getAtomContainerCount();i++) {
	// Ring r=(Ring)rs.getAtomContainer(i);
	// r.setFlag(CDKConstants.ISAROMATIC,false);
	// }

	public void DetectAromaticity(IAtomContainer m, IRingSet rs) {

		/*
		 * http://cdk.github.io/cdk/2.0/docs/api/org/openscience/cdk/aromaticity
		 * /Aromaticity.html
		 * 
		 * A configurable model to perceive aromatic systems. Aromaticity is
		 * useful as both a chemical property indicating stronger stabilisation
		 * and as a way to treat different resonance forms as equivalent. Each
		 * has its own implications the first in physicochemical attributes and
		 * the second in similarity, depiction and storage. To address the
		 * resonance forms, several simplified (sometimes conflicting) models
		 * have arisen. Generally the models loosely follow Hückel's rule for
		 * determining aromaticity. A common omission being that planarity is
		 * not tested and chemical compounds which are non-planar can be
		 * perceived as aromatic. An example of one such compound is,
		 * cyclodeca-1,3,5,7,9-pentaene. Although there is not a single
		 * universally accepted model there are models which may better suited
		 * for a specific use (Cheminformatics Toolkits: A Personal Perspective,
		 * Roger Sayle). The different models are often ill-defined or
		 * unpublished but it is important to acknowledge that there are
		 * differences (see. Aromaticity Perception Differences, Blue Obelisk).
		 * Although models may get more complicated (e.g. considering tautomers)
		 * normally the reasons for differences are:
		 * 
		 * the atoms allowed and how many electrons each contributes the
		 * rings/cycles are tested
		 * 
		 * This implementation allows configuration of these via an
		 * ElectronDonation model and CycleFinder. To obtain an instance of the
		 * electron donation model use one of the factory methods,
		 * ElectronDonation.cdk(), ElectronDonation.cdkAllowingExocyclic(),
		 * ElectronDonation.daylight() or ElectronDonation.piBonds().
		 * 
		 * TODO- see if other ElectronDonation models work better (i.e. piBonds)
		 * 
		 * 
		 * // mimics the old CDKHuckelAromaticityDetector which uses the CDK
		 * atom types ElectronDonation model = ElectronDonation.cdk();
		 * CycleFinder cycles = Cycles.cdkAromaticSet(); Aromaticity aromaticity
		 * = new Aromaticity(model, cycles);
		 * 
		 * // apply our configured model to each molecule, the CDK model //
		 * requires that atom types are perceived for (IAtomContainer molecule :
		 * molecules) {
		 * AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(
		 * molecule); aromaticity.apply(molecule); }
		 * 
		 */

		double time3 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.print("Finding aromaticity...");

		try {

			// ********************************************************************************************

			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(m);

			ElectronDonation model = ElectronDonation.cdk();
			// ElectronDonation model = ElectronDonation.piBonds();

			CycleFinder cycles = Cycles.cdkAromaticSet();

			// Aromaticity aromaticity = new Aromaticity(model, cycles);
			// Aromaticity aromaticity=new Aromaticity(ElectronDonation.cdk(),
			// Cycles.or(Cycles.all(), Cycles.relevant()));

			// aromaticity.apply(m);

			// For now use the legacy code:
			Aromaticity.cdkLegacy().apply(m);// TODO compare legacy method to
												// new method

			// ********************************************************************************************
			// if all atoms are aromatic, set aromaticity of ring to true:
			for (int i = 0; i < rs.getAtomContainerCount(); i++) {
				Ring r = (Ring) rs.getAtomContainer(i);

				boolean allAromatic = true;

				for (int j = 0; j < r.getAtomCount(); j++) {
					Atom atomj = (Atom) r.getAtom(j);
					if (!atomj.isAromatic()) {
						allAromatic = false;
						break;
					}
				}
				r.setFlag(CDKConstants.ISAROMATIC, allAromatic);

				// aromaticity.apply(r);
				// Aromaticity.cdkLegacy().apply(r);
				// System.out.println("Ring size = "+r.getAtomCount()+",
				// aromaticity =
				// "+r.getFlag(CDKConstants.ISAROMATIC)+"\t"+allAromatic);
			}

			// ********************************************************************************************
			// set bonds to be aromatic if both atoms are aromatic and in same
			// ring:

			for (int i = 0; i < m.getBondCount(); i++) {
				IBond bondi = m.getBond(i);

				IAtom atom0 = bondi.getAtom(0);
				IAtom atom1 = bondi.getAtom(1);

				if (atom0.getFlag(CDKConstants.ISAROMATIC) && atom1.getFlag(CDKConstants.ISAROMATIC)) {
					// check if in same ring:
					boolean SameRing = EStateFragmentDescriptor.InSameAromaticRing(m, atom0, atom1, rs);

					if (SameRing) {
						// System.out.println(i+"\t"+bondi.getFlag(CDKConstants.ISAROMATIC));
						bondi.setFlag(CDKConstants.ISAROMATIC, true);
					}
				}
			}

			// int count=0;
			// for (int i=0;i<m.getAtomCount();i++) {
			// if (m.getAtom(i).getFlag(CDKConstants.ISAROMATIC) &&
			// m.getAtom(i).getSymbol().equals("C")) {
			// count++;
			// }
			// System.out.println(m.getAtom(i).getSymbol()+"\t"+m.getAtom(i).getFlag(CDKConstants.ISAROMATIC));
			// }
			// System.out.println("aromaticCarbonCount="+count);

		} catch (Exception e) {
			logger.catching(e);
			errorMsg = "DetectAromaticity";
		}

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");

	}

	/**
	 * @param m
	 * @param rs
	 * @return
	 */
	private IRingSet FindRings(IAtomContainer m) {

		IRingSet ringSet = null;

		double time3 = System.currentTimeMillis() / 1000.0;
		if (debug)
			System.out.print("Finding rings...");

		try {
			AllRingsFinder arf = new AllRingsFinder();

			arf.usingThreshold(AllRingsFinder.Threshold.PubChem_994);

			// arf.setTimeout(100000);
			// TODO: set the threshold instead of timeout
			ringSet = arf.findAllRings(m);

			// ringSet = Cycles.sssr(m).toRingSet();//TODO how does compare to
			// method above???
		} catch (Exception e) {
			this.errorMsg = "FindRings";
			logger.catching(e);
		}
		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.println("done-" + df.format(time4 - time3) + " secs");

		return ringSet;

	}

	private void FindPaths(IAtomContainer m) {
		// find all paths up to length 10:

		// TODO add code to figure out if it will run out of memory

		double time3 = System.currentTimeMillis() / 1000.0;

		if (debug)
			System.out.print("Finding paths...");

		// paths = PathFinder.FindPaths3(10, m);
		paths = PathFinder.FindPaths4(10, m);

		// System.out.println(paths.length);

		// if (paths==null) done=true;

		if (paths == null) {
			errorMsg = "FindPaths";
		}

		double time4 = System.currentTimeMillis() / 1000.0;

		if (debug) {
			System.out.println("done-" + df.format(time4 - time3) + " secs");
		}

	}

	
	

	private String GetFragmentNames(String Delimiter, DescriptorData dd) {
		String Line = "";

		try {
			int count = 0;

			for (int i = 0; i < dd.strFragments.length; i++) {

				String fragname = dd.strFragments[i];

				if (Delimiter.equals(",") && fragname.indexOf(",") > -1) {
					Line += "\"" + fragname + "\"";
				} else {
					Line += fragname;
				}

				if (i < dd.strFragments.length - 1)
					Line += Delimiter;

				count++;
			}
		} catch (Exception e) {
			logger.catching(e);
		}
		return Line;
	}



	private String GetDescriptorNames(DescriptorData dd, String[] varlist, String Delimiter) {

		String Line = "";
		try {

			int count = 0;

			for (int i = 0; i <= varlist.length - 1; i++) {
				// System.out.println(dd.varlist2d[i]);

				Field myField = dd.getClass().getField(varlist[i]);
				String[] names = (String[]) myField.get(dd);

				for (int j = 0; j <= names.length - 1; j++) {
					count++;
					Field myField2 = dd.getClass().getField(names[j]);

					String sname = names[j];

					if (Delimiter.equals(",") && sname.indexOf(",") > -1) {
						sname = "\"" + sname + "\"";
					}

					// System.out.println(i+"\t"+j);
					if (i != varlist.length - 1 || j != names.length - 1) {
						Line += sname + Delimiter;
					} else {
						Line += sname + "";
					}
				}

			}

			// System.out.println(count);

		} catch (Exception e) {
			logger.catching(e);
		}
		return Line;
	}

	private String GetFragmentValues(String Delimiter, DescriptorData dd) {
		String Line = "";

		DecimalFormat df = new DecimalFormat("0");

		try {
			int count = 0;
			for (java.util.Enumeration e = dd.FragmentList.keys(); e.hasMoreElements();) {
				String strVar = (String) e.nextElement();
				double Val = (Double) dd.FragmentList.get(strVar);

				Line += df.format(Val);

				if (e.hasMoreElements())
					Line += (Delimiter);

				count++;
			}
		} catch (Exception e) {
			logger.catching(e);
		}

		return Line;
	}

	//	/**
//	 * Writes out descriptors to readable file for looking at all descriptors in
//	 * a single vertical table
//	 * 
//	 * @param m
//	 * @param dd
//	 */
//	public void WriteOut2DDescriptors(IAtomContainer m, DescriptorData dd, String folderpath, String Delimiter) {
//
//		try {
//			FileWriter fw = new FileWriter(folderpath + "/2d_descriptors.txt");
//
//			DecimalFormat myD8 = new DecimalFormat("0.########");
//
//			String[] varlist = dd.varlist2d;
//
//			fw.write("Descriptor" + Delimiter + "Value\r\n");
//			for (int i = 0; i <= varlist.length - 1; i++) {
//				// System.out.println(dd.varlist2d[i]);
//
//				Field myField = dd.getClass().getField(varlist[i]);
//				String[] names = (String[]) myField.get(dd);
//				for (int j = 0; j <= names.length - 1; j++) {
//					Field myField2 = dd.getClass().getField(names[j]);
//
//					String val = myD8.format(myField2.getDouble(dd));
//
//					fw.write(names[j] + Delimiter + val + "\r\n");
//				}
//			}
//
//			// for (java.util.Enumeration e = dd.FragmentList.keys(); e
//			// .hasMoreElements();) {
//			// String strVar = (String) e.nextElement();
//			// double Val = (Double) dd.FragmentList.get(strVar);
//			//
//			// fw.write(strVar+Delimiter+Val + "\r\n");
//			// }
//
//			for (int i = 0; i < dd.strFragments.length; i++) {
//				String strVar = dd.strFragments[i];
//				double Val = (Double) dd.FragmentList.get(strVar);
//				fw.write(strVar + Delimiter + Val + "\r\n");
//			}
//
//			fw.close();
//
//		} catch (Exception e) {
//			logger.catching(e);
//		}
//	}

	
	
	
	
	

	

	private String GetDescriptorValues(DescriptorData dd, String[] varlist, String Delimiter) {

		String Line = "";
		DecimalFormat myD8 = new DecimalFormat("0.########");

		try {

			int countwritten = 0;

			for (int i = 0; i <= varlist.length - 1; i++) {
				// System.out.println(dd.varlist2d[i]);

				Field myField = dd.getClass().getField(varlist[i]);
				String[] names = (String[]) myField.get(dd);
				for (int j = 0; j <= names.length - 1; j++) {
					Field myField2 = dd.getClass().getField(names[j]);

					String val = myD8.format(myField2.getDouble(dd));

					countwritten++;
					// System.out.println(i+"\t"+j);
					if (i != varlist.length - 1 || j != names.length - 1) {
						Line += (val + Delimiter);
					} else {
						Line += (val + "");
					}
				}
			}
			// System.out.println(countwritten);

		} catch (Exception e) {
			logger.catching(e);
		}

		return Line;
	}
}

	