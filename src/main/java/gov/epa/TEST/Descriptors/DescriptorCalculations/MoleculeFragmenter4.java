package gov.epa.TEST.Descriptors.DescriptorCalculations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
//import org.openscience.cdk.Molecule;
import org.openscience.cdk.Ring;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.AtomicProperties;

//import ToxPredictor.Utilities.CDKUtilities;
//import ToxPredictor.Utilities.SmilesParser;

/**
 * 
 * This class finds counts of molecular fragments in a molecule
 * 
 * 
 * TODO: Unaccounted for fragments: 
 * AN+: 113359-04-9 (need SaaaNp Estate frag)
 * As=S: 2533-82-6
 * I=O [aromatic attach]: 304-91-6, 536-80-1
 * [C-]#[Si-]: 409-21-2
 * S(=O)(=S): 41206-16-0
 * >S=N-: 51528-03-1: 
 * [C-]#[O+]: 630-08-0
 * I(=O)=O [aromatic attach]: 696-33-3
 * Carbon: 7440-44-0
 * CH4: 74-82-8
 * HC#N: 74-90-8
 * 
 * @author TMARTI02
 * 
 */

public class MoleculeFragmenter4 {

	/**
	 * @param args
	 */

	DescriptorData dd;

//	IRingSet rs;

	IRingSet rs2; // ringset for molecule w/ hydrogens
	
	//IAtomContainer m;
	IAtomContainer m2; //molecule w/ hydrogens

	String[] Fragment;

	public boolean[] Assigned;

	public String[] AssignedFragment;

	Hashtable htFragVal;

//	SmilesParser sp = new SmilesParser();

	Isotopes isf;
	
	boolean Debug=false;
	
	AtomicProperties ap;
	
//	public double MWcalc=0;
//	public double MW=0;
	
	public boolean CalculateReactiveFlags=true;
	
	
	UniversalIsomorphismTester uit=new UniversalIsomorphismTester();

	public MoleculeFragmenter4() {

		try {
//			isf = IsotopeFactory.getInstance(new Molecule().getBuilder()); //Needed to get isotope information for SMILES elements
			
//			isf=new IsotopeFactory();
			
			isf=Isotopes.getInstance();
//			
			
			ap=AtomicProperties.getInstance();
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	/**
	 * adds hydrogens to m to create m2
	 * This method assumes that the hydrogen counts have already been set by CalcFragmentEStates method  
	 * @param m
	 * @return
	 */
	private IAtomContainer Create_m2(IAtomContainer m) {
//  
		
		m2 = null;

		try {
			m2 = (IAtomContainer) m.clone();

			for (int i = 0; i < m.getAtomCount(); i++) {
				IIsotope isotope = isf.getMajorIsotope("H");
				IAtom atom = m2.getAtom(i);

				for (int j = 1; j <= m.getAtom(i).getImplicitHydrogenCount(); j++) {
					Atom hydrogen = new Atom("H");
					isf.configure(hydrogen, isotope);
					m2.addAtom(hydrogen);
					Bond newBond = new Bond(atom, hydrogen, IBond.Order.SINGLE);
					m2.addBond(newBond);
				}

			}

			for (int i = 0; i < m2.getAtomCount(); i++) {
				//System.out.println(m2.getAtom(i).getSymbol()+"\t"+m2.getAtom(i).getHydrogenCount());				
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
		return m2;
	}

	
	private void FindFrags() {
		
		//TODO- add generic [N-] fragment
		
		this.getCount("->S<-","S","Any","Any","S","SssssssS","","Any");

		this.getCount("[C-]#[N+]","[C-]#[N+]","Any","Any","N","StsNp","","Any");
		
		this.getCount("S=C=S","S=C=S","Any","Any","Any","Any","","Any");
		this.getCount("HN=C=N","[H]N=C=N","Any","Any","Any","Any","","Any");
		this.getCount("N=C=N","N=C=N","Any","Any","Any","Any","","Any");
		
		
		this.getCount("C=C=O","C=C=O","Any","Any","C","Any","","Any");
		this.getCount("C=[N+][O-]","C=[N+][O-]","Any","Any","N","SdssNp","","Any");
		this.getCount("C=[N+]=[N-]","C=[N+]=[N-]","Any","Any","C","Any","","Any");


		this.getCount("-N=S=O","N=S=O","Any","Any","N","Any","","Any");
		
//		added by TMM-11/14/07:
		this.getCount("-N=[N+][O-]","N=[N+][O-]","Any","Any","N","SdssNp","","Any");
		
		this.getCount("-N=[N+]=[N-]","N=[N+]=[N-]","Any","Any","N","Any","","Any");
		this.getCount("-[N+]#N","[N+]#N","Any","Any","N","StsNp","","Any");
		this.getCount("A[N+][O-]","[N+][O-]","2Aromatic","Any","N","SsaaNp","","Any");
		this.getCount("N+ [four single bonds]","[N+]","Any","Any","N","SssssNp","","Any");

		this.getCount("CH2=C(CH3)C(=O)O-","C([H])([H])=C(C([H])([H])([H]))C(=O)O","Any","Any","O","SssO","","Any");
		this.getCount("CH2=CHC(=O)O-","C([H])([H])=C([H])C(=O)O","Any","Any","O","SssO","","Any");
		this.getCount("-N=C=S","N=C=S","Any","Any","N","Any","","Any");

		this.getCount("-SC(=S)S-","SC(=S)S","Any","Any","S","Any","","Any");
		this.getCount("-SC(=S)O-","SC(=S)O","Any","Any","S","Any","","Any");
		this.getCount("-SC(=O)S-","SC(=O)S","Any","Any","S","Any","","Any");
		this.getCount("-SC(=O)O-","SC(=O)O","Any","Any","S","Any","","Any");
		this.getCount("-OC(=S)O-","OC(=S)O","Any","Any","O","Any","","Any");

		this.getCount("-C(=O)OC(=O)-","C(=O)OC(=O)","Any","Any","C","Any","","Any");

		this.getCount("-NH- [3 membered ring]","C1CN([H])1","Any","Any","N","SssNH","C","Any");
		this.getCount("-N< [3 membered ring]","C1CN1","Any","Any","N","SsssN","C","Any");

		this.getCount("-NH2 [attached to P]","N([H])([H])","Any","Any","N","SsNH2","","P");
		this.getCount("-NH- [attached to P]","N([H])","Any","Any","N","SssNH","","XP");
		this.getCount("-N< [attached to P]","N","Any","Any","N","SsssN","","XP");

		this.getCount("Sn=O","[Sn]=O","Any","Any","Sn","Any","","Any");
		this.getCount("Sn [oxygen attach]","[Sn]","Any","Any","Sn","Any","","XO");
		this.getCount("Sn [aromatic attach]","[Sn]","Aromatic","Any","Sn","Any","","Any");
		this.getCount("Sn","[Sn]","Any","Any","Sn","Any","","Any");

		this.getCount("Si [aromatic attach]","[Si]","Aromatic","Any","Si","Any","","Any");
		this.getCount("Si [oxygen, aromatic attach]","[Si]","Aromatic","Any","Si","Any","","XO");
		this.getCount("Si [oxygen attach]","[Si]","Any","Any","Si","Any","","XO");
		this.getCount("Si","[Si]","Any","Any","Si","Any","","Any");

		this.getCount("-S(=O)(=O)- [2 nitrogen attach]","S(=O)(=O)","Any","Any","S","SddssS","","2N");
		this.getCount("-S(=O)(=O)- [nitrogen, aromatic attach]","S(=O)(=O)","Aromatic","Any","S","SddssS","","N");
		this.getCount("-S(=O)(=O)- [aromatic attach]","S(=O)(=O)","Aromatic","Any","S","SddssS","","Any");
		this.getCount("-S(=O)(=O)- [olefinic attach]","S(=O)(=O)","Olefinic","Any","S","SddssS","","Any");
		this.getCount("-S(=O)(=O)- [nitrogen, aliphatic attach]","S(=O)(=O)","Any","Any","S","SddssS","","N");
		this.getCount("-S(=O)(=O)- [aliphatic attach]","S(=O)(=O)","Any","Any","S","SddssS","","Any");

		this.getCount("-S(=O)- [2 nitrogen attach]","S(=O)","Any","Any","S","SdssS","","2N");
		this.getCount("-S(=O)- [nitrogen, aromatic attach]","S(=O)","Aromatic","Any","S","SdssS","","N");
		this.getCount("-S(=O)- [aromatic attach]","S(=O)","Aromatic","Any","S","SdssS","","Any");
		this.getCount("-S(=O)- olefinic attach]","S(=O)","Olefinic","Any","S","SdssS","","Any");
		this.getCount("-S(=O)- [nitrogen, aliphatic attach]","S(=O)","Any","Any","S","SdssS","","N");
		this.getCount("-S(=O)- [aliphatic attach]","S(=O)","Any","Any","S","SdssS","","Any");

		this.getCount("-OH [aromatic attach]","O([H])","Aromatic","Any","O","SsOH","","Any");
		this.getCount("-SH [aromatic attach]","S([H])","Aromatic","Any","S","SsSH","","Any");
		this.getCount("-N=O [aromatic attach]","N=O","Aromatic","Any","N","SdsN","","Any");
		this.getCount("-ONO2","O[N+](=O)[O-]","Any","Any","O","Any","","Any");
		this.getCount("-NO2 [aromatic attach]","[N+](=O)[O-]","Aromatic","Any","N","SdssNp","","Any");
		this.getCount("-NO2 [olefinic attach]","[N+](=O)[O-]","Olefinic","Any","N","SdssNp","","Any");
		this.getCount("-NO2 [nitrogen attach]","[N+](=O)[O-]","Any","Any","N","SdssNp","","N");
		this.getCount("-NO2 [aliphatic attach]","[N+](=O)[O-]","Any","Any","N","SdssNp","","Any");
		this.getCount("-N=C=O [aromatic attach]","N=C=O","Aromatic","Any","N","SdsN","","Any");
		this.getCount("-N=C=O [aliphatic attach]","N=C=O","Any","Any","N","Any","","Any");
		this.getCount("-COOH [aromatic attach]","C(=O)O([H])","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-OC(=O)O-","OC(=O)O","Any","Any","O","Any","","Any");

//		need for IGC50 chemicals like 137-40-6:
		this.getCount("-C(=O)[O-]","C(=O)[O-]","Any","Any","O","SsOm","","Any");
		
		this.getCount("-C(=O)O- [aromatic attach]","C(=O)O","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=O)SH [aromatic attach]","C(=O)S[H]","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=O)S- [aromatic attach]","C(=O)S","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=S)SH [aromatic attach]","C(=S)S[H]","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=S)S- [aromatic attach]","C(=S)S","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=S)OH [aromatic attach]","C(=S)O[H]","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=S)O- [aromatic attach]","C(=S)O","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=S)- [nitrogen, aromatic attach]","C(=S)","Aromatic","Any","C","SdssC","","N");
		this.getCount("-C(=S)- [aromatic attach]","C(=S)","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-CHO [aromatic attach]","C(=O)([H])","Aromatic","Any","C","SdsCH","","Any");

		this.getCount("-COOH [aliphatic attach]","C(=O)O([H])","Any","Any","C","Any","","Any");
		this.getCount("HC(=O)O-","[H]C(=O)O","Any","Any","C","Any","","Any");
		this.getCount("-C(=O)O- [nitrogen attach]","C(=O)O","Any","Any","C","Any","","N");
		this.getCount("-C(=O)O- [cyclic]","C(=O)O","Any","InRing","C","Any","","Any");
		this.getCount("-C(=O)O- [olefinic attach]","C(=O)O","Olefinic","Any","C","Any","","Any");
		this.getCount("-C(=O)O- [aliphatic attach]","C(=O)O","Any","Any","C","Any","","Any");
		this.getCount("-C(=O)SH [nitrogen attach]","C(=O)S[H]","Any","Any","C","Any","","N");
		this.getCount("-C(=O)SH [aliphatic attach]","C(=O)S[H]","Any","Any","C","Any","","Any");
		this.getCount("-C(=O)S- [nitrogen attach]","C(=O)S","Any","Any","C","Any","","N");
		this.getCount("-C(=O)S- [aliphatic attach]","C(=O)S","Any","Any","C","Any","","Any");
		this.getCount("-C(=S)SH [nitrogen attach]","C(=S)S[H]","Any","Any","C","Any","","N");
		this.getCount("-C(=S)SH [aliphatic attach]","C(=S)S[H]","Any","Any","C","Any","","Any");
		this.getCount("-C(=S)OH [aliphatic attach]","C(=S)O[H]","Any","Any","C","Any","","Any");
		this.getCount("-C(=S)O- [nitrogen attach]","C(=S)O","Any","Any","C","Any","","N");
		this.getCount("-C(=S)O- [aliphatic attach]","C(=S)O","Any","Any","C","Any","","Any");
		this.getCount("-C(=S)S- [nitrogen attach]","C(=S)S","Any","Any","C","Any","","N");
		this.getCount("-C(=S)S- [aliphatic attach]","C(=S)S","Any","Any","C","Any","","Any");
		this.getCount("-C(=S)- [2 nitrogen attach]","C(=S)","Any","Any","C","Any","","2N");
		this.getCount("-C(=S)- [nitrogen, aliphatic attach]","C(=S)","Any","Any","C","Any","","N");
		this.getCount("-C(=S)- [aliphatic attach]","C(=S)","Any","Any","C","Any","","Any");

		this.getCount("-CF3 [aromatic attach]","C(F)(F)F","Aromatic","Any","C","SssssC","","Any");
		this.getCount("-CCl3 [aromatic attach]","C(Cl)(Cl)Cl","Aromatic","Any","C","SssssC","","Any");
		this.getCount("-C#N [aromatic attach]","C#N","Aromatic","Any","C","StsC","","Any");

		this.getCount("-CH=NOH","C([H])=NO([H])","Any","Any","C","Any","","Any");
		this.getCount("-CH=NO","C([H])=NO","Any","Any","C","Any","","Any");

		this.getCount(">C=NOH","C=NO([H])","Any","Any","C","Any","","Any");
		this.getCount(">C=NO","C=NO","Any","Any","C","Any","","Any");

		this.getCount("-C([H])=N[H] [Nitrogen attach]","C([H])=N[H]","Any","Any","C","SdsCH","","N");
		this.getCount("-C([H])=N [Nitrogen attach]","C([H])=N","Any","Any","C","SdsCH","","N");
		this.getCount(">C=N[H] [2 Nitrogen attach]","C=N[H]","Any","Any","C","SdssC","","2N");
		this.getCount(">C=N[H] [Nitrogen attach]","C=N[H]","Any","Any","C","SdssC","","N");
		this.getCount(">C(=N) [2 Nitrogen attach]","C=N","Any","Any","C","SdssC","","2N");
		this.getCount(">C(=N) [Nitrogen attach]","C=N","Any","Any","C","SdssC","","N");

		this.getCount("-CH=NH","C([H])=N([H])","Any","Any","C","SdsCH","","Any");
		this.getCount(">C=NH","C=N([H])","Any","Any","C","SdssC","","Any");
		this.getCount("CH2=N","C([H])([H])=N","Any","Any","C","SdCH2","","Any");
		this.getCount("-CH=N","C([H])=N","Any","Any","C","SdsCH","","Any");
		this.getCount(">C=N","C=N","Any","Any","C","SdssC","","Any");

		this.getCount("-NHN=O","N([H])N=O","Any","Any","N","SssNH","","Any");
		this.getCount(">NN=O","NN=O","Any","Any","N","SsssN","","Any");
		this.getCount("-N(=O)","N(=O)","Any","Any","N","SdsN","","Any");
		this.getCount("-N=NH","N=N([H])","Any","Any","N","SdsN","","Any");
		this.getCount("-N=N-","N=N","Any","Any","N","SdsN","","Any");

		this.getCount("-NH2 [nitrogen attach]","N([H])([H])","Any","Any","N","SsNH2","","N");
		this.getCount("-NH- [nitrogen attach]","N([H])","Any","Any","N","SssNH","","XN");
		this.getCount("-N< [nitrogen attach]","N","Any","Any","N","SsssN","","XN");

		this.getCount("-CHO [nitrogen attach]","C(=O)([H])","Any","Any","C","SdsCH","","N");
		this.getCount("-CHO [aliphatic attach]","C(=O)([H])","Any","Any","C","Any","","Any");

		this.getCount("-C(=O)- [halogen attach]","C(=O)","Any","Any","C","SdssC","","F");
		this.getCount("-C(=O)- [halogen attach]","C(=O)","Any","Any","C","SdssC","","Cl");
		this.getCount("-C(=O)- [halogen attach]","C(=O)","Any","Any","C","SdssC","","Br");
		this.getCount("-C(=O)- [halogen attach]","C(=O)","Any","Any","C","SdssC","","I");
		this.getCount("-C(=O)- [2 aromatic attach]","C(=O)","2Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=O)- [nitrogen, aromatic attach]","C(=O)","Aromatic","Any","C","SdssC","","N");
		this.getCount("-C(=O)- [aromatic attach]","C(=O)","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C(=O)- [2 nitrogen attach]","C(=O)","Any","Any","C","SdssC","","2N");
		this.getCount("-C(=O)- [nitrogen, aliphatic attach]","C(=O)","Any","Any","C","SdssC","","N");
		this.getCount("-C(=O)- [phosphorus attach]","C(=O)","Any","Any","C","SdssC","","P");
		this.getCount("-C(=O)- [olefinic attach]","C(=O)","Olefinic","Any","C","SdssC","","2C");
		this.getCount("C=O(ketone, aliphatic attach)","C(=O)","Any","Any","C","SdssC","","2C");
		this.getCount("C=O(non-ketone, aliphatic attach)","C(=O)","Any","Any","C","SdssC","","Any");

		this.getCount("-CCl3 [aliphatic attach]","C(Cl)(Cl)(Cl)","Any","Any","C","SssssC","","Any");
		this.getCount("-CF3 [aliphatic attach]","C(F)(F)F","Any","Any","C","SssssC","","Any");

		this.getCount("-C#N [aliphatic nitrogen attach]","C#N","Any","Any","C","StsC","","N");
		this.getCount("-C#N [aliphatic sulfur attach]","C#N","Any","Any","C","StsC","","S");
		this.getCount("-C#N [olefinic attach]","C#N","Olefinic","Any","C","StsC","","Any");
		this.getCount("-C#N [aliphatic attach]","C#N","Any","Any","C","StsC","","Any");

		this.getCount("P=NH","P=N[H]","Any","Any","P","Any","","Any");
		this.getCount("P=N","P=N","Any","Any","P","Any","","Any");
		this.getCount("P=S","P=S","Any","Any","P","Any","","Any");
		this.getCount("P(=O)","P(=O)","Any","Any","P","Any","","Any");

		this.getCount("As(=O)","[As](=O)","Any","Any","As","Any","","Any");

		this.getCount("-SH [aliphatic attach]","S([H])","Any","Any","S","SsSH","","Any");
		this.getCount("-S- [three membered ring]","C1CS1","Any","Any","S","SssS","C","Any");
		this.getCount("-S- [sulfur attach]","S","Any","Any","S","SssS","","XS");
		this.getCount("-S- [arsenic attach]","S","Any","Any","S","SssS","","As");
		this.getCount("-S- [phosphorus attach]","S","Any","Any","S","SssS","","XP");
		this.getCount("-S- [nitrogen attach]","S","Any","Any","S","SssS","","XN");
		this.getCount("-S- [2 aromatic attach]","S","2Aromatic","Any","S","SssS","","Any");
		this.getCount("-S- [aromatic attach]","S","Aromatic","Any","S","SssS","","Any");

		this.getCount("-OH [sulfur attach]","O[H]","Any","Any","O","SsOH","","S");
		this.getCount("-OH [oxygen attach]","O[H]","Any","Any","O","SsOH","","O");
		this.getCount("-OH [arsenic attach]","O[H]","Any","Any","O","SsOH","","As");
		this.getCount("-OH [phosphorus attach]","O[H]","Any","Any","O","SsOH","","P");
		this.getCount("-OH [nitrogen attach]","O[H]","Any","Any","O","SsOH","","N");
		this.getCount("-OH [aliphatic attach]","O([H])","Any","Any","O","SsOH","","Any");
		this.getCount("-O- (epoxide)","C1CO1","Any","Any","O","SssO","C","Any");
		this.getCount("-O- [2 aromatic attach]","O","2Aromatic","Any","O","SssO","","Any");
		this.getCount("-O- [phosphorus, aromatic attach]","O","Aromatic","Any","O","SssO","","P");
		this.getCount("-O- [2 phosphorus attach]","O","Any","Any","O","SssO","","2P");
		this.getCount("-O- [phosphorus attach]","O","Any","Any","O","SssO","","P");
		this.getCount("-O- [arsenic attach]","O","Any","Any","O","SssO","","XAs");
		this.getCount("-O- [nitrogen attach]","O","Any","Any","O","SssO","","XN");
		this.getCount("-O- [sulfur attach]","O","Any","Any","O","SssO","","XS");
		this.getCount("-O- [oxygen attach]","O","Any","Any","O","SssO","","XO");
		this.getCount("-O- [aromatic attach]","O","Aromatic","Any","O","SssO","","Any");

		this.getCount("-CH= [aromatic attach]","C([H])","Aromatic","Any","C","SdsCH","","Any");
		this.getCount(">C= [aromatic attach]","C","Aromatic","Any","C","SdssC","","Any");
		this.getCount("-C# [aromatic attach]","C","Aromatic","Any","C","StsC","","Any");
		this.getCount("-CH3 [aromatic attach]","C([H])([H])([H])","Aromatic","Any","C","SsCH3","","Any");
		this.getCount("-CH2- [aromatic attach]","C([H])([H])","Aromatic","Any","C","SssCH2","","Any");
		this.getCount("-CH< [aromatic attach]","C([H])","Aromatic","Any","C","SsssCH","","Any");
		this.getCount(">C< [aromatic attach]","C","Aromatic","Any","C","SssssC","","Any");

		this.getCount("-CH3 [aliphatic attach]","C([H])([H])([H])","Any","Any","C","SsCH3","","Any");
		this.getCount("-CH2- [aliphatic attach]","C([H])[H]","Any","Any","C","SssCH2","","Any");
		this.getCount("-CH< [aliphatic attach]","C([H])","Any","Any","C","SsssCH","","Any");
		this.getCount(">C< [aliphatic attach]","C","Any","Any","C","SssssC","","Any");
		this.getCount("=C=","C","Any","Any","C","SddC","","Any");
		this.getCount("=CH2 [aliphatic attach]","C([H])([H])","Any","Any","C","SdCH2","","Any");
		this.getCount("=CH [aliphatic attach]","C([H])","Any","Any","C","SdsCH","","Any");
		this.getCount("=C [aliphatic attach]","C","Any","Any","C","SdssC","","Any");
		this.getCount("#CH [aliphatic attach]","C([H])","Any","Any","C","StCH","","Any");
		this.getCount("#C [aliphatic attach]","C","Any","Any","C","StsC","","Any");

		this.getCount("-NH2 [aromatic attach]","N([H])([H])","Aromatic","Any","N","SsNH2","","Any");
		this.getCount("-NH- [aromatic attach]","N([H])","Aromatic","Any","N","SssNH","","Any");
		this.getCount("-N< [aromatic attach]","N","Aromatic","Any","N","SsssN","","Any");
		this.getCount("-NH2 [aliphatic attach]","N([H])([H])","Any","Any","N","SsNH2","","Any");
		this.getCount("-NH- [aliphatic attach]","N([H])","Any","Any","N","SssNH","","Any");
		this.getCount("-N< [aliphatic attach]","N","Any","Any","N","SsssN","","Any");

		this.getCount("Fused aromatic carbon","C","3Aromatic","Any","C","SaaaC","","Any");
		this.getCount("Fused aromatic nitrogen","N","3Aromatic","Any","N","SaaaN","","Any");
		this.getCount("ACH","C[H]","Aromatic","Any","C","SaaCH","","Any");
		this.getCount("ANH","N([H])","Aromatic","Any","N","SaaNH","","Any");
		this.getCount("AC","C","Aromatic","Any","C","SsaaC","","Any");
		this.getCount("AN","N","Aromatic","Any","N","SaaN","","Any");
		this.getCount("AN","N","Aromatic","Any","N","SsaaN","","Any");
		this.getCount("AN+","N","Aromatic","Any","N","SsaaNp","","Any");
		this.getCount("AN+","N","Aromatic","Any","N","SaaaNp","","Any");//TODO: wont work unless SaaaNp Estate frag is added!
		this.getCount("AO","O","Aromatic","Any","O","SaaO","","Any");
		this.getCount("AS","S","Aromatic","Any","S","SaaS","","Any");

		this.getCount("=O [other]","O","Any","Any","O","SdO","","Any");
		this.getCount("=S [other]","S","Any","Any","O","SdS","","Any");
		this.getCount("-[O-]","[O-]","Any","Any","O","SsOm","","Any");
		this.getCount("-O- [aliphatic attach]","O","Any","Any","O","SssO","","Any");
		this.getCount("-S- [aliphatic attach]","S","Any","Any","S","SssS","","Any");

		this.getCount("-F [aromatic attach]","F","Aromatic","Any","F","SsF","","Any");
		this.getCount("-Br [aromatic attach]","Br","Aromatic","Any","Br","SsBr","","Any");
		this.getCount("-Cl [aromatic attach]","Cl","Aromatic","Any","Cl","SsCl","","Any");
		this.getCount("-I [aromatic attach]","I","Aromatic","Any","I","SsI","","Any");

		this.getCount("-Br [olefinic attach]","Br","Olefinic","Any","Br","SsBr","","Any");
		this.getCount("-Cl [olefinic attach]","Cl","Olefinic","Any","Cl","SsCl","","Any");
		this.getCount("-F [olefinic attach]","F","Olefinic","Any","F","SsF","","Any");
		this.getCount("-I [olefinic attach]","I","Olefinic","Any","I","SsI","","Any");

		this.getCount("Halogen [Nitrogen attach]","Br","Any","Any","Br","SsBr","","N");
		this.getCount("Halogen [Nitrogen attach]","Cl","Any","Any","Cl","SsCl","","N");
		this.getCount("Halogen [Nitrogen attach]","F","Any","Any","F","SsF","","N");
		this.getCount("Halogen [Nitrogen attach]","I","Any","Any","I","SsI","","N");

		this.getCount("-Br [aliphatic attach]","Br","Any","Any","Br","SsBr","","Any");
		this.getCount("-Cl [aliphatic attach]","Cl","Any","Any","Cl","SsCl","","Any");
		this.getCount("-F [aliphatic attach]","F","Any","Any","F","SsF","","Any");
		this.getCount("-I [aliphatic attach]","I","Any","Any","I","SsI","","Any");

		this.getCount("H [phosphorus attach]","[H]","Any","Any","H","Any","","P");
		this.getCount("H [nitrogen attach]","[H]","Any","Any","H","Any","","N");
		this.getCount("H [silicon attach]","[H]","Any","Any","H","Any","","Si");
		this.getCount("H [carbon attach]","[H]","Any","Any","H","Any","","C");
		this.getCount("H [other]","[H]","Any","Any","H","Any","","Any");

		this.getCount("-Hg","[Hg]","Any","Any","Hg","SsHg","","Any");
		this.getCount("-Hg-","[Hg]","Any","Any","Hg","SssHg","","Any");
		this.getCount("-Hg+","[Hg]","Any","Any","Hg","SsHgp","","Any");

		this.getCount("Hg","[Hg]","Any","Any","Hg","Any","","Any");
		this.getCount("Pb","[Pb]","Any","Any","Pb","Any","","Any");
		this.getCount("B","[B]","Any","Any","B","Any","","Any");

		this.getCount("P [+3 valence, all single bonds]","P","Any","Any","P","SsssP","","Any");
		this.getCount("P [+3 valence, all single bonds]","P","Any","Any","P","SssPH","","Any");
		this.getCount("P [+3 valence, all single bonds]","P","Any","Any","P","SsPH2","","Any");
		this.getCount("P [+5 valence, one double bond]","P","Any","Any","P","SdsssP","","Any");
		this.getCount("P [+5 valence, two double bonds]","P","Any","Any","P","SddsP","","Any");
		this.getCount("P [+5 valence, all single bonds]","P","Any","Any","P","SsssssP","","Any");
		this.getCount("P","P","Any","Any","P","Any","","Any");

		this.getCount("As [+3 valence, one double bond]","[As]","Any","Any","As","SdsAs","","Any");
		this.getCount("As [+3 valence, all single bonds]","[As]","Any","Any","As","SsAsH2","","Any");
		this.getCount("As [+3 valence, all single bonds]","[As]","Any","Any","As","SssAsH","","Any");
		this.getCount("As [+3 valence, all single bonds]","[As]","Any","Any","As","SsssAs","","Any");
		this.getCount("As [+5 valence, one double bond]","[As]","Any","Any","As","SdsssAs","","Any");
		this.getCount("As [+5 valence, two double bonds]","[As]","Any","Any","As","SddsAs","","Any");
		this.getCount("As [+5 valence, all single bonds]","[As]","Any","Any","As","SsssssAs","","Any");
		this.getCount("As","[As]","Any","Any","As","Any","","Any");//TODO change name of frag to "[As]" to avoid name conflict with WHIM "As" descriptor
		
//		if (!CalculateReactiveFlags) return;
//
//		//find flag variables for reactive substructures:
//		this.getCount2("C(=O)O alpha to primary or secondary halogen", "COC(=O)C(Br)([H])", "Any","Any", "C", "SdssC", "", "Any");
//		this.getCount2("C(=O)O alpha to primary or secondary halogen", "COC(=O)C(Cl)([H])", "Any","Any", "C", "SdssC", "", "Any");		
//		this.getCount2("C(=O)O alpha to primary or secondary halogen", "COC(=O)C(F)([H])", "Any","Any", "C", "SdssC", "", "Any");
//		this.getCount2("C(=O)O alpha to primary or secondary halogen", "COC(=O)C(I)([H])", "Any","Any", "C", "SdssC", "", "Any");
//		
//		this.getCount2("C(=O) alpha to primary or secondary halogen", "C(=O)(C)C(Br)([H])", "Any","Any", "C", "SdssC", "", "Any");
//		this.getCount2("C(=O) alpha to primary or secondary halogen", "C(=O)(C)C(Cl)([H])", "Any","Any", "C", "SdssC", "", "Any");		
//		this.getCount2("C(=O) alpha to primary or secondary halogen", "C(=O)(C)C(F)([H])", "Any","Any", "C", "SdssC", "", "Any");
//		this.getCount2("C(=O) alpha to primary or secondary halogen", "C(=O)(C)C(I)([H])", "Any","Any", "C", "SdssC", "", "Any");
//
//		//need at least 1 H attached to carbon that halogen is attached to for primary or secondary
//		this.getCount2("C#N alpha to primary or secondary halogen", "N#CC(Br)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C#N alpha to primary or secondary halogen", "N#CC(Cl)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C#N alpha to primary or secondary halogen", "N#CC(I)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C#N alpha to primary or secondary halogen", "N#CC(F)([H])", "Any","Any", "C", "Any", "", "Any");
//
//		this.getCount2("C#C alpha to primary  or secondary OH", "C#CC(O)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C=C alpha to primary  or secondary OH", "C=CC(O)([H])", "Any","Any", "C", "Any", "", "Any");
//		
//		
//		this.getCount2("C#C alpha to primary or secondary halogen", "C#CC(F)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C#C alpha to primary or secondary halogen", "C#CC(Cl)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C#C alpha to primary or secondary halogen", "C#CC(Br)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C#C alpha to primary or secondary halogen", "C#CC(I)([H])", "Any","Any", "C", "Any", "", "Any");
//
//		this.getCount2("C=C alpha to primary or secondary halogen", "C=CC(F)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C=C alpha to primary or secondary halogen", "C=CC(Cl)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C=C alpha to primary or secondary halogen", "C=CC(Br)([H])", "Any","Any", "C", "Any", "", "Any");
//		this.getCount2("C=C alpha to primary or secondary halogen", "C=CC(I)([H])", "Any","Any", "C", "Any", "", "Any");
		
				
//		this.GetNitroSubstitutionCount(1, 2, -1);
//		this.GetNitroSubstitutionCount(1, 3, -1);
//		this.GetNitroSubstitutionCount(1, 4, -1);
//		
//		this.GetNitroSubstitutionCount(1, 2, 3);
//		this.GetNitroSubstitutionCount(1, 2, 4);
//		this.GetNitroSubstitutionCount(1, 3, 5);
				
//		this.GetCOOHxCH2InteractionCount();
//		this.GetCHOxCH2InteractionCount();
//		this.GetDiesterxCH2InteractionCount();
		
		
}
	
	private void GetCOOHxCH2InteractionCount() {
		String fragdesc1="-COOH [aliphatic attach]";
		String fragdesc2="-CH2- [aliphatic attach]";
					
		double val1=(Double)this.htFragVal.get(fragdesc1);
		double val2=(Double)this.htFragVal.get(fragdesc2);
		
		String fragdesc="COOHxCH2";
		
		double val=val1*val2;
		
		if (val>0) System.out.println(dd.ID+"\t"+fragdesc);
		
		this.htFragVal.put(fragdesc,new Double (val));			
		
	}

	private void GetCHOxCH2InteractionCount() {
		String fragdesc1="-CHO [aliphatic attach]";
		String fragdesc2="-CH2- [aliphatic attach]";
					
		double val1=(Double)this.htFragVal.get(fragdesc1);
		double val2=(Double)this.htFragVal.get(fragdesc2);
		
		String fragdesc="CHOxCH2";
		
		double val=val1*val2;
		
		if (val>0) System.out.println(dd.ID+"\t"+fragdesc);

		this.htFragVal.put(fragdesc,new Double (val));			
		
	}

	private void GetDiesterxCH2InteractionCount() {
		String fragdesc1a="-C(=O)O- [aliphatic attach]";
		String fragdesc1b="-C(=O)O- [aromatic attach]";
		String fragdesc2="-CH2- [aliphatic attach]";
					
		double val1a=(Double)this.htFragVal.get(fragdesc1a);
		double val1b=(Double)this.htFragVal.get(fragdesc1b);
		
		double val2=(Double)this.htFragVal.get(fragdesc2);
		
		String fragdesc="COOxCH2 (diesters)";
		
		double val1=val1a+val1b;
		
		if (val1<2) val1=0;
		
		double val=val1*val2;
		
		if (val>0) System.out.println(dd.ID+"\t"+fragdesc);
			
		this.htFragVal.put(fragdesc,new Double (val));	
		
		
	}

	
	private void FindACAC(IAtomContainer m,IRingSet rs) {
		
		double count=0;
		
		for (int i=0;i<=m.getAtomCount()-1;i++) {
			
			if (Assigned[i]) continue;
			if (!Fragment[i].equals("SsaaC")) continue;
			
			List ca=m.getConnectedAtomsList(m.getAtom(i));
			
			for (int j=0;j<ca.size();j++) {
				
				IAtom caj=(IAtom)ca.get(j);
				
				if (Assigned[m.getAtomNumber(caj)]) {
					continue;
				}
				
				if (!Fragment[m.getAtomNumber(caj)].equals("SsaaC")) continue;
				
				if (!EStateFragmentDescriptor.InSameAromaticRing(m,m.getAtom(i),caj,rs)) {
					Assigned[i]=true;
					AssignedFragment[i] = "ACAC";
					
					Assigned[m.getAtomNumber(caj)]=true;
					AssignedFragment[m.getAtomNumber(caj)] = "ACAC";
					
					count++;
				}
			}

		}
		
		String fragdesc="ACAC";
		this.htFragVal.put(fragdesc, new Double(count));
		dd.MW_Frag+=12.01115*2*count;

		
	}
	
	private void FindANANSameRing(IAtomContainer m,IAtomContainer m2,IRingSet rs) {
		
		double countAN=0;
		double countANH=0;
		
		iloop:
		for (int i=0;i<=m.getAtomCount()-1;i++) {	
			if (Assigned[i]) continue;
			if (!Fragment[i].equals("SsaaN") && !Fragment[i].equals("SaaN") && !Fragment[i].equals("SaaNH")) continue;
			
			List ca=m.getConnectedAtomsList(m.getAtom(i));
			
			jloop:
			for (int j=0;j<ca.size();j++) {
				IAtom caj=(IAtom)ca.get(j);								
				if (!Fragment[m.getAtomNumber(caj)].equals("SsaaN") && !Fragment[m.getAtomNumber(caj)].equals("SaaN") && !Fragment[m.getAtomNumber(caj)].equals("SaaNH")) continue;
				
				if (EStateFragmentDescriptor.InSameAromaticRing(m,m.getAtom(i),caj,rs)) {					
					Assigned[i]=true;
					
					if (Fragment[i].equals("SaaNH")) {
						countANH++;
						AssignedFragment[i] = "ANH [attached to AN in same ring]";
					} else if (Fragment[i].equals("SaaN") || Fragment[i].equals("SsaaN")) {
						countAN++;			
						AssignedFragment[i] = "AN [attached to AN in same ring]";
					}					
					break;
				}
			}
			
						
			if (Fragment[i].equals("SaaNH") && Assigned[i]) {
				// assign hydrogen:
				List ca2=m2.getConnectedAtomsList(m2.getAtom(i));
				
				for (int j=0;j<ca2.size();j++) {					
					IAtom ca2j=(IAtom)ca2.get(j);
					if (ca2j.getSymbol().equals("H")) {
						Assigned[m2.getAtomNumber(ca2j)]=true;						
						break;
					}					
				}				
			}
		}
		
		String fragdesc="AN [attached to AN in same ring]";
		this.htFragVal.put(fragdesc, new Double(countAN));
		dd.MW_Frag+=14.00674*countAN;

		fragdesc="ANH [attached to AN in same ring]";
		this.htFragVal.put(fragdesc, new Double(countANH));
		dd.MW_Frag+=(14.00674+1.00797)*countANH;

		
	}
	
	
	
	private void getCount(String fragdesc, String smiles, String AttachType, String RingType,
			String AttachAtomSymbol, String EStateAttachAtom, 
			String DontAssignSymbol,
			String AttachedToSymbol) {
		
		String checkdesc = "-Cl [aromatic attach]";
		
		IAtomContainer molecule;
		IRingSet ringset;

		molecule=m2;
		ringset=rs2;

		double MWfrag=0;
		
		
		if (!(htFragVal.get(fragdesc) instanceof Double)) {
			//if not assigned a value yet, put in 0:
			this.htFragVal.put(fragdesc, new Double(0));
		}
		
		

		boolean AllAssigned = true;

		for (int i = 0; i <= molecule.getAtomCount() - 1; i++) {

			//if (!molecule.getAtom(i).getSymbol().equals("H")) {
				if (!Assigned[i]) {
					AllAssigned = false;
					break;
				}
			//}
		}

		if (AllAssigned)
			return;


		
		int count = 0;

		try {
			double time3 = System.currentTimeMillis() / 1000.0;
			
			
			SmilesParser   sp  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			AtomContainer fragment = (AtomContainer)sp.parseSmiles(smiles);
			
			double time4 = System.currentTimeMillis() / 1000.0;
			
			
//			if (fragdesc.equals("-NO2 [aromatic attach]")) {
//				System.out.println("here1");
//				
//				for (int i=0;i<fragment.getAtomCount();i++) {
//					IAtom atomi=fragment.getAtom(i);
//					System.out.println(atomi.getSymbol()+"\t"+atomi.getFormalCharge());
//				}
//				System.out.println("");
//				
//				for (int i=0;i<molecule.getAtomCount();i++) {
//					IAtom atomi=molecule.getAtom(i);
//					System.out.println(atomi.getSymbol()+"\t"+atomi.getFormalCharge());
//				}				
//			}


			//			System.out.println((time4 - time3) + " secs");

			for (int i = 0; i < fragment.getAtomCount(); i++) {
				IAtom a = fragment.getAtom(i);
				isf.configure(a);
				
				//System.out.println(a.getSymbol()+"\t"+a.getHydrogenCount());
			}

			

			
			//			System.out.println("");

			List u = null;

			
			try {
				if (fragment.getAtomCount() > 2) {
					List ll = uit.getSubgraphAtomsMaps(molecule, fragment);
					u = MoleculeFragmenter4.getUniqueSubgraphs(ll); // Remove
																	// duplicates
																	// from the
																	// list

				} else if (fragment.getAtomCount() == 2) {
					u = this.GetTwoAtomFragmentLists(molecule, fragment);
				} else { // one atom
					u = this.GetSingleAtomFragmentLists(molecule, fragment);
				}

			} catch (Exception e) {
				return;
			}
			
			//			System.out.println(smiles+"\t"+u.size());

			if (fragdesc.equals(checkdesc)) {
				if (Debug) System.out.println("Count=" + u.size());
			}


			
			if (u.size() == 0)
				return;
			
			MWfrag=0;
			
			iloop: 
			for (int i = 0; i < u.size(); i++) {
				List mlist = (List) u.get(i);
				
				if (i==0) {
					for (int j = 0; j < mlist.size(); j++) { //Remove graphs that contain used atoms
						int atomnumber = (Integer) mlist.get(j);
						if (atomnumber>-1) {
							String s=molecule.getAtom(atomnumber).getSymbol();
							if (!s.equals(DontAssignSymbol)) {
								MWfrag+=ap.GetMass(s);
							}
							
//							MWfrag+=m.getAtom(atomnumber).getExactMass();
						}
						
					}
					
//					System.out.println(fragdesc+"\t"+MWfrag);	
				}
				

				if (mlist.size() == 0)
					continue;

				for (int j = 0; j < mlist.size(); j++) { //Remove graphs that contain used atoms

					int atomnumber = (Integer) mlist.get(j);

//					if (molecule.getAtom(atomnumber).getSymbol().equals("H"))
//						continue;

					if ((Assigned[atomnumber] == true)) {
						continue iloop;
					}

				}//for j

				//					if (smiles.equals("S(=O)(=O)N")) System.out.println("bob0");

				int AttachmentAtomNumber=this.GetAttachmentAtomNumber(AttachAtomSymbol,EStateAttachAtom,molecule,mlist);
				
						
				if (AttachmentAtomNumber==-1 && !EStateAttachAtom.equals("Any")) { // added !EStateAttachAtom.equals("Any") on 3/18/10
					if (Debug) {
						if (fragdesc.equals(checkdesc)) {
							System.out.println("failed estate check");
						}
					}
					continue;
				}

//				if (smiles.equals("S=C=S")) {
//					System.out.println(u.size());
//				}


				String at = this.GetAttachmentType(AttachmentAtomNumber, molecule,mlist);
				String InARing = this.InARing(AttachmentAtomNumber, molecule,
						ringset, mlist);


				boolean HaveAttached = false;

				if (!AttachedToSymbol.equals("Any")) {
					HaveAttached = this
							.IsAttachedToElement(AttachmentAtomNumber,molecule,
									mlist, AttachedToSymbol);
				}

				
				
				if (!AttachedToSymbol.equals("Any") && !HaveAttached) {
					if (fragdesc.equals(checkdesc))
						if (Debug) System.out.println("failed attached to symbol");
					continue;
				}

				if (!RingType.equals("Any") && !RingType.equals(InARing)) {
					if (fragdesc.equals(checkdesc))
						if (Debug) System.out.println("failed ring type");
					continue;
				}

				
				// note Aromatic can be Aromatic,TwoAromatic, ThreeAromatic, AromaticOlefinic, etc.
				
				
				if (!AttachType.equals("Any") && !at.equals(AttachType)) {
					
					if (AttachType.equals("Aromatic") && at.indexOf("Aromatic")>-1) {
						// do nothing
					} else {
						if (fragdesc.equals(checkdesc)) {
							if (Debug) System.out.println("failed attach type2");
							if (Debug) System.out.println(at + "\t" + AttachType);
						}

						continue;
						
					}
					
				}

								
				for (int j = 0; j < mlist.size(); j++) { //Remove graphs that contain used atoms            			
					int atomnumber = (Integer) mlist.get(j);
//					if (molecule.getAtom(atomnumber).getSymbol().equals("H"))
//						continue;
//					if (molecule.getAtom(atomnumber).getFlag(
//							CDKConstants.ISAROMATIC)) {
//						continue iloop;
//					}
					if (molecule.getAtom(atomnumber).getSymbol().equals(
							DontAssignSymbol))
						continue;
				}//for j

				count++;

				// assign atoms        			
				for (int j = 0; j < mlist.size(); j++) { //Remove graphs that contain used atoms            			
					int atomnumber = (Integer) mlist.get(j);
//					if (molecule.getAtom(atomnumber).getSymbol().equals("H"))
//						continue;
					if (molecule.getAtom(atomnumber).getSymbol().equals(
							DontAssignSymbol))
						continue;
					Assigned[atomnumber] = true;
					AssignedFragment[atomnumber] = fragdesc;
				}//for j

			}// for i
			
			if (count==0) return;
			
//			System.out.println("Count="+count);
			
			if (htFragVal.get(fragdesc) instanceof Double) {
			
				dd.MW_Frag+=count*MWfrag;

				double OldCount = (Double)(htFragVal.get(fragdesc));
				count += OldCount;
				this.htFragVal.put(fragdesc, new Double(count));
			} else {
				dd.MW_Frag+=count*MWfrag;
				this.htFragVal.put(fragdesc, new Double(count));
			}
			
			
//			if (Debug) System.out.println(fragdesc+"\t"+count+"\t"+MWfrag);
//			System.out.println(fragdesc+"\t"+count+"\t"+MWfrag);

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	private void GetNitroSubstitutionCount(int pos1, int pos2,int pos3) {
		String fragdesc="";
		
		int NumNitroGroups=-1;

		if (pos3>-1) {
			fragdesc=pos1+","+pos2+","+pos3+"-trinitro";
			NumNitroGroups=3;
		} else {
			fragdesc=pos1+","+pos2+"-dinitro";
			NumNitroGroups=2;
		}

		if (!(htFragVal.get(fragdesc) instanceof Double)) {
			//if not assigned a value yet, put in 0:
			this.htFragVal.put(fragdesc, new Double(0));
		}

		
		IAtomContainer molecule;
		IRingSet ringset;

		molecule=m2;
		ringset=rs2;
		
		try {

			int nitrocount=0;
			
			if (htFragVal.get("-NO2 [aromatic attach]") instanceof Double) {
				double dnc=(Double)htFragVal.get("-NO2 [aromatic attach]");
				nitrocount=(int)dnc;
			}
			
					
			if (nitrocount<2) return;
			
			if (pos3>-1) {
				if (nitrocount<3) return;
			}
			

			for (int i=0;i<rs2.getAtomContainerCount();i++) {
				IAtomContainer ring=rs2.getAtomContainer(i);
				
				if (ring.getAtomCount()!=6) continue;
				
				java.util.List AtomsWithNitro=this.GetRingAtomsWithNitroGroups(ring,m2);
				
				Collections.sort(AtomsWithNitro);
				
				int count=AtomsWithNitro.size();
				
				if (count!=NumNitroGroups) continue;
				

				if (NumNitroGroups==2) {
					int val1=(Integer)AtomsWithNitro.get(0);
					int val2=(Integer)AtomsWithNitro.get(1);
					
					if (val1==pos1 && val2==pos2) {
						this.htFragVal.put(fragdesc, new Double(1.0));
						System.out.println(dd.ID+"\t"+fragdesc);
					}
				}

				if (NumNitroGroups==3) {
					int val1=(Integer)AtomsWithNitro.get(0);
					int val2=(Integer)AtomsWithNitro.get(1);
					int val3=(Integer)AtomsWithNitro.get(2);
					
					if (val1==pos1 && val2==pos2  && val3==pos3) {
						this.htFragVal.put(fragdesc, new Double(1.0));
						System.out.println(dd.ID+"\t"+fragdesc);
					}
				}
				
//				System.out.println(count);
				
				
				
			}
			
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}
	
	private List GetRingAtomsWithNitroGroups(IAtomContainer ring,IAtomContainer molecule) {
		
		java.util.LinkedList<Integer> list=new java.util.LinkedList(); 
		
		iloop:
		for (int i=0;i<ring.getAtomCount();i++) {
			IAtom ai=ring.getAtom(i);
			
			java.util.List ca=molecule.getConnectedAtomsList(ai);
			
			for (int j=0;j<ca.size();j++) {
				IAtom caj=(IAtom)ca.get(j);
				
				if (!ring.contains(caj)) {
					if (caj.getSymbol().equals("N")) {
						boolean IsNitro=IsNitro(caj,molecule);
						if (IsNitro) {
							list.add(new Integer(molecule.getAtomNumber(ai)));
							continue iloop;
						}
					}
				}
			}
			
			
		}
		
		Collections.sort(list);
		
		int Min=list.get(0);
		
		for (int k=0;k<list.size();k++) {
			int val=list.get(k)-Min+1;
			list.set(k, new Integer(val));
		}
		
		return list;
	}
	
	boolean IsNitro(IAtom atom,IAtomContainer molecule) {
		
		if (atom.getFormalCharge()!=1) return false;
		
		boolean Have1=false;
		boolean Have2=false;
		
		java.util.List ca=molecule.getConnectedAtomsList(atom);
		
		for (int j=0;j<ca.size();j++) {
			IAtom caj=(IAtom)ca.get(j);
			if (!caj.getSymbol().equals("O")) continue;
			
			if (caj.getFormalCharge()==-1) {
				if (molecule.getBond(atom, caj).getOrder().equals(IBond.Order.SINGLE)) {
					Have1=true;
				}
			} else if (caj.getFormalCharge()==0) {
				if (molecule.getBond(atom, caj).getOrder().equals(IBond.Order.DOUBLE)) {
					Have2=true;
				} 
			}
		}

		return Have1 && Have2;
		
	}
	
	//used for finding things like propargylic alcohols (flag variables not fragments)
	
	private void getCount2(String fragdesc, String smiles, String AttachType, String RingType,
			String AttachAtomSymbol, String EStateAttachAtom, 
			String DontAssignSymbol,
			String AttachedToSymbol) {
		
		String checkdesc = "-Cl [aromatic attach]";
		
		IAtomContainer molecule;
		IRingSet ringset;

		molecule=m2;
		ringset=rs2;

		double MWfrag=0;
		
		
		if (!(htFragVal.get(fragdesc) instanceof Double)) {
			//if not assigned a value yet, put in 0:
			this.htFragVal.put(fragdesc, new Double(0));
		}
		

//		boolean AllAssigned = true;
//
//		for (int i = 0; i <= molecule.getAtomCount() - 1; i++) {
//
//			//if (!molecule.getAtom(i).getSymbol().equals("H")) {
//				if (!Assigned[i]) {
//					AllAssigned = false;
//					break;
//				}
//			//}
//		}
//
//		if (AllAssigned)
//			return;

		int count = 0;

		try {
			double time3 = System.currentTimeMillis() / 1000.0;
			
			SmilesParser   sp  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			IAtomContainer fragment = (IAtomContainer)sp.parseSmiles(smiles);
			double time4 = System.currentTimeMillis() / 1000.0;


			for (int i = 0; i < fragment.getAtomCount(); i++) {
				IAtom a = fragment.getAtom(i);
				isf.configure(a);
			}

			List u = null;

			try {
				if (fragment.getAtomCount() > 2) {
					List ll = uit.getSubgraphAtomsMaps(	molecule, fragment);
					u = MoleculeFragmenter4.getUniqueSubgraphs(ll); // Remove
																	// duplicates
																	// from the
																	// list

				} else if (fragment.getAtomCount() == 2) {
					u = this.GetTwoAtomFragmentLists(molecule, fragment);
				} else { // one atom
					u = this.GetSingleAtomFragmentLists(molecule, fragment);
				}

			} catch (Exception e) {
				return;
			}
			
			//			System.out.println(smiles+"\t"+u.size());

			if (fragdesc.equals(checkdesc)) {
				if (Debug) System.out.println("Count=" + u.size());
			}

			if (u.size() == 0)
				return;

			MWfrag=0;
			
			iloop: 
			for (int i = 0; i < u.size(); i++) {
				List mlist = (List) u.get(i);
				
//				if (i==0) {
//					for (int j = 0; j < mlist.size(); j++) { //Remove graphs that contain used atoms
//						int atomnumber = (Integer) mlist.get(j);
//						if (atomnumber>-1) {
//							String s=molecule.getAtom(atomnumber).getSymbol();
//							if (!s.equals(DontAssignSymbol)) {
//								MWfrag+=ap.GetMass(s);
//							}
//							
////							MWfrag+=m.getAtom(atomnumber).getExactMass();
//						}
//						
//					}
//					
////					System.out.println(fragdesc+"\t"+MWfrag);	
//				}
				

				if (mlist.size() == 0)
					continue;

//				for (int j = 0; j < mlist.size(); j++) { //Remove graphs that contain used atoms
//
//					int atomnumber = (Integer) mlist.get(j);
//
////					if (molecule.getAtom(atomnumber).getSymbol().equals("H"))
////						continue;
//
////					if ((Assigned[atomnumber] == true)) {
////						continue iloop;
////					}
//
//				}//for j

				//					if (smiles.equals("S(=O)(=O)N")) System.out.println("bob0");

				int AttachmentAtomNumber=this.GetAttachmentAtomNumber(AttachAtomSymbol,EStateAttachAtom,molecule,mlist);
				
				if (AttachmentAtomNumber==-1) {
					if (Debug) {
						if (fragdesc.equals(checkdesc)) {
							System.out.println("failed estate check");
						}
					}
					continue;
				}
				
				
				String at = this.GetAttachmentType(AttachmentAtomNumber, molecule,mlist);
				String InARing = this.InARing(AttachmentAtomNumber, molecule,
						ringset, mlist);

				
				boolean HaveAttached = false;

				if (!AttachedToSymbol.equals("Any")) {
					HaveAttached = this
							.IsAttachedToElement(AttachmentAtomNumber,molecule,
									mlist, AttachedToSymbol);
				}

				
				
				if (!AttachedToSymbol.equals("Any") && !HaveAttached) {
					if (fragdesc.equals(checkdesc))
						if (Debug) System.out.println("failed attached to symbol");
					continue;
				}

				if (!RingType.equals("Any") && !RingType.equals(InARing)) {
					if (fragdesc.equals(checkdesc))
						if (Debug) System.out.println("failed ring type");
					continue;
				}

				
				// note Aromatic can be Aromatic,TwoAromatic, ThreeAromatic, AromaticOlefinic, etc.
				
				
				if (!AttachType.equals("Any") && !at.equals(AttachType)) {
					
					if (AttachType.equals("Aromatic") && at.indexOf("Aromatic")>-1) {
						// do nothing
					} else {
						if (fragdesc.equals(checkdesc)) {
							if (Debug) System.out.println("failed attach type2");
							if (Debug) System.out.println(at + "\t" + AttachType);
						}

						continue;
						
					}
					
				}

				count++;


			}// for i
			
			if (count==0) return;
			else {
				this.htFragVal.put(fragdesc, new Double(1));	
			}
			
			
//			if (Debug) System.out.println(fragdesc+"\t"+1+"\t"+MWfrag);
			System.out.println(dd.ID+"\t"+fragdesc);
			

		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	int GetAttachmentAtomNumber(String symbol,String EState, IAtomContainer molecule, List mlist) {
		
		
		for(int j=0;j<mlist.size();j++){ 
			int tempatomnumber=(Integer)mlist.get(j);
			
			if (molecule.getAtom(tempatomnumber).getSymbol().equals(symbol)) {
				if (EState.equals("Any")) {
					return tempatomnumber;
				} else if (Fragment[tempatomnumber].equals(EState)) {					
					return tempatomnumber;
					
				}
			}			
    	}//for j
		return -1;
	}
	
	
	
	
	 String GetAttachmentAtomEStateFragment(String symbol,IAtomContainer molecule, List mlist) {
			for(int j=0;j<mlist.size();j++){ //Remove graphs that contain used atoms
				int tempatomnumber=(Integer)mlist.get(j);
				if (molecule.getAtom(tempatomnumber).getSymbol().equals(symbol)) {
					return Fragment[tempatomnumber];				
				}			
	    	}//for j
			return "";
		}
	
	boolean IsAttachedToElement(int AttachmentAtomNumber,IAtomContainer molecule, List mlist,String attachsymbol) {
			
		IAtom ai=molecule.getAtom(AttachmentAtomNumber);
		List ca=molecule.getConnectedAtomsList(ai);
		
		int attachcount=1;
		
		int count=0;
		
		
		if (attachsymbol.indexOf("2") >-1) {
			attachcount=2;
			attachsymbol=attachsymbol.substring(1,attachsymbol.length());
		} else if (attachsymbol.indexOf("3") >-1) {
			attachcount=3;
			attachsymbol=attachsymbol.substring(1,attachsymbol.length());
		} else if (attachsymbol.indexOf("X") >-1) { // any number of this element attached
			attachcount=-1;
			attachsymbol=attachsymbol.substring(1,attachsymbol.length());
		}
		
		
		//TODO: if attached atom is aromatic should we not return true?
		
			for (int j = 0; j < ca.size(); j++) {
				IAtom caj=(IAtom)ca.get(j);
				
				if (!InFragment(mlist, molecule.getAtomNumber(caj))) {

					if (caj.getSymbol().toUpperCase().equals(attachsymbol.toUpperCase())) {// have
						if (!caj.getFlag(CDKConstants.ISAROMATIC)) {
							count++;
						}
						
						// lower
//						if (ca[j].getFlag(CDKConstants.ISAROMATIC)) {
//							if (ca[j].getSymbol().equals(
//									attachsymbol.toUpperCase()))
//								count++;
//						}
//					} else {
//						if (ca[j].getSymbol().equals(attachsymbol))
//							count++;
//					}
						
					}
						
				} // end if (!fragment.contains(ca[j])) {
			} // and ca for loop
		
			
//		if (attachsymbol.equals("Cl")) System.out.println(count+"\t"+attachcount);	
		
		if (count>0 && attachcount==-1) return true;
			
		if (count==attachcount) return true;		
		else return false;
	}

	boolean InFragment(List atomlist,int AtomNumber) {
		for (int ii=0;ii<atomlist.size();ii++) {
			int num=(Integer)atomlist.get(ii);
			if (num==AtomNumber) {
				return true;				
			}
		}
		
		return false;
		
	}
	private String InARing(int AttachmentAtomNumber,IAtomContainer molecule,IRingSet ringset,List mlist) {		
		
		if (AttachmentAtomNumber==-1) return "NotInRing";
//		System.out.println("enter in a ring");
		
		IAtom ai=molecule.getAtom(AttachmentAtomNumber);
		
		//System.out.println("atomnumber = "+atomnumber);
		
		
		for (int i=0;i<=ringset.getAtomContainerCount()-1;i++) {		
			IRing r=(IRing)ringset.getAtomContainer(i);
							
			if (r.contains(ai)) return "InRing";
		}
		
		
		return "NotInRing";
		
	}
	
	private List GetTwoAtomFragmentLists(IAtomContainer molecule, IAtomContainer fragment) {

		ArrayList u = new ArrayList();

		IAtom a0 = fragment.getAtom(0);
		IAtom a1 = fragment.getAtom(1);
		IBond.Order BondOrder = fragment.getBond(a0, a1).getOrder();

		//System.out.println(a0.getSymbol()+"\t"+a1.getSymbol()+"\t"+BondOrder);

		iloop: for (int i = 0; i < molecule.getAtomCount(); i++) {
			IAtom ai = molecule.getAtom(i);

			int AtomNum1 = molecule.getAtomNumber(ai);

			if (ai.getSymbol().equals(a0.getSymbol())) {
				List ca = molecule.getConnectedAtomsList(ai);

				for (int j = 0; j < ca.size(); j++) {
					IAtom caj=(IAtom)ca.get(j);
					
					int AtomNum2 = molecule.getAtomNumber(caj);
					if (caj.getSymbol().equals(a1.getSymbol())) {
						if (molecule.getBond(ai, caj).getOrder() == BondOrder) {
							if (AtomNum1 < AtomNum2) {
								ArrayList al = new ArrayList();
								al.add(new Integer(AtomNum1));
								al.add(new Integer(AtomNum2));
								u.add(al);
							} else {
								ArrayList al = new ArrayList();
								al.add(new Integer(AtomNum2));
								al.add(new Integer(AtomNum1));
								u.add(al);
							}

							continue iloop;

						}
					}
				}

			}
		}

		//TODO- add code for removing possible duplicates

		//		for (int i=0;i<=u.size()-1;i++) {
		//			ArrayList al=(ArrayList)u.get(i);
		//			
		//			for (int j=0;j<=al.size()-1;j++) {
		//				System.out.println(al.get(j)+"\t");
		//			}
		//			System.out.println("");
		//		}

		return u;
	}

	private List GetSingleAtomFragmentLists(IAtomContainer molecule, IAtomContainer fragment) {

		ArrayList u = new ArrayList();

		IAtom a0 = fragment.getAtom(0);

		iloop: for (int i = 0; i < molecule.getAtomCount(); i++) {
			IAtom ai = molecule.getAtom(i);

			int AtomNum1 = molecule.getAtomNumber(ai);

			if (ai.getSymbol().equals(a0.getSymbol())) {
				ArrayList al = new ArrayList();
				al.add(new Integer(AtomNum1));
				u.add(al);
			}
		}

		return u;
	}

	public static List getUniqueSubgraphs(List subgraphs) {
        ArrayList unique = new ArrayList();
//        System.out.println("getting unique subgraphs");
        for (int i = 0; i < subgraphs.size(); i++) {
            List current = (List) subgraphs.get(i);
            ArrayList ids = new ArrayList();
//            System.out.print(i+": ");
            for (Iterator iter = current.iterator(); iter.hasNext();) {
                RMap rmap = (RMap) iter.next();
                ids.add( new Integer(rmap.getId1()));
//                System.out.print((rmap.getId1()+1)+" ");
            }
//            System.out.println();
            Collections.sort(ids);
            unique.add(ids);
        }
        HashSet hs = new HashSet(unique);
        unique = new ArrayList();
        for (Iterator iter = hs.iterator(); iter.hasNext();) {
            unique.add(iter.next());
        }
        return unique;

    }
	
	
	
	public Hashtable Calculate(IAtomContainer m, IAtomContainer m2,IRingSet rs,IRingSet rs2,DescriptorData dd, String[] Fragment) {

//		System.out.println(Fragment.length);
		
		//this.m=molecule;
		dd.MW_Frag=0;
		
		this.m2 = m2;
		this.dd = dd;
		this.Fragment = Fragment;
		this.rs2=rs2;
		

		Assigned = new boolean[m2.getAtomCount()];
		AssignedFragment = new String[m2.getAtomCount()];

		for (int i = 0; i <= m2.getAtomCount() - 1; i++) {
			Assigned[i] = false;
			AssignedFragment[i] = "";
		}

		htFragVal = new Hashtable();

		this.FindACAC(m,rs);
		this.FindANANSameRing(m,m2,rs);
		this.FindFrags();
		
		double MW=this.Calculate_mw(m);
		double errMW=Math.abs(MW-dd.MW_Frag);
		java.text.DecimalFormat df=new java.text.DecimalFormat("0.00");
		
		if (errMW>0.1) {
			System.out.println(dd.ID+"\t"+df.format(MW)+"\t"+df.format(dd.MW_Frag)+"\t"+df.format(errMW));
		}

		
		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			if (Debug)  {
				System.out.println((i + 1) + "\t" + Fragment[i] + "\t"
					+ Assigned[i] + "\t" + AssignedFragment[i]);
			}
		}
		
		if (Debug) System.out.println("MW calc\t"+dd.MW_Frag);
		
		
		
		if (Debug) {
			for (java.util.Enumeration e = htFragVal.keys(); e.hasMoreElements();) {
				String strVar = (String) e.nextElement();
				double val = (Double) htFragVal.get(strVar);
				System.out.println(strVar+"\t"+val);
			}
		}
		
		return htFragVal;

	}
	
	private double Calculate_mw(IAtomContainer m) {
		// tried to use CDK built in methods but they suck
		// alternative method would be to use m2 which includes the hydrogens
		
		try {
			ap=AtomicProperties.getInstance();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		double MW=0;
				
		for (int i=0;i<=m.getAtomCount()-1;i++) {			
			IAtom a=m.getAtom(i);
			MW+=ap.GetMass(a.getSymbol());
			MW+=a.getImplicitHydrogenCount()*ap.GetMass("H");
			
		}
		return MW;			
	}

	private String GetAttachmentType(int AttachmentAtomNumber, IAtomContainer molecule, List mlist) {

		String at = "";

		if (AttachmentAtomNumber==-1) return "";// 3/18/10
		
		IAtom ai = molecule.getAtom(AttachmentAtomNumber);
		List ca = molecule.getConnectedAtomsList(ai);

		int AromaticCount = 0;
		int OlefinicCount = 0;

		for (int j = 0; j < ca.size(); j++) {

			//int cajatomnumber = molecule.getAtomNumber(ca[j]);
			IAtom caj=(IAtom)ca.get(j);
			
			if (!InFragment(mlist, molecule.getAtomNumber(caj))) {
				if (caj.getFlag(CDKConstants.ISAROMATIC)) {
					AromaticCount++;
				} else {
					if (caj.getSymbol().equals("C")) {
						List ca2 = molecule.getConnectedAtomsList(caj);

						for (int k = 0; k < ca2.size(); k++) {
							IAtom ca2k=(IAtom)ca2.get(k);
							if (!ca2k.equals(ai)) {
								if (ca2k.getSymbol().equals("C")
										&& molecule.getBond(caj, ca2k)
												.getOrder().equals(IBond.Order.DOUBLE)) {
									OlefinicCount++;
								}
							}
						}

					}

				}// end else

			} // end if (!fragment.contains(ca[j])) {
		} // and ca for loop
		
//		System.out.println(AttachmentAtomNumber+"\t"+AromaticCount);
		
		

		if (AromaticCount >=2) 
			return AromaticCount+"Aromatic";
		else if (AromaticCount == 1) {
			if (OlefinicCount == 0)
				return "Aromatic";
			else
				return "AromaticOlefinic";
		} else {
			if (OlefinicCount == 0)
				return "Aliphatic";
			else if (OlefinicCount >= 1)
				return "Olefinic";
			
		}

		return "";

	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
