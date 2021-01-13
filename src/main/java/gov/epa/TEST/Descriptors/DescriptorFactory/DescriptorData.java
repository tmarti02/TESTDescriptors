package gov.epa.TEST.Descriptors.DescriptorFactory;

import java.lang.reflect.Field;
import java.util.*;
import gov.epa.TEST.Descriptors.DescriptorUtilities.Utilities;
import java.io.*;
import java.text.*;


public class DescriptorData {

	//TODO- store as JsonObject directly in the future? each descriptor class could be a separate object

	// acceptable elements:
	// C, H, O, N, F, Cl, Br, I, S, P, Si, As, Hg, Sn

	public boolean LoadRaghu3dDescriptions=false;
	public boolean LoadMy3dDescriptions=false;
	public boolean Load2dDescriptions=true;
	
	public String ID=""; // chemical abstract number
	public String Error;
	public String InChiKey;
	public String InChi;
	public String InChi_Warning;

	public boolean ThreeD=true;
	public String SoftwareVersion;

	// ***************************************************************************
	// Constitutional descriptors
	public double MW;  // molecular weight
	public double AMW; // average molecular weight
	public double Sv;  // sum of atomic van der Waals volumes (scaled on Carbon atom)
	public double Se;  // sum of atomic Sanderson electronegativities (scaled on Carbon atom)
	public double Sp;  // sum of atomic polarizabilities (scaled on Carbon atom)
	public double Ss;  // sum of Kier-Hall electrotopological states
	public double Mv;  // mean atomic van der Waals volume (scaled on Carbon atom)
	public double Me;  // mean atomic Sanderson electronegativity (scaled on Carbon atom)
	public double Mp;  // mean atomic polarizability (scaled on Carbon atom)
	public double Ms;  // mean electrotopological state
	public int nAT;    // number of atoms
	public int nSK;    // number of non-H atoms
	public int nBT;    // number of bonds
	public int nBO;    // number of non-H bonds
	public int nBM;    // number of multiple bonds
	public double SCBO;// sum of conventional bond orders (H-depleted)
	public double ARR; // aromatic ratio
	public int nCIC;   // number of rings
	public int nCIR;   // number of circuits
	public int nDB;    // number of double bonds
	public int nTB;    // number of triple bonds
	public int nAB;    // number of aromatic bonds
	public int nH;     // number of Hydrogen atoms
	public int nC;     // number of Carbon atoms
	public int nN;     // number of Nitrogen atoms
	public int nO;     // number of Oxygen atoms
	public int nP;     // number of Phosphorous atoms
	public int nS;     // number of Sulfur atoms
	public int nF;     // number of Fluorine atoms
	public int nCL;    // number of Chlorine atoms
	public int nBR;    // number of Bromine atoms
	public int nI;     // number of Iodine atoms
	public int nB;     // number of Boron atoms
	public int nX;     // number of halogen atoms
	public int nR03;   // number of 3-membered rings
	public int nR04;   // number of 4-membered rings
	public int nR05;   // number of 5-membered rings
	public int nR06;   // number of 6-membered rings
	public int nR07;   // number of 7-membered rings
	public int nR08;   // number of 8-membered rings
	public int nR09;   // number of 9-membered rings
	public int nR10;   // number of 10-membered rings
	public int nR11;   // number of 11-membered rings
	public int nR12;   // number of 12-membered rings
	public int nBnz;   // number of benzene-like rings


	// *************************************************************************
	//	Chi Connectivity Indices (46)

	 //A. Simple
	public double x0;     // Simple zero order chi index
	public double x1;     // Simple 1st order chi index
	public double x2;     // Simple 2nd order chi index
	public double xp3;    // Simple 3rd order path chi index
	public double xp4;    // Simple 4th order path chi index
	public double xp5;    // Simple 5th order path chi index
	public double xp6;    // Simple 6th order path chi index
	public double xp7;    // Simple 7th order path chi index
	public double xp8;    // Simple 8th order path chi index
	public double xp9;    // Simple 9th order path chi index
	public double xp10;   // Simple 10th order path chi index
	public double xc3;    // Simple 3rd order cluster chi index
	public double xc4;    // Simple 4th order cluster chi index
	public double xpc4;   // Simple 4th order path/cluster chi index
	public double xch3;   // Simple 3rd order chain chi index
	public double xch4;   // Simple 4th order chain chi index
	public double xch5;   // Simple 5th order chain chi index
	public double xch6;   // Simple 6th order chain chi index
	public double xch7;   // Simple 7th order chain chi index
	public double xch8;   // Simple 8th order chain chi index
	public double xch9;   // Simple 9th order chain chi index
	public double xch10;  // Simple 10th order chain chi index
	public double knotp;  // Difference between chi cluster-3 and chi path/cluster-4

//	B. Valence
	public double xv0;    // Valence zero order chi index
	public double xv1;    // Valence 1st order chi index
	public double xv2;    // Valence 2nd order chi index
	public double xvp3;   // Valence 3rd order path chi index
	public double xvp4;   // Valence 4th order path chi index
	public double xvp5;   // Valence 5th order path chi index
	public double xvp6;   // Valence 6th order path chi index
	public double xvp7;   // Valence 7th order path chi index
	public double xvp8;   // Valence 8th order path chi index
	public double xvp9;   // Valence 9th order path chi index
	public double xvp10;  // Valence 10th order path chi index
	public double xvc3;   // Valence 3rd order cluster chi index
	public double xvc4;   // Valence 4th order cluster chi index
	public double xvpc4;  // Valence 4th order path/cluster chi index
	public double xvch3;  // Valence 3rd order chain chi index
	public double xvch4;  // Valence 4th order chain chi index
	public double xvch5;  // Valence 5th order chain chi index
	public double xvch6;  // Valence 6th order chain chi index
	public double xvch7;  // Valence 7th order chain chi index
	public double xvch8;  // Valence 8th order chain chi index
	public double xvch9;  // Valence 9th order chain chi index
	public double xvch10; // Valence 10th order chain chi index
	public double knotpv; // Difference between chi valence cluster-3 and chi valence path/cluster-4.



//	 *****************************************************************************
	//Kappa Shape Indices (7)

	public double k1;    // First Order Kappa Shape Index
	public double k2;    // Second Order Kappa Shape Index
	public double k3;    // Third Order Kappa Shape Index
	public double ka1;   // First Order Kappa  Alpha Shape Index
	public double ka2;   // Second Order Kappa Alpha Shape Index
	public double ka3;   // Third Order Kappa Alpha Shape Index
	public double phia;  // Kappa flexibility index

// *****************************************************************************
//	Electrotopological State Indices (62)

	//TODO: add SaaaNp

	//A.	Non-Hydrogen Indices (E-state sum)
	public double SsCH3;   // Sum of all ( – CH3 ) E-state values in the molecule
	public double SdCH2;   // Sum of all ( = CH2 ) E-state values in the molecule
	public double SssCH2;  // Sum of all ( – CH2 – ) E-state values in the molecule
	public double StCH;    // Sum of all ( tCH ) E-state values in the molecule
	public double SdsCH;   // Sum of all ( = CH – ) E-state values in the molecule
	public double SaaCH;   // Sum of all ( aaCH  ) E-state values in the molecule
	public double SsssCH;  // Sum of all ( > CH – ) E-state values in the molecule
	public double SddC;    // Sum of all ( = C = ) E-state values in the molecule
	public double StsC;    // Sum of all ( tC – ) E-state values in the molecule
	public double SdssC;   // Sum of all ( = C < ) E-state values in the molecule
	public double SsaaC;   // Sum of all ( saaC  ) E-state values in the molecule (SaasC in MDL QSAR)
	public double SaaaC;   // Sum of all ( aaaC ) E-state values in the molecule
	public double SssssC;  // Sum of all ( > C < ) E-state values in the molecule
	public double SsNH2;   // Sum of all ( – NH2 ) E-state values in the molecule
	public double SdNH;    // Sum of all ( = NH ) E-state values in the molecule
	public double SssNH;   // Sum of all ( – NH – ) E-state values in the molecule
	public double SaaNH;   // Sum of all ( aaNH ) E-state values in the molecule
	public double StN;     // Sum of all ( tN ) E-state values in the molecule
	public double SdsN;    // Sum of all ( = N – ) E-state values in the molecule
	public double SaaN;    // Sum of all ( aaN ) E-state values in the molecule
	public double SsssN;   // Sum of all ( > N – ) E-state values in the molecule
	public double SaaaN;   // Sum of all ( aaaN ) E-state values in the molecule (not in MDL QSAR)
	public double SsaaN;   // Sum of all ( saaN ) E-state values in the molecule (not in MDL QSAR)
//	public double StdN;    // Sum of all ( # N = ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SssaaN;  // Sum of all ( ssaaN ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SdsssN;  // Sum of all ( dsssN ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SdaaN;   // Sum of all ( daaN ) E-state values in the molecule
//	public double SddsN;   // Sum of all ( ddsN ) E-state values in the molecule
//	public double SsssssN; // Sum of all ( sssssN ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SssssNH; // Sum of all ( ssssNH ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SdssNH;  // Sum of all ( dssNH ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SssaaNH; // Sum of all ( ssaaNH ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SsaaNH;  // Sum of all ( saaNH ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SdsNH2;  // Sum of all ( = NH2 – ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SdNH3;   // Sum of all ( = NH3 ) E-state values in the molecule (not in MDL QSAR, omit?)
//	public double SddNH;   // Sum of all ( = NH = ) E-state values in the molecule (not in MDL QSAR, omit?)
	public double SsNH3p;  // Sum of all ( – NH3+ ) E-state values in the molecule
	public double SssNH2p; // Sum of all ( – NH2+ – ) E-state values in the molecule
	public double SsssNHp; // Sum of all ( > NH+ – ) E-state values in the molecule
	public double SddNp;   // Sum of all ( = N+ = ) E-state values in the molecule (not in MDL QSAR)
	public double StsNp;   // Sum of all ( # N+ – ) E-state values in the molecule (not in MDL QSAR)
	public double SdNm;    // Sum of all ( = N - ) E-state values in the molecule (not in MDL QSAR)
	public double SdssNp;  // Sum of all ( = N+ < ) E-state values in the molecule (not in MDL QSAR)
	public double SssssNp; // Sum of all ( > N+ < ) E-state values in the molecule
	public double SsaaNp;  // Sum of all ( saaN+ ) E-state values in the molecule (not in MDL QSAR)
	public double SsOH;    // Sum of all ( – OH ) E-state values in the molecule
	public double SdO;     // Sum of all ( = O ) E-state values in the molecule
	public double SssO;    // Sum of all ( – O – ) E-state values in the molecule
	public double SaaO;    // Sum of all ( aaO ) E-state values in the molecule
	public double SsOm;    // Sum of all ( – O- ) E-state values in the molecule (not in MDL QSAR)
	public double SsF;     // Sum of all ( – F ) E-state values in the molecule
	public double SsCl;    // Sum of all ( – Cl ) E-state values in the molecule
	public double SsBr;    // Sum of all ( – Br ) E-state values in the molecule
	public double SsI;     // Sum of all ( – I ) E-state values in the molecule
	public double SsSiH3;  // Sum of all ( – SiH3 ) E-state values in the molecule
	public double SssSiH2; // Sum of all ( – SiH2 – ) E-state values in the molecule
	public double SsssSiH; // Sum of all ( > SiH – ) E-state values in the molecule
	public double SssssSi; // Sum of all ( > Si < ) E-state values in the molecule
	public double SsssP;   // Sum of all ( > P – ) E-state values in the molecule
	public double SdsssP;  // Sum of all (dsssP ) E-state values in the molecule
	public double SddsP;   // Sum of all (ddsP ) E-state values in the molecule (not in MDL QSAR)
	public double SsssssP; // Sum of all ( sssssP ) E-state values in the molecule
	public double SdssPH;  // Sum of all ( =PH< ) E-state values in the molecule (not in MDL QSAR)
	public double SsPH2;   // Sum of all ( – PH2 ) E-state values in the molecule
	public double SssPH;   // Sum of all ( – PH – ) E-state values in the molecule
	public double SsSH;    // Sum of all ( – SH ) E-state values in the molecule
	public double SdS;     // Sum of all ( = S ) E-state values in the molecule
	public double SssS;    // Sum of all ( – S – ) E-state values in the molecule
	public double SaaS;    // Sum of all ( aaS ) E-state values in the molecule
	public double SdssS;   // Sum of all ( = S < ) E-state values in the molecule
	public double SddssS;  // Sum of all ( ddssS ) E-state values in the molecule
	public double SssssssS;// Sum of all ( ssssssS ) E-state values in the molecule
//	public double SdaaS;   // Sum of all ( daaS ) E-state values in the molecule (not in MDL QSAR, omit???)
	public double SsSnH3;  // Sum of all ( – SnH3 ) E-state values in the molecule
	public double SssSnH2; // Sum of all ( – SnH2 – ) E-state values in the molecule
	public double SsssSnH;  // Sum of all ( > SnH – ) E-state values in the molecule
	public double SssssSn; // Sum of all ( > Sn < ) E-state values in the molecule
	public double SsPbH3;  // Sum of all ( – PbH3 ) E-state values in the molecule
	public double SssPbH2; // Sum of all ( – PbH2 – ) E-state values in the molecule
	public double SsssPbH; // Sum of all ( > PbH – ) E-state values in the molecule
	public double SssssPb; // Sum of all ( > Pb < ) E-state values in the molecule
	public double SsssAs;  // Sum of all ( – As < ) E-state values in the molecule
	public double SdsssAs; // Sum of all ( dsssAs ) E-state values in the molecule
	public double SddsAs;  // Sum of all ( ddsAs ) E-state values in the molecule (not in MDL QSAR, omit???)
	public double SsssssAs;// Sum of all ( sssssAs ) E-state values in the molecule
	public double SsAsH2;  // Sum of all (  – AsH2 ) E-state values in the molecule
	public double SssAsH;  // Sum of all ( – AsH – ) E-state values in the molecule
	public double SsHgp;   // Sum of all ( – Hg+) E-state values in the molecule
	public double SssHg;   // Sum of all ( – Hg –) E-state values in the molecule
	public double SsHg;    // Sum of all ( – Hg) E-state values in the molecule
	public double StCm; // sum of all ->S<- E-state values

	//B.	Non-Hydrogen Indices (fragment count)
	public int SsCH3_acnt;   // Count of all ( – CH3 ) groups in molecule
	public int SdCH2_acnt;   // Count of all ( = CH2 ) groups in molecule
	public int SssCH2_acnt;  // Count of all ( – CH2 – ) groups in molecule
	public int StCH_acnt;    // Count of all ( tCH ) groups in molecule
	public int SdsCH_acnt;   // Count of all ( = CH – ) groups in molecule
	public int SaaCH_acnt;   // Count of all ( aaCH  ) groups in molecule
	public int SsssCH_acnt;  // Count of all ( > CH – ) groups in molecule
	public int SddC_acnt;    // Count of all ( = C = ) groups in molecule
	public int StsC_acnt;    // Count of all ( tC – ) groups in molecule
	public int SdssC_acnt;   // Count of all ( = C < ) groups in molecule
	public int SsaaC_acnt;   // Count of all ( saaC  ) groups in molecule (SaasC in MDL QSAR)
	public int SaaaC_acnt;   // Count of all ( aaaC ) groups in molecule
	public int SssssC_acnt;  // Count of all ( > C < ) groups in molecule
	public int SsNH2_acnt;   // Count of all ( – NH2 ) groups in molecule
	public int SdNH_acnt;    // Count of all ( = NH ) groups in molecule
	public int SssNH_acnt;   // Count of all ( – NH – ) groups in molecule
	public int SaaNH_acnt;   // Count of all ( aaNH ) groups in molecule
	public int StN_acnt;     // Count of all ( tN ) groups in molecule
	public int SdsN_acnt;    // Count of all ( = N – ) groups in molecule
	public int SaaN_acnt;    // Count of all ( aaN )  groups in molecule
	public int SsssN_acnt;   // Count of all ( > N – ) groups in molecule
	public int SaaaN_acnt;   // Count of all ( aaaN )  groups in molecule (not in MDL QSAR)
	public int SsaaN_acnt;   // Count of all ( saaN )  groups in molecule (not in MDL QSAR)
//	public int StdN_acnt;    // Count of all ( # N = ) groups in molecule (not in MDL QSAR, omit?)
//	public int SssaaN_acnt;  // Count of all ( ssaaN ) groups in molecule (not in MDL QSAR, omit?)
//	public int SdsssN_acnt;  // Count of all ( dsssN ) groups in molecule (not in MDL QSAR, omit?)
//	public int SdaaN_acnt;   // Count of all ( daaN ) groups in molecule
//	public int SddsN_acnt;   // Count of all ( ddsN ) groups in molecule
//	public int SsssssN_acnt; // Count of all ( sssssN ) groups in molecule (not in MDL QSAR, omit?)
//	public int SssssNH_acnt; // Count of all ( ssssNH ) groups in molecule (not in MDL QSAR, omit?)
//	public int SdssNH_acnt;  // Count of all ( dssNH ) groups in molecule (not in MDL QSAR, omit?)
//	public int SssaaNH_acnt; // Count of all ( ssaaNH ) groups in molecule (not in MDL QSAR, omit?)
//	public int SsaaNH_acnt;  // Count of all ( saaNH ) groups in molecule (not in MDL QSAR, omit?)
//	public int SdsNH2_acnt;  // Count of all ( = NH2 – ) groups in molecule (not in MDL QSAR, omit?)
//	public int SdNH3_acnt;   // Count of all ( = NH3 ) groups in molecule (not in MDL QSAR, omit?)
//	public int SddNH_acnt;   // Count of all ( = NH = ) groups in molecule (not in MDL QSAR, omit?)
	public int SsNH3p_acnt;  // Count of all ( – NH3+ ) groups in molecule
	public int SssNH2p_acnt; // Count of all ( – NH2+ – ) groups in molecule
	public int SsssNHp_acnt; // Count of all ( > NH+ – ) groups in molecule
	public int SddNp_acnt;   // Count of all ( = N+ = ) groups in molecule (not in MDL QSAR)
	public int StsNp_acnt;   // Count of all ( # N+ – ) groups in molecule (not in MDL QSAR)
	public int SdNm_acnt;    // Count of all ( = N - ) groups in molecule (not in MDL QSAR)
	public int SdssNp_acnt;  // Count of all ( = N+ < ) groups in molecule (not in MDL QSAR)
	public int SssssNp_acnt; // Count of all ( > N+ < ) groups in molecule
	public int SsaaNp_acnt;  // Count of all ( saaN+ ) groups in molecule (not in MDL QSAR)
	public int SsOH_acnt;    // Count of all ( – OH ) groups in molecule
	public int SdO_acnt;     // Count of all ( = O ) groups in molecule
	public int SssO_acnt;    // Count of all ( – O – ) groups in molecule
	public int SaaO_acnt;    // Count of all ( aaO )   groups in molecule
	public int SsOm_acnt;    // Count of all ( – O- ) groups in molecule (not in MDL QSAR)
	public int SsF_acnt;     // Count of all ( – F ) groups in molecule
	public int SsCl_acnt;    // Count of all ( – Cl ) groups in molecule
	public int SsBr_acnt;    // Count of all ( – Br ) groups in molecule
	public int SsI_acnt;     // Count of all ( – I ) groups in molecule
	public int SsSiH3_acnt;  // Count of all ( – SiH3 ) groups in molecule
	public int SssSiH2_acnt; // Count of all ( – SiH2 – ) groups in molecule
	public int SsssSiH_acnt; // Count of all ( > SiH – ) groups in molecule
	public int SssssSi_acnt; // Count of all ( > Si < ) groups in molecule
	public int SsssP_acnt;   // Count of all ( > P – ) groups in molecule
	public int SdsssP_acnt;  // Count of all (dsssP ) groups in molecule
	public int SddsP_acnt;   // Count of all (ddsP ) groups in molecule (not in MDL QSAR)
	public int SsssssP_acnt; // Count of all ( sssssP ) groups in molecule
	public int SdssPH_acnt;  // Count of all ( =PH< ) groups in molecule (not in MDL QSAR)
	public int SsPH2_acnt;   // Count of all ( – PH2 ) groups in molecule
	public int SssPH_acnt;   // Count of all ( – PH – ) groups in molecule
	public int SsSH_acnt;    // Count of all ( – SH ) groups in molecule
	public int SdS_acnt;     // Count of all ( = S ) groups in molecule
	public int SssS_acnt;    // Count of all ( – S – ) groups in molecule
	public int SaaS_acnt;    // Count of all ( aaS ) groups in molecule
	public int SdssS_acnt;   // Count of all ( = S < ) groups in molecule
	public int SddssS_acnt;  // Count of all ( ddssS ) groups in molecule
	public int SssssssS_acnt;// Count of all ( ssssssS ) groups in molecule
//	public int SdaaS_acnt;   // Count of all ( daaS ) groups in molecule (not in MDL QSAR, omit???)
	public int SsSnH3_acnt;  // Count of all ( – SnH3 ) groups in molecule
	public int SssSnH2_acnt; // Count of all ( – SnH2 – ) groups in molecule
	public int SsssSnH_acnt;  // Count of all ( > SnH – ) groups in molecule
	public int SssssSn_acnt; // Count of all ( > Sn < ) groups in molecule
	public int SsPbH3_acnt;  // Count of all ( – PbH3 ) groups in molecule
	public int SssPbH2_acnt; // Count of all ( – PbH2 – ) groups in molecule
	public int SsssPbH_acnt; // Count of all ( > PbH – ) groups in molecule
	public int SssssPb_acnt; // Count of all ( > Pb < ) groups in molecule
	public int SsssAs_acnt;  // Count of all ( – As < ) groups in molecule
	public int SdsssAs_acnt; // Count of all ( dsssAs ) groups in molecule
	public int SddsAs_acnt;  // Count of all ( ddsAs ) groups in molecule (not in MDL QSAR, omit???)
	public int SsssssAs_acnt;// Count of all ( sssssAs ) groups in molecule
	public int SsAsH2_acnt;  // Count of all (  – AsH2 ) groups in molecule
	public int SssAsH_acnt;  // Count of all ( – AsH – ) groups in molecule
	public int SsHgp_acnt;   // Count of all (-Hg+) groups in a molecule
	public int SssHg_acnt;   // Count of all (-Hg-) groups in a molecule
	public int SsHg_acnt;    // Count of all (-Hg) groups in a molecule
	public int StCm_acnt; // Count of all ->S<- E-state values

	// C. Hydrogen Indices
	public double SHsOH;  // Sum of the Hydrogen E-State values for all [– OH ] groups in the molecule.
	public double SHdNH;  // Sum of  all Hydrogen E-State values for all [= NH ] groups in the molecule.
	public double SHsSH;  // Sum of  all Hydrogen E-State values for all [– SH ] groups in the molecule.
	public double SHsNH2; // Sum of  all Hydrogen E-State values for all [– NH2 ] groups in the molecule.
	public double SHssNH; // Sum of  all Hydrogen E-State values for all [– NH – ] groups in the molecule.
	public double SHtCH;  // Sum of  all Hydrogen E-State values for all [ tCH ] groups in the molecule.
	public double SHother;// Sum of all Hydrogen E-State values for all atom types in the molecule other than specified in the list of H E-state indices.
	public double SHCHnX; // Sum of the Hydrogen E-State values for all CH or CH2 groups with a -F or -Cl also bonded to the carbon.

//	D.	Maximum and Minimum Values
	public double Hmax;   // Maximum hydrogen E-State value in molecule.
	public double Gmax;   // Maximum E-State value in molecule.
	public double Hmin;   // Minimum hydrogen E-State value in molecule.
	public double Gmin;   // Minimum E-State value in molecule.
	public double Hmaxpos;// Maximum positive hydrogen E-State value in molecule.
	public double Hminneg;// Minimum negative hydrogen E-State value in molecule.

	// **************************************************************************
	// Topological Descriptors
	public double ZM1;    // First Zagreb index M1
	public double ZM1V;   // First Zagreb index by valence vertex degrees
	public double ZM2;    // Second Zagreb index M2
	public double ZM2V;   // Second Zagreb index by valence vertex degrees
	public double J;      // Balaban distance connectivity index (pg 21 Todeschini book)
	public double Jt;     // Jt index (pg 21 Todeschini book)
	public double BAC;    // Balaban centric index
	public double Lop;    // Lopping centric index
	public double ICR;    // radial centric information index
	public double TIE;    // E-state topological parameter
	public double MAXDN;  // maximal electrotopological negative variation
	public double MAXDP;  // maximal electrotopological positive variation
	public double DELS;   // molecular electrotopological variation
	public double W;      // Wiener Index
	public double WA;     // Mean Wiener Index

	// **************************************************************************
	// Molecular fragments
	public Hashtable <String,Double>FragmentList;//used to store name and count of all fragments in a molecule

	public static String [] strFragments={//used to output fragment in organized fashion
			//Arsenic:
			"As [+5 valence, one double bond]",
			"As [+5 valence, all single bonds]",
			"As [+5 valence, two double bonds]",
			"As [+3 valence, one double bond]",
			"As [+3 valence, all single bonds]",
			"As",

			//Boron:
			"B",

			//Bromine:
			"-Br [aliphatic attach]",
			"-Br [aromatic attach]",
			"-Br [olefinic attach]",

			//Carbon:
			"-CH3 [aliphatic attach]",
			"-CH3 [aromatic attach]",
			"-CH2- [aliphatic attach]",
			"-CH2- [aromatic attach]",
			"-CH< [aliphatic attach]",
			"-CH< [aromatic attach]",
			">C< [aliphatic attach]",
			">C< [aromatic attach]",
			"=CH2 [aliphatic attach]",
			"-CH= [aromatic attach]",
			"=CH [aliphatic attach]",
			"=C [aliphatic attach]",
			">C= [aromatic attach]",
			"#C [aliphatic attach]",
			"-C# [aromatic attach]",
			"#CH [aliphatic attach]",
			"=C=",
			"ACH",
			"AC",
			"Fused aromatic carbon",
			"ACAC",

			//chlorine:
			"-Cl [aliphatic attach]",
			"-Cl [aromatic attach]",
			"-Cl [olefinic attach]",

			//Fluorine:
			"-F [aliphatic attach]",
			"-F [aromatic attach]",
			"-F [olefinic attach]",

			//Hydrogen:
			"H [carbon attach]",
			"H [phosphorus attach]",
			"H [nitrogen attach]",
			"H [silicon attach]",
			"H [other]",

			//Iodine:
			"-I [aliphatic attach]",
			"-I [aromatic attach]",
			"-I [olefinic attach]",

			//lead:
			"Pb",

			//mercury:
			"-Hg",
			"-Hg-",
			"-Hg+",
			"Hg",

			//Oxygen
			"-O- [nitrogen attach]",
			"-O- [phosphorus attach]",
			"-O- [arsenic attach]",
			"-O- [aromatic attach]",
			"-O- [oxygen attach]",
			"-O- [sulfur attach]",
			"-O- (epoxide)",
			"-O- [2 aromatic attach]",
			"-O- [phosphorus, aromatic attach]",
			"-O- [2 phosphorus attach]",
			"-O- [aliphatic attach]",
			"=O [other]",
			"-OH [phosphorus attach]",
			"-OH [sulfur attach]",
			"-OH [arsenic attach]",
			"-OH [aliphatic attach]",
			"-OH [nitrogen attach]",
			"-OH [oxygen attach]",
			"-OH [aromatic attach]",
			"AO",
			"-[O-]",

			//Nitrogen:
			"-NH2 [aliphatic attach]",
			"-NH2 [aromatic attach]",
			"-NH2 [nitrogen attach]",
			"-NH2 [attached to P]",
			"-NH- [aliphatic attach]",
			"-NH- [aromatic attach]",
			"-NH- [3 membered ring]",
			"-NH- [attached to P]",
			"-NH- [nitrogen attach]",
			"-N< [aliphatic attach]",
			"-N< [aromatic attach]",
			"-N< [nitrogen attach]",
			"-N< [attached to P]",
			"-N< [3 membered ring]",
			"N+ [four single bonds]",
			"AN+",
			"AN",
			"ANH",
			"ANH [attached to AN in same ring]",
			"AN [attached to AN in same ring]",
			"Fused aromatic nitrogen",
			"-[N+]#N",
			"-N=[N+]=[N-]",
			"-N=NH",
			"-N=N-",

			//Phosphorus:
			"P [+5 valence, two double bonds]",
			"P [+5 valence, one double bond]",
			"P [+5 valence, all single bonds]",
			"P [+3 valence, all single bonds]",
			"P",

			//Sulfur
			"-SH [aliphatic attach]",
			"-SH [aromatic attach]",
			"-S- [aromatic attach]",
			"-S- [three membered ring]",
			"-S- [2 aromatic attach]",
			"-S- [nitrogen attach]",
			"-S- [sulfur attach]",
			"-S- [aliphatic attach]",
			"-S- [arsenic attach]",
			"-S- [phosphorus attach]",
			"=S [other]",
			"->S<-",
			"AS",

			//Silicon:
			"Si [aromatic attach]",
			"Si [oxygen, aromatic attach]",
			"Si [oxygen attach]",
			"Si",

			//Tin:
			"Sn [oxygen attach]",
			"Sn [aromatic attach]",
			"Sn",

			// *****************************************
			//Multi atom fragments:
			//carbon + nitrogen
			"-C#N [aliphatic attach]",
			"-C#N [aromatic attach]",
			"-C#N [olefinic attach]",
			"CH2=N",
			"-CH=NH",
			"-CH=N",
			">C=NH",
			">C=N",
			">C(=N) [Nitrogen attach]",
			"[C-]#[N+]",
			"C=[N+]=[N-]",
			">C(=N) [2 Nitrogen attach]",
			"-C([H])=N [Nitrogen attach]",
			"-C([H])=N[H] [Nitrogen attach]",
			">C=N[H] [2 Nitrogen attach]",
			">C=N[H] [Nitrogen attach]",
			"N=C=N",
			"HN=C=N",
			"-C#N [aliphatic sulfur attach]",
			"-C#N [aliphatic nitrogen attach]",

			//Carbon + oxygen:
			"C=O(ketone, aliphatic attach)",
			"-C(=O)- [aromatic attach]",
			"-C(=O)- [phosphorus attach]",
			"-C(=O)- [halogen attach]",
			"-C(=O)- [2 nitrogen attach]",
			"-C(=O)- [nitrogen, aliphatic attach]",
			"-C(=O)- [olefinic attach]",
			"-C(=O)- [nitrogen, aromatic attach]",
			"-C(=O)- [2 aromatic attach]",
			"C=O(non-ketone, aliphatic attach)",
			"-COOH [aliphatic attach]",
			"-COOH [aromatic attach]",
			"-CHO [aliphatic attach]",
			"-CHO [aromatic attach]",
			"-CHO [nitrogen attach]",
			"-C(=O)O- [olefinic attach]",
			"-C(=O)O- [aromatic attach]",
			"-C(=O)O- [aliphatic attach]",
			"-C(=O)O- [nitrogen attach]",
			"-C(=O)O- [cyclic]",
			"HC(=O)O-",
			"-C(=O)OC(=O)-",
			"-C(=O)[O-]",
			"-OC(=O)O-",
			"C=C=O",
			"CH2=C(CH3)C(=O)O-",
			"CH2=CHC(=O)O-",

			//carbon + sulfur:
			"-C(=S)- [aliphatic attach]",
			"-C(=S)- [aromatic attach]",
			"-C(=S)- [nitrogen, aliphatic attach]",
			"-C(=S)- [nitrogen, aromatic attach]",
			"-C(=S)- [2 nitrogen attach]",
			"-C(=S)S- [aliphatic attach]",
			"-C(=S)S- [aromatic attach]",
			"-C(=S)S- [nitrogen attach]",
			"-C(=S)SH [aliphatic attach]",
			"-C(=S)SH [aromatic attach]",
			"-C(=S)SH [nitrogen attach]",
			"-SC(=S)S-",
			"S=C=S",

			//carbon + oxygen + sulfur:
			"-C(=S)OH [aliphatic attach]",
			"-C(=S)OH [aromatic attach]",
			"-C(=S)O- [aliphatic attach]",
			"-C(=S)O- [aromatic attach]",
			"-C(=S)O- [nitrogen attach]",
			"-C(=O)SH [aliphatic attach]",
			"-C(=O)SH [aromatic attach]",
			"-C(=O)SH [nitrogen attach]",
			"-C(=O)S- [aliphatic attach]",
			"-C(=O)S- [aromatic attach]",
			"-C(=O)S- [nitrogen attach]",
			"-OC(=S)O-",
			"-SC(=S)O-",
			"-SC(=O)S-",
			"-SC(=O)O-",

			//carbon + nitrogen + oxygen
			">C=NOH",
			">C=NO",
			"-CH=NOH",
			"-CH=NO",
			"-N=C=O [aliphatic attach]",
			"-N=C=O [aromatic attach]",
			"C=[N+][O-]",

			//nitrogen + oxygen:
			"-ONO2",
			"-N=[N+][O-]",
			"-NO2 [aliphatic attach]",
			"-NO2 [aromatic attach]",
			"-NO2 [nitrogen attach]",
			"-NO2 [olefinic attach]",
			"A[N+][O-]",
			"-NHN=O",
			">NN=O",
			"-N=O [aromatic attach]",
			"-N(=O)",

			//oxygen + sulfur:
			"-S(=O)- [aliphatic attach]",
			"-S(=O)- [aromatic attach]",
			"-S(=O)- [nitrogen, aromatic attach]",
			"-S(=O)- olefinic attach]",
			"-S(=O)- [2 nitrogen attach]",
			"-S(=O)- [nitrogen, aliphatic attach]",
			"-S(=O)(=O)- [nitrogen, aromatic attach]",
			"-S(=O)(=O)- [nitrogen, aliphatic attach]",
			"-S(=O)(=O)- [olefinic attach]",
			"-S(=O)(=O)- [2 nitrogen attach]",
			"-S(=O)(=O)- [aromatic attach]",
			"-S(=O)(=O)- [aliphatic attach]",

			//double bonded phosphorus:
			"P(=O)",
			"P=N",
			"P=NH",
			"P=S",

			//misc halogen
			"-CF3 [aliphatic attach]",
			"-CF3 [aromatic attach]",
			"-CCl3 [aromatic attach]",
			"-CCl3 [aliphatic attach]",
			"Halogen [Nitrogen attach]",

			//misc
			"As(=O)",
			"-N=C=S",
			"Sn=O",
			"-N=S=O"
	};

	public double MW_Frag=0;//Molecular weight from MoleculeFragmenter4

	/*

	*/

	// **************************************************************************
	// 2-D Molecular Properties (9)

	//nvx removed since it is already done in nSK:
//	public double nvx;    // Number of Graph Vertices (number of nonhydrogen atoms)

	public double numHBa; // Number of Hydrogen bond acceptors
	public double numHBd; // Number of Hydrogen bond donors
	public double numwHBd;// Number of weak Hydrogen bond donors (i.e. –CHX, where X = Cl, F)
	public double SHHBa;  // Sum of E-state indices for Hydrogen bond acceptors
	public double SHHBd;  // Sum of E-state indices for Hydrogen bond donors
	public double SHwHBd; // Sum of E-state indices for weak Hydrogen bond donors
	public double Qs;     // Molecular and Group Polarity Index
	public double Qv;     // Molecular and Group Polarity Index
	public double Qsv;    // Average of Qs and Qv


	//Information Indexes (34)
	//A.	Connectivity Based (8)
	public double ic;     // Information Content
	public double si;     // Shannon Information
	public double ib;     // Information Bond Index
	public double I;      // Total Information Content
	public double maxic;  // Maximum Information Content, maxIC
	public double ssi;    // Standardized Shannon Information
	public double R;      // Brillouin Redundancy Index
	public double eim;    // Mean Information Content on the Edge Magnitude

	//B.	Topologically Based (26)
	public double iadje;  // Total Information Content on the Adjacency Equality (iadje)
	public double iadjem; // Mean Information Content on the Adjacency Equality
	public double iadjm;  // Total Information Content on the Adjacency Magnitude
	public double iadjmm; // Mean Information Content on the Adjacency Magnitude
	public double ivdem;  // Mean Information Content of the Vertex Degree Equality
	public double ivdmm;  // Mean Information Content of the Vertex Degree Magnitude
	public double ieadje; // Total Information Content on the Edge Adjacency Equality
	public double ieadjem;// Mean Information Content on the Edge Adjacency Equality
	public double ieadjm; // Total Information Content on the Edge Adjacency Magnitude
	public double ieadjmm;// Mean Information Content on the Edge Adjacency Magnitude
	public double ide;    // Total Information Content on the Distance Equality
	public double idem;   // Mean Information Content on the Distance Equality
	public double idm;    // Total Information Content on the Distance Magnitude
	public double idmm;   // Mean Information Content on the Distance Magnitude
	public double iddem;  // Mean information content on the distance degree equality
	public double iddmm;  // Mean information content on the distance degree magnitude
	public double iede;   // Total Information Content on the Edge Distance Equality
	public double iedem;  // Mean Information Content on the Edge Distance Equality
	public double We;     // Edge Wiener Index
	public double iedm;   // Total Information Content on the Edge Distance Magnitude
	public double iedmm;  // Mean Information Content on the Edge Distance Magnitude
	public double tvc;    // Total Vertex Cyclicity
	public double icyce;  // Total Information on the Vertex Cycle Matrix Equality
	public double icycem; // Mean Information on the Vertex Cycle Matrix Equality
	public double icycm;  // Total Information on the Vertex Cycle Matrix Magnitude
	public double icycmm; // Mean Information on the Vertex Cycle Matrix Magnitude
	public double tti;    // Total Topological Index
	public double ttvi;   // Total Topological Valence Index

	// **************************************************************************
	//Molecular Distance-Edge Vector (19)
	public double MDEC11; //  molecular distance edge between all primary carbons
	public double MDEC12; //  molecular distance edge between all primary and secondary carbons
	public double MDEC13; //  molecular distance edge between all primary and tertiary carbons
	public double MDEC14; //  molecular distance edge between all primary and quaternary carbons
	public double MDEC22; //  molecular distance edge between all secondary carbons
	public double MDEC23; //  molecular distance edge between all secondary and tertiary carbons
	public double MDEC24; //  molecular distance edge between all secondary and quaternary carbons
	public double MDEC33; //  molecular distance edge between all tertiary carbons
	public double MDEC34; //  molecular distance edge between all tertiary and quaternary carbons
	public double MDEC44; //  molecular distance edge between all quaternary carbons
	public double MDEO11; //  molecular distance edge between all primary oxygens
	public double MDEO12; //  molecular distance edge between all primary and secondary oxygens
	public double MDEO22; //  molecular distance edge between all secondary oxygens
	public double MDEN11; //  molecular distance edge between all primary nitrogens
	public double MDEN12; //  molecular distance edge between all primary and secondary nitrogens
	public double MDEN13; //  molecular distance edge between all primary and tertiary nitrogens
	public double MDEN22; //  molecular distance edge between all secondary nitrogens
	public double MDEN23; //  molecular distance edge between all secondary and tertiary nitrogens
	public double MDEN33; //  molecular distance edge between all tertiary nitrogens


	// **************************************************************************
	// Burden eigenvalue descriptors (64)
	public double BEHm1; //  highest eigenvalue n. 1 of Burden matrix / weighted by atomic masses
	public double BEHm2; //  highest eigenvalue n. 2 of Burden matrix / weighted by atomic masses
	public double BEHm3; //  highest eigenvalue n. 3 of Burden matrix / weighted by atomic masses
	public double BEHm4; //  highest eigenvalue n. 4 of Burden matrix / weighted by atomic masses
	public double BEHm5; //  highest eigenvalue n. 5 of Burden matrix / weighted by atomic masses
	public double BEHm6; //  highest eigenvalue n. 6 of Burden matrix / weighted by atomic masses
	public double BEHm7; //  highest eigenvalue n. 7 of Burden matrix / weighted by atomic masses
	public double BEHm8; //  highest eigenvalue n. 8 of Burden matrix / weighted by atomic masses
	public double BELm1; //  lowest eigenvalue n. 1 of Burden matrix / weighted by atomic masses
	public double BELm2; //  lowest eigenvalue n. 2 of Burden matrix / weighted by atomic masses
	public double BELm3; //  lowest eigenvalue n. 3 of Burden matrix / weighted by atomic masses
	public double BELm4; //  lowest eigenvalue n. 4 of Burden matrix / weighted by atomic masses
	public double BELm5; //  lowest eigenvalue n. 5 of Burden matrix / weighted by atomic masses
	public double BELm6; //  lowest eigenvalue n. 6 of Burden matrix / weighted by atomic masses
	public double BELm7; //  lowest eigenvalue n. 7 of Burden matrix / weighted by atomic masses
	public double BELm8; //  lowest eigenvalue n. 8 of Burden matrix / weighted by atomic masses
	public double BEHv1; //  highest eigenvalue n. 1 of Burden matrix / weighted by atomic van der Waals volumes
	public double BEHv2; //  highest eigenvalue n. 2 of Burden matrix / weighted by atomic van der Waals volumes
	public double BEHv3; //  highest eigenvalue n. 3 of Burden matrix / weighted by atomic van der Waals volumes
	public double BEHv4; //  highest eigenvalue n. 4 of Burden matrix / weighted by atomic van der Waals volumes
	public double BEHv5; //  highest eigenvalue n. 5 of Burden matrix / weighted by atomic van der Waals volumes
	public double BEHv6; //  highest eigenvalue n. 6 of Burden matrix / weighted by atomic van der Waals volumes
	public double BEHv7; //  highest eigenvalue n. 7 of Burden matrix / weighted by atomic van der Waals volumes
	public double BEHv8; //  highest eigenvalue n. 8 of Burden matrix / weighted by atomic van der Waals volumes
	public double BELv1; //  lowest eigenvalue n. 1 of Burden matrix / weighted by atomic van der Waals volumes
	public double BELv2; //  lowest eigenvalue n. 2 of Burden matrix / weighted by atomic van der Waals volumes
	public double BELv3; //  lowest eigenvalue n. 3 of Burden matrix / weighted by atomic van der Waals volumes
	public double BELv4; //  lowest eigenvalue n. 4 of Burden matrix / weighted by atomic van der Waals volumes
	public double BELv5; //  lowest eigenvalue n. 5 of Burden matrix / weighted by atomic van der Waals volumes
	public double BELv6; //  lowest eigenvalue n. 6 of Burden matrix / weighted by atomic van der Waals volumes
	public double BELv7; //  lowest eigenvalue n. 7 of Burden matrix / weighted by atomic van der Waals volumes
	public double BELv8; //  lowest eigenvalue n. 8 of Burden matrix / weighted by atomic van der Waals volumes
	public double BEHe1; //  highest eigenvalue n. 1 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BEHe2; //  highest eigenvalue n. 2 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BEHe3; //  highest eigenvalue n. 3 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BEHe4; //  highest eigenvalue n. 4 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BEHe5; //  highest eigenvalue n. 5 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BEHe6; //  highest eigenvalue n. 6 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BEHe7; //  highest eigenvalue n. 7 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BEHe8; //  highest eigenvalue n. 8 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BELe1; //  lowest eigenvalue n. 1 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BELe2; //  lowest eigenvalue n. 2 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BELe3; //  lowest eigenvalue n. 3 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BELe4; //  lowest eigenvalue n. 4 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BELe5; //  lowest eigenvalue n. 5 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BELe6; //  lowest eigenvalue n. 6 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BELe7; //  lowest eigenvalue n. 7 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BELe8; //  lowest eigenvalue n. 8 of Burden matrix / weighted by atomic Sanderson electronegativities
	public double BEHp1; //  highest eigenvalue n. 1 of Burden matrix / weighted by atomic polarizabilities
	public double BEHp2; //  highest eigenvalue n. 2 of Burden matrix / weighted by atomic polarizabilities
	public double BEHp3; //  highest eigenvalue n. 3 of Burden matrix / weighted by atomic polarizabilities
	public double BEHp4; //  highest eigenvalue n. 4 of Burden matrix / weighted by atomic polarizabilities
	public double BEHp5; //  highest eigenvalue n. 5 of Burden matrix / weighted by atomic polarizabilities
	public double BEHp6; //  highest eigenvalue n. 6 of Burden matrix / weighted by atomic polarizabilities
	public double BEHp7; //  highest eigenvalue n. 7 of Burden matrix / weighted by atomic polarizabilities
	public double BEHp8; //  highest eigenvalue n. 8 of Burden matrix / weighted by atomic polarizabilities
	public double BELp1; //  lowest eigenvalue n. 1 of Burden matrix / weighted by atomic polarizabilities
	public double BELp2; //  lowest eigenvalue n. 2 of Burden matrix / weighted by atomic polarizabilities
	public double BELp3; //  lowest eigenvalue n. 3 of Burden matrix / weighted by atomic polarizabilities
	public double BELp4; //  lowest eigenvalue n. 4 of Burden matrix / weighted by atomic polarizabilities
	public double BELp5; //  lowest eigenvalue n. 5 of Burden matrix / weighted by atomic polarizabilities
	public double BELp6; //  lowest eigenvalue n. 6 of Burden matrix / weighted by atomic polarizabilities
	public double BELp7; //  lowest eigenvalue n. 7 of Burden matrix / weighted by atomic polarizabilities
	public double BELp8; //  lowest eigenvalue n. 8 of Burden matrix / weighted by atomic polarizabilities

	// ***********************************************************************
	// Walk and path counts
	public double MWC01;    // molecular walk count of order 1
	public double MWC02;    // molecular walk count of order 2
	public double MWC03;    // molecular walk count of order 3
	public double MWC04;    // molecular walk count of order 4
	public double MWC05;    // molecular walk count of order 5
	public double MWC06;    // molecular walk count of order 6
	public double MWC07;    // molecular walk count of order 7
	public double MWC08;    // molecular walk count of order 8
	public double MWC09;    // molecular walk count of order 9
	public double MWC10;    // molecular walk count of order 10

	public double SRW01;    // self-returning walk count of order 1
	public double SRW02;    // self-returning walk count of order 2
	public double SRW03;    // self-returning walk count of order 3
	public double SRW04;    // self-returning walk count of order 4
	public double SRW05;    // self-returning walk count of order 5
	public double SRW06;    // self-returning walk count of order 6
	public double SRW07;    // self-returning walk count of order 7
	public double SRW08;    // self-returning walk count of order 8
	public double SRW09;    // self-returning walk count of order 9
	public double SRW10;    // self-returning walk count of order 10

	public double MPC01;    // molecular path count of order 1
	public double MPC02;    // molecular path count of order 2
	public double MPC03;    // molecular path count of order 3
	public double MPC04;    // molecular path count of order 4
	public double MPC05;    // molecular path count of order 5
	public double MPC06;    // molecular path count of order 6
	public double MPC07;    // molecular path count of order 7
	public double MPC08;    // molecular path count of order 8
	public double MPC09;    // molecular path count of order 9
	public double MPC10;    // molecular path count of order 10

	public double piPC01;   // molecular multiple path count of order 01 (sum of conventional bond orders, SCBO)
	public double piPC02;   // molecular multiple path count of order 02
	public double piPC03;   // molecular multiple path count of order 03
	public double piPC04;   // molecular multiple path count of order 04
	public double piPC05;   // molecular multiple path count of order 05
	public double piPC06;   // molecular multiple path count of order 06
	public double piPC07;   // molecular multiple path count of order 07
	public double piPC08;   // molecular multiple path count of order 08
	public double piPC09;   // molecular multiple path count of order 09
	public double piPC10;   // molecular multiple path count of order 10
	public double piID;     // conventional bond order id number

	public double TWC;      // total walk count
	public double TPC;      // total path count
	public double CID;      // Randic Connectivity ID number
	public double CID2;     // Average Randic Connectivity ID number
	public double BID;      // Balaban ID number


	//	 ***********************************************************************
	// 2D Autocorrelation descriptors

	public double ATS1m;  // Broto-Moreau autocorrelation of a topological structure - lag 1 / weighted by atomic masses
	public double ATS2m;  // Broto-Moreau autocorrelation of a topological structure - lag 2 / weighted by atomic masses
	public double ATS3m;  // Broto-Moreau autocorrelation of a topological structure - lag 3 / weighted by atomic masses
	public double ATS4m;  // Broto-Moreau autocorrelation of a topological structure - lag 4 / weighted by atomic masses
	public double ATS5m;  // Broto-Moreau autocorrelation of a topological structure - lag 5 / weighted by atomic masses
	public double ATS6m;  // Broto-Moreau autocorrelation of a topological structure - lag 6 / weighted by atomic masses
	public double ATS7m;  // Broto-Moreau autocorrelation of a topological structure - lag 7 / weighted by atomic masses
	public double ATS8m;  // Broto-Moreau autocorrelation of a topological structure - lag 8 / weighted by atomic masses
	public double ATS1v;  // Broto-Moreau autocorrelation of a topological structure - lag 1 / weighted by atomic van der Waals volumes
	public double ATS2v;  // Broto-Moreau autocorrelation of a topological structure - lag 2 / weighted by atomic van der Waals volumes
	public double ATS3v;  // Broto-Moreau autocorrelation of a topological structure - lag 3 / weighted by atomic van der Waals volumes
	public double ATS4v;  // Broto-Moreau autocorrelation of a topological structure - lag 4 / weighted by atomic van der Waals volumes
	public double ATS5v;  // Broto-Moreau autocorrelation of a topological structure - lag 5 / weighted by atomic van der Waals volumes
	public double ATS6v;  // Broto-Moreau autocorrelation of a topological structure - lag 6 / weighted by atomic van der Waals volumes
	public double ATS7v;  // Broto-Moreau autocorrelation of a topological structure - lag 7 / weighted by atomic van der Waals volumes
	public double ATS8v;  // Broto-Moreau autocorrelation of a topological structure - lag 8 / weighted by atomic van der Waals volumes
	public double ATS1e;  // Broto-Moreau autocorrelation of a topological structure - lag 1 / weighted by atomic Sanderson electronegativities
	public double ATS2e;  // Broto-Moreau autocorrelation of a topological structure - lag 2 / weighted by atomic Sanderson electronegativities
	public double ATS3e;  // Broto-Moreau autocorrelation of a topological structure - lag 3 / weighted by atomic Sanderson electronegativities
	public double ATS4e;  // Broto-Moreau autocorrelation of a topological structure - lag 4 / weighted by atomic Sanderson electronegativities
	public double ATS5e;  // Broto-Moreau autocorrelation of a topological structure - lag 5 / weighted by atomic Sanderson electronegativities
	public double ATS6e;  // Broto-Moreau autocorrelation of a topological structure - lag 6 / weighted by atomic Sanderson electronegativities
	public double ATS7e;  // Broto-Moreau autocorrelation of a topological structure - lag 7 / weighted by atomic Sanderson electronegativities
	public double ATS8e;  // Broto-Moreau autocorrelation of a topological structure - lag 8 / weighted by atomic Sanderson electronegativities
	public double ATS1p;  // Broto-Moreau autocorrelation of a topological structure - lag 1 / weighted by atomic polarizabilities
	public double ATS2p;  // Broto-Moreau autocorrelation of a topological structure - lag 2 / weighted by atomic polarizabilities
	public double ATS3p;  // Broto-Moreau autocorrelation of a topological structure - lag 3 / weighted by atomic polarizabilities
	public double ATS4p;  // Broto-Moreau autocorrelation of a topological structure - lag 4 / weighted by atomic polarizabilities
	public double ATS5p;  // Broto-Moreau autocorrelation of a topological structure - lag 5 / weighted by atomic polarizabilities
	public double ATS6p;  // Broto-Moreau autocorrelation of a topological structure - lag 6 / weighted by atomic polarizabilities
	public double ATS7p;  // Broto-Moreau autocorrelation of a topological structure - lag 7 / weighted by atomic polarizabilities
	public double ATS8p;  // Broto-Moreau autocorrelation of a topological structure - lag 8 / weighted by atomic polarizabilities

	public double MATS1m; // Moran autocorrelation - lag 1 / weighted by atomic masses
	public double MATS2m; // Moran autocorrelation - lag 2 / weighted by atomic masses
	public double MATS3m; // Moran autocorrelation - lag 3 / weighted by atomic masses
	public double MATS4m; // Moran autocorrelation - lag 4 / weighted by atomic masses
	public double MATS5m; // Moran autocorrelation - lag 5 / weighted by atomic masses
	public double MATS6m; // Moran autocorrelation - lag 6 / weighted by atomic masses
	public double MATS7m; // Moran autocorrelation - lag 7 / weighted by atomic masses
	public double MATS8m; // Moran autocorrelation - lag 8 / weighted by atomic masses
	public double MATS1v; // Moran autocorrelation - lag 1 / weighted by atomic van der Waals volumes
	public double MATS2v; // Moran autocorrelation - lag 2 / weighted by atomic van der Waals volumes
	public double MATS3v; // Moran autocorrelation - lag 3 / weighted by atomic van der Waals volumes
	public double MATS4v; // Moran autocorrelation - lag 4 / weighted by atomic van der Waals volumes
	public double MATS5v; // Moran autocorrelation - lag 5 / weighted by atomic van der Waals volumes
	public double MATS6v; // Moran autocorrelation - lag 6 / weighted by atomic van der Waals volumes
	public double MATS7v; // Moran autocorrelation - lag 7 / weighted by atomic van der Waals volumes
	public double MATS8v; // Moran autocorrelation - lag 8 / weighted by atomic van der Waals volumes
	public double MATS1e; // Moran autocorrelation - lag 1 / weighted by atomic Sanderson electronegativities
	public double MATS2e; // Moran autocorrelation - lag 2 / weighted by atomic Sanderson electronegativities
	public double MATS3e; // Moran autocorrelation - lag 3 / weighted by atomic Sanderson electronegativities
	public double MATS4e; // Moran autocorrelation - lag 4 / weighted by atomic Sanderson electronegativities
	public double MATS5e; // Moran autocorrelation - lag 5 / weighted by atomic Sanderson electronegativities
	public double MATS6e; // Moran autocorrelation - lag 6 / weighted by atomic Sanderson electronegativities
	public double MATS7e; // Moran autocorrelation - lag 7 / weighted by atomic Sanderson electronegativities
	public double MATS8e; // Moran autocorrelation - lag 8 / weighted by atomic Sanderson electronegativities
	public double MATS1p; // Moran autocorrelation - lag 1 / weighted by atomic polarizabilities
	public double MATS2p; // Moran autocorrelation - lag 2 / weighted by atomic polarizabilities
	public double MATS3p; // Moran autocorrelation - lag 3 / weighted by atomic polarizabilities
	public double MATS4p; // Moran autocorrelation - lag 4 / weighted by atomic polarizabilities
	public double MATS5p; // Moran autocorrelation - lag 5 / weighted by atomic polarizabilities
	public double MATS6p; // Moran autocorrelation - lag 6 / weighted by atomic polarizabilities
	public double MATS7p; // Moran autocorrelation - lag 7 / weighted by atomic polarizabilities
	public double MATS8p; // Moran autocorrelation - lag 8 / weighted by atomic polarizabilities

	public double GATS1m; // Geary autocorrelation - lag 1 / weighted by atomic masses
	public double GATS2m; // Geary autocorrelation - lag 2 / weighted by atomic masses
	public double GATS3m; // Geary autocorrelation - lag 3 / weighted by atomic masses
	public double GATS4m; // Geary autocorrelation - lag 4 / weighted by atomic masses
	public double GATS5m; // Geary autocorrelation - lag 5 / weighted by atomic masses
	public double GATS6m; // Geary autocorrelation - lag 6 / weighted by atomic masses
	public double GATS7m; // Geary autocorrelation - lag 7 / weighted by atomic masses
	public double GATS8m; // Geary autocorrelation - lag 8 / weighted by atomic masses
	public double GATS1v; // Geary autocorrelation - lag 1 / weighted by atomic van der Waals volumes
	public double GATS2v; // Geary autocorrelation - lag 2 / weighted by atomic van der Waals volumes
	public double GATS3v; // Geary autocorrelation - lag 3 / weighted by atomic van der Waals volumes
	public double GATS4v; // Geary autocorrelation - lag 4 / weighted by atomic van der Waals volumes
	public double GATS5v; // Geary autocorrelation - lag 5 / weighted by atomic van der Waals volumes
	public double GATS6v; // Geary autocorrelation - lag 6 / weighted by atomic van der Waals volumes
	public double GATS7v; // Geary autocorrelation - lag 7 / weighted by atomic van der Waals volumes
	public double GATS8v; // Geary autocorrelation - lag 8 / weighted by atomic van der Waals volumes
	public double GATS1e; // Geary autocorrelation - lag 1 / weighted by atomic Sanderson electronegativities
	public double GATS2e; // Geary autocorrelation - lag 2 / weighted by atomic Sanderson electronegativities
	public double GATS3e; // Geary autocorrelation - lag 3 / weighted by atomic Sanderson electronegativities
	public double GATS4e; // Geary autocorrelation - lag 4 / weighted by atomic Sanderson electronegativities
	public double GATS5e; // Geary autocorrelation - lag 5 / weighted by atomic Sanderson electronegativities
	public double GATS6e; // Geary autocorrelation - lag 6 / weighted by atomic Sanderson electronegativities
	public double GATS7e; // Geary autocorrelation - lag 7 / weighted by atomic Sanderson electronegativities
	public double GATS8e; // Geary autocorrelation - lag 8 / weighted by atomic Sanderson electronegativities
	public double GATS1p; // Geary autocorrelation - lag 1 / weighted by atomic polarizabilities
	public double GATS2p; // Geary autocorrelation - lag 2 / weighted by atomic polarizabilities
	public double GATS3p; // Geary autocorrelation - lag 3 / weighted by atomic polarizabilities
	public double GATS4p; // Geary autocorrelation - lag 4 / weighted by atomic polarizabilities
	public double GATS5p; // Geary autocorrelation - lag 5 / weighted by atomic polarizabilities
	public double GATS6p; // Geary autocorrelation - lag 6 / weighted by atomic polarizabilities
	public double GATS7p; // Geary autocorrelation - lag 7 / weighted by atomic polarizabilities
	public double GATS8p; // Geary autocorrelation - lag 8 / weighted by atomic polarizabilities


	// ***********************************************************************
	// Molecular Properties

	public double ALOGP;  // Ghose-Crippen octanol water coefficient
	public double ALOGP2; // Ghose-Crippen octanol water coefficient squared
	public double MLOGP;  // Moriguchi octanol water partition coefficient
	public double XLOGP;  // Wang octanol water partition coefficient
	public double XLOGP2; // Wang octanol water partition coefficient squared
	public double KLOGP;  // KOWWIN octanol water partition coefficient
	public double KLOGP2; // KOWWIN octanol water partition coefficient squared
	public double AMR;    // Ghose-Crippen molar refractivity
	public double Hy;     // hydrophilic factor
	public double Ui;     // unsaturation index



// *****************************************************************************
	//3-D Molecular Properties (by QM Calculations) (26)
	//A.  General (8)
/*	public double Volume ;  // Volume
	public double SA;       // Surface area
	public double Ovality;  // Ovality
	public double Dipole;   // Dipole
	public double Polar;    // Polarizability
	public double SpcPolar; // Molar or specific polarizability
	public double HOMO;     // Energy of the highest occupied molecular orbital
	public double LUMO;     // Energy of the lowest unoccupied molecular orbital

	//B.  Charge Descriptors (18)
	public double qpmax;    // Maximum positive charge on any atom
	public double qnmax;    // Maximum negative charge on any atom
	public double Qpos;     // Total positive charge
	public double Qneg;     // Total negative charge
	public double Qtot;     // Sum of absolute charges of all atoms
	public double Qmean;    // Charge polarization or mean absolute charge
	public double Q2;       // Total squared atomic charge
	public double SPP;      // Submolecular polarity parameter
	public double TE1;      // Topological electronic descriptor
	public double  TE2;     // Topological electronic descriptor (bond restricted)
	public double PCWTe;    // Partial charge weighted topological electronic descriptor
	public double LDI;      // Local dipole index
	public double et0;      // 0th order electronic topological index, 0ET
	public double et1;      // 1st order electronic topological index, 1ET
	public double et2;      // 2nd order electronic topological index, 2ET
	public double et3;      // 3rd order electronic topological index, 3ET
	public double MaxHp;    // Maximum positive charge on Hydrogen atom
	public double ABSQon;   // Sum of absolute charges of all Oxygen and Nitrogen atoms
*/
	//	 ********************************************************************************
	// XIV.  Principal moment of inertia (7) (see DMOMI in ADAPT, see pg 352 and 405 Todeschini book
	public double MOMI1;  // X-principal component
	public double MOMI2;  // Y-principal component
	public double MOMI3;  // Z-principal component
	public double MOMI4;  // X/Y ratio
	public double MOMI5;  // X/Z ratio
	public double MOMI6;  // Y/Z ratio
	public double MOMI7;  // radius of gyration of


	// ***********************************************************************
	//Charged Partial Surface Area descriptors (29) (pg 52-55 Todeschini book)
	public double PNSA1;    // Partial negative surface area
	public double PPSA1;    // Partial positive surface area
	public double PNSA2;    // Total charge weighted negative surface area
	public double PPSA2;    // Total charge weighted positive surface area
	public double PNSA3;    // Atomic charge weighted negative surface area
	public double PPSA3;    // Atomic charge weighted positive surface area
	public double DPSA1;    // Difference in charged partial surface area
	public double DPSA2;    // Difference in total charge weighted surface area
	public double DPSA3;    // Difference in atomic charge weighted surface area
	public double FNSA1;    // Fractional charged partial negative surface area 1
	public double FNSA2;    // Fractional charged partial negative surface area 2
	public double FNSA3;    // Fractional charged partial negative surface area 3
	public double FPSA1;    // Fractional charged partial positive surface area 1
	public double FFSA2;    // Fractional charged partial positive surface area 2
	public double FPSA3;    // Fractional charged partial positive surface area 3
	public double WNSA1;    // Surface weighted charged partial negative surface area 1
	public double WNSA2;    // Surface weighted charged partial negative surface area 2
	public double WNSA3;    // Surface weighted charged partial negative surface area 3
	public double WPSA1;    // Surface weighted charged partial positive surface area 1
	public double WPSA2;    // Surface weighted charged partial positive surface area 2
	public double WPSA3;    // Surface weighted charged partial positive surface area 3
	public double RNCG;     // Relative negative charge
	public double RPCG;     // Relative positive charge
	public double RNCS;     // Relative negative charge surface area
	public double RPCS;     // Relative positive charge surface area
	public double THSA;     // Total hydrophobic surface area
	public double TPSA;     // Total polar surface area
	public double RHSA;     // Relative hydrophobic surface area
	public double RPSA;     // Relative polar surface area

	// ***********************************************************************
	// Randic molecular profiles
	public double DP01;   // molecular profile no. 01
	public double DP02;   // molecular profile no. 02
	public double DP03;   // molecular profile no. 03
	public double DP04;   // molecular profile no. 04
	public double DP05;   // molecular profile no. 05
	public double DP06;   // molecular profile no. 06
	public double DP07;   // molecular profile no. 07
	public double DP08;   // molecular profile no. 08
	public double DP09;   // molecular profile no. 09
	public double DP10;   // molecular profile no. 10
	public double DP11;   // molecular profile no. 11
	public double DP12;   // molecular profile no. 12
	public double DP13;   // molecular profile no. 13
	public double DP14;   // molecular profile no. 14
	public double DP15;   // molecular profile no. 15
	public double DP16;   // molecular profile no. 16
	public double DP17;   // molecular profile no. 17
	public double DP18;   // molecular profile no. 18
	public double DP19;   // molecular profile no. 19
	public double DP20;   // molecular profile no. 20
	public double SP01;   // shape profile no. 01
	public double SP02;   // shape profile no. 02
	public double SP03;   // shape profile no. 03
	public double SP04;   // shape profile no. 04
	public double SP05;   // shape profile no. 05
	public double SP06;   // shape profile no. 06
	public double SP07;   // shape profile no. 07
	public double SP08;   // shape profile no. 08
	public double SP09;   // shape profile no. 09
	public double SP10;   // shape profile no. 10
	public double SP11;   // shape profile no. 11
	public double SP12;   // shape profile no. 12
	public double SP13;   // shape profile no. 13
	public double SP14;   // shape profile no. 14
	public double SP15;   // shape profile no. 15
	public double SP16;   // shape profile no. 16
	public double SP17;   // shape profile no. 17
	public double SP18;   // shape profile no. 18
	public double SP19;   // shape profile no. 19
	public double SP20;   // shape profile no. 20
	public double SHP2;   // average shape profile index of order 2 (SHP2)

//	 ***********************************************************************
	//X.  WHIM Descriptors
	public double L1u; // 1st component size directional WHIM index / unweighted
	public double L2u; // 2nd component size directional WHIM index / unweighted
	public double L3u; // 3rd component size directional WHIM index / unweighted
	public double P1u; // 1st component shape directional WHIM index / unweighted
	public double P2u; // 2nd component shape directional WHIM index / unweighted
	public double G1u; // 1st component symmetry directional WHIM index / unweighted
	public double G2u; // 2st component symmetry directional WHIM index / unweighted
	public double G3u; // 3st component symmetry directional WHIM index / unweighted
	public double E1u; // 1st component accessibility directional WHIM index / unweighted
	public double E2u; // 2nd component accessibility directional WHIM index / unweighted
	public double E3u; // 3rd component accessibility directional WHIM index / unweighted
	public double L1m; // 1st component size directional WHIM index / weighted by atomic masses
	public double L2m; // 2nd component size directional WHIM index / weighted by atomic masses
	public double L3m; // 3rd component size directional WHIM index / weighted by atomic masses
	public double P1m; // 1st component shape directional WHIM index / weighted by atomic masses
	public double P2m; // 2nd component shape directional WHIM index / weighted by atomic masses
	public double G1m; // 1st component symmetry directional WHIM index / weighted by atomic masses
	public double G2m; // 2st component symmetry directional WHIM index / weighted by atomic masses
	public double G3m; // 3st component symmetry directional WHIM index / weighted by atomic masses
	public double E1m; // 1st component accessibility directional WHIM index / weighted by atomic masses
	public double E2m; // 2nd component accessibility directional WHIM index / weighted by atomic masses
	public double E3m; // 3rd component accessibility directional WHIM index / weighted by atomic masses
	public double L1v; // 1st component size directional WHIM index / weighted by atomic van der Waals volumes
	public double L2v; // 2nd component size directional WHIM index / weighted by atomic van der Waals volumes
	public double L3v; // 3rd component size directional WHIM index / weighted by atomic van der Waals volumes
	public double P1v; // 1st component shape directional WHIM index / weighted by atomic van der Waals volumes
	public double P2v; // 2nd component shape directional WHIM index / weighted by atomic van der Waals volumes
	public double G1v; // 1st component symmetry directional WHIM index / weighted by atomic van der Waals volumes
	public double G2v; // 2st component symmetry directional WHIM index / weighted by atomic van der Waals volumes
	public double G3v; // 3st component symmetry directional WHIM index / weighted by atomic van der Waals volumes
	public double E1v; // 1st component accessibility directional WHIM index / weighted by atomic van der Waals volumes
	public double E2v; // 2nd component accessibility directional WHIM index / weighted by atomic van der Waals volumes
	public double E3v; // 3rd component accessibility directional WHIM index / weighted by atomic van der Waals volumes
	public double L1e; // 1st component size directional WHIM index / weighted by atomic Sanderson electronegativities
	public double L2e; // 2nd component size directional WHIM index / weighted by atomic Sanderson electronegativities
	public double L3e; // 3rd component size directional WHIM index / weighted by atomic Sanderson electronegativities
	public double P1e; // 1st component shape directional WHIM index / weighted by atomic Sanderson electronegativities
	public double P2e; // 2nd component shape directional WHIM index / weighted by atomic Sanderson electronegativities
	public double G1e; // 1st component symmetry directional WHIM index / weighted by atomic Sanderson electronegativities
	public double G2e; // 2st component symmetry directional WHIM index / weighted by atomic Sanderson electronegativities
	public double G3e; // 3st component symmetry directional WHIM index / weighted by atomic Sanderson electronegativities
	public double E1e; // 1st component accessibility directional WHIM index / weighted by atomic Sanderson electronegativities
	public double E2e; // 2nd component accessibility directional WHIM index / weighted by atomic Sanderson electronegativities
	public double E3e; // 3rd component accessibility directional WHIM index / weighted by atomic Sanderson electronegativities
	public double L1p; // 1st component size directional WHIM index / weighted by atomic polarizabilities
	public double L2p; // 2nd component size directional WHIM index / weighted by atomic polarizabilities
	public double L3p; // 3rd component size directional WHIM index / weighted by atomic polarizabilities
	public double P1p; // 1st component shape directional WHIM index / weighted by atomic polarizabilities
	public double P2p; // 2nd component shape directional WHIM index / weighted by atomic polarizabilities
	public double G1p; // 1st component symmetry directional WHIM index / weighted by atomic polarizabilities
	public double G2p; // 2st component symmetry directional WHIM index / weighted by atomic polarizabilities
	public double G3p; // 3st component symmetry directional WHIM index / weighted by atomic polarizabilities
	public double E1p; // 1st component accessibility directional WHIM index / weighted by atomic polarizabilities
	public double E2p; // 2nd component accessibility directional WHIM index / weighted by atomic polarizabilities
	public double E3p; // 3rd component accessibility directional WHIM index / weighted by atomic polarizabilities
	public double L1s; // 1st component size directional WHIM index / weighted by atomic electrotopological states
	public double L2s; // 2nd component size directional WHIM index / weighted by atomic electrotopological states
	public double L3s; // 3rd component size directional WHIM index / weighted by atomic electrotopological states
	public double P1s; // 1st component shape directional WHIM index / weighted by atomic electrotopological states
	public double P2s; // 2nd component shape directional WHIM index / weighted by atomic electrotopological states
	public double G1s; // 1st component symmetry directional WHIM index / weighted by atomic electrotopological states
	public double G2s; // 2st component symmetry directional WHIM index / weighted by atomic electrotopological states
	public double G3s; // 3st component symmetry directional WHIM index / weighted by atomic electrotopological states
	public double E1s; // 1st component accessibility directional WHIM index / weighted by atomic electrotopological states
	public double E2s; // 2nd component accessibility directional WHIM index / weighted by atomic electrotopological states
	public double E3s; // 3rd component accessibility directional WHIM index / weighted by atomic electrotopological states
	public double Tu;  // T total size index / unweighted
	public double Tm;  // T total size index / weighted by atomic masses
	public double Tv;  // T total size index / weighted by atomic van der Waals volumes
	public double Te;  // T total size index / weighted by atomic Sanderson electronegativities
	public double Tp;  // T total size index / weighted by atomic polarizabilities
	public double Ts;  // T total size index / weighted by atomic electrotopological states
	public double Au;  // A total size index / unweighted
	public double Am;  // A total size index / weighted by atomic masses
	public double Av;  // A total size index / weighted by atomic van der Waals volumes
	public double Ae;  // A total size index / weighted by atomic Sanderson electronegativities
	public double Ap;  // A total size index / weighted by atomic polarizabilities
	public double As;  // A total size index / weighted by atomic electrotopological states
	public double Vu;  // V total size index / unweighted
	public double Vm;  // V total size index / weighted by atomic masses
	public double Vv;  // V total size index / weighted by atomic van der Waals volumes
	public double Ve;  // V total size index / weighted by atomic Sanderson electronegativities
	public double Vp;  // V total size index / weighted by atomic polarizabilities
	public double Vs;  // V total size index / weighted by atomic electrotopological states
	public double Gu;  // G total symmetry index / unweighted
	public double Gm;  // G total symmetry index / weighted by atomic masses
	public double Gs;  // G total symmetry index / weighted by atomic electrotopological states
	public double Ku;  // K global shape index / unweighted
	public double Km;  // K global shape index / weighted by atomic masses
	public double Kv;  // K global shape index / weighted by atomic van der Waals volumes
	public double Ke;  // K global shape index / weighted by atomic Sanderson electronegativities
	public double Kp;  // K global shape index / weighted by atomic polarizabilities
	public double Ks;  // K global shape index / weighted by atomic electrotopological states
	public double Du;  // D total accessibility index / unweighted
	public double Dm;  // D total accessibility index / weighted by atomic masses
	public double Dv;  // D total accessibility index / weighted by atomic van der Waals volumes
	public double De;  // D total accessibility index / weighted by atomic Sanderson electronegativities
	public double Dp;  // D total accessibility index / weighted by atomic polarizabilities
	public double Ds;  // D total accessibility index / weighted by atomic electrotopological states

//	 ***********************************************************************
	// XII.  GETAWAY descriptors

	public double HGM;    // geometric mean on the leverage magnitude
	public double ITH;    // total information content on the leverage equality
	public double ISH;    // standardized information content on the leverage equality
	public double HIC;    // mean information content on the leverage magnitude
	public double H0u;    // H autocorrelation of lag 0 / unweighted
	public double H1u;    // H autocorrelation of lag 1 / unweighted
	public double H2u;    // H autocorrelation of lag 2 / unweighted
	public double H3u;    // H autocorrelation of lag 3 / unweighted
	public double H4u;    // H autocorrelation of lag 4 / unweighted
	public double H5u;    // H autocorrelation of lag 5 / unweighted
	public double H6u;    // H autocorrelation of lag 6 / unweighted
	public double H7u;    // H autocorrelation of lag 7 / unweighted
	public double H8u;    // H autocorrelation of lag 8 / unweighted
	public double HTu;    // H total index / unweighted
	public double HATS0u; // leverage-weighted autocorrelation of lag 0 / unweighted
	public double HATS1u; // leverage-weighted autocorrelation of lag 1 / unweighted
	public double HATS2u; // leverage-weighted autocorrelation of lag 2 / unweighted
	public double HATS3u; // leverage-weighted autocorrelation of lag 3 / unweighted
	public double HATS4u; // leverage-weighted autocorrelation of lag 4 / unweighted
	public double HATS5u; // leverage-weighted autocorrelation of lag 5 / unweighted
	public double HATS6u; // leverage-weighted autocorrelation of lag 6 / unweighted
	public double HATS7u; // leverage-weighted autocorrelation of lag 7 / unweighted
	public double HATS8u; // leverage-weighted autocorrelation of lag 8 / unweighted
	public double HATSu;  // leverage-weighted total index / unweighted
	public double H0m;    // H autocorrelation of lag 0 / weighted by atomic masses
	public double H1m;    // H autocorrelation of lag 1 / weighted by atomic masses
	public double H2m;    // H autocorrelation of lag 2 / weighted by atomic masses
	public double H3m;    // H autocorrelation of lag 3 / weighted by atomic masses
	public double H4m;    // H autocorrelation of lag 4 / weighted by atomic masses
	public double H5m;    // H autocorrelation of lag 5 / weighted by atomic masses
	public double H6m;    // H autocorrelation of lag 6 / weighted by atomic masses
	public double H7m;    // H autocorrelation of lag 7 / weighted by atomic masses
	public double H8m;    // H autocorrelation of lag 8 / weighted by atomic masses
	public double HTm;    // H total index / weighted by atomic masses
	public double HATS0m; // leverage-weighted autocorrelation of lag 0 / weighted by atomic masses
	public double HATS1m; // leverage-weighted autocorrelation of lag 1 / weighted by atomic masses
	public double HATS2m; // leverage-weighted autocorrelation of lag 2 / weighted by atomic masses
	public double HATS3m; // leverage-weighted autocorrelation of lag 3 / weighted by atomic masses
	public double HATS4m; // leverage-weighted autocorrelation of lag 4 / weighted by atomic masses
	public double HATS5m; // leverage-weighted autocorrelation of lag 5 / weighted by atomic masses
	public double HATS6m; // leverage-weighted autocorrelation of lag 6 / weighted by atomic masses
	public double HATS7m; // leverage-weighted autocorrelation of lag 7 / weighted by atomic masses
	public double HATS8m; // leverage-weighted autocorrelation of lag 8 / weighted by atomic masses
	public double HATSm;  // leverage-weighted total index / weighted by atomic masses
	public double H0v;    // H autocorrelation of lag 0 / weighted by atomic van der Waals volumes
	public double H1v;    // H autocorrelation of lag 1 / weighted by atomic van der Waals volumes
	public double H2v;    // H autocorrelation of lag 2 / weighted by atomic van der Waals volumes
	public double H3v;    // H autocorrelation of lag 3 / weighted by atomic van der Waals volumes
	public double H4v;    // H autocorrelation of lag 4 / weighted by atomic van der Waals volumes
	public double H5v;    // H autocorrelation of lag 5 / weighted by atomic van der Waals volumes
	public double H6v;    // H autocorrelation of lag 6 / weighted by atomic van der Waals volumes
	public double H7v;    // H autocorrelation of lag 7 / weighted by atomic van der Waals volumes
	public double H8v;    // H autocorrelation of lag 8 / weighted by atomic van der Waals volumes
	public double HTv;    // H total index / weighted by atomic van der Waals volumes
	public double HATS0v; // leverage-weighted autocorrelation of lag 0 / weighted by atomic van der Waals volumes
	public double HATS1v; // leverage-weighted autocorrelation of lag 1 / weighted by atomic van der Waals volumes
	public double HATS2v; // leverage-weighted autocorrelation of lag 2 / weighted by atomic van der Waals volumes
	public double HATS3v; // leverage-weighted autocorrelation of lag 3 / weighted by atomic van der Waals volumes
	public double HATS4v; // leverage-weighted autocorrelation of lag 4 / weighted by atomic van der Waals volumes
	public double HATS5v; // leverage-weighted autocorrelation of lag 5 / weighted by atomic van der Waals volumes
	public double HATS6v; // leverage-weighted autocorrelation of lag 6 / weighted by atomic van der Waals volumes
	public double HATS7v; // leverage-weighted autocorrelation of lag 7 / weighted by atomic van der Waals volumes
	public double HATS8v; // leverage-weighted autocorrelation of lag 8 / weighted by atomic van der Waals volumes
	public double HATSv;  // leverage-weighted total index / weighted by atomic van der Waals volumes
	public double H0e;    // H autocorrelation of lag 0 / weighted by atomic Sanderson electronegativities
	public double H1e;    // H autocorrelation of lag 1 / weighted by atomic Sanderson electronegativities
	public double H2e;    // H autocorrelation of lag 2 / weighted by atomic Sanderson electronegativities
	public double H3e;    // H autocorrelation of lag 3 / weighted by atomic Sanderson electronegativities
	public double H4e;    // H autocorrelation of lag 4 / weighted by atomic Sanderson electronegativities
	public double H5e;    // H autocorrelation of lag 5 / weighted by atomic Sanderson electronegativities
	public double H6e;    // H autocorrelation of lag 6 / weighted by atomic Sanderson electronegativities
	public double H7e;    // H autocorrelation of lag 7 / weighted by atomic Sanderson electronegativities
	public double H8e;    // H autocorrelation of lag 8 / weighted by atomic Sanderson electronegativities
	public double HTe;    // H total index / weighted by atomic Sanderson electronegativities
	public double HATS0e; // leverage-weighted autocorrelation of lag 0 / weighted by atomic Sanderson electronegativities
	public double HATS1e; // leverage-weighted autocorrelation of lag 1 / weighted by atomic Sanderson electronegativities
	public double HATS2e; // leverage-weighted autocorrelation of lag 2 / weighted by atomic Sanderson electronegativities
	public double HATS3e; // leverage-weighted autocorrelation of lag 3 / weighted by atomic Sanderson electronegativities
	public double HATS4e; // leverage-weighted autocorrelation of lag 4 / weighted by atomic Sanderson electronegativities
	public double HATS5e; // leverage-weighted autocorrelation of lag 5 / weighted by atomic Sanderson electronegativities
	public double HATS6e; // leverage-weighted autocorrelation of lag 6 / weighted by atomic Sanderson electronegativities
	public double HATS7e; // leverage-weighted autocorrelation of lag 7 / weighted by atomic Sanderson electronegativities
	public double HATS8e; // leverage-weighted autocorrelation of lag 8 / weighted by atomic Sanderson electronegativities
	public double HATSe;  // leverage-weighted total index / weighted by atomic Sanderson electronegativities
	public double H0p;    // H autocorrelation of lag 0 / weighted by atomic polarizabilities
	public double H1p;    // H autocorrelation of lag 1 / weighted by atomic polarizabilities
	public double H2p;    // H autocorrelation of lag 2 / weighted by atomic polarizabilities
	public double H3p;    // H autocorrelation of lag 3 / weighted by atomic polarizabilities
	public double H4p;    // H autocorrelation of lag 4 / weighted by atomic polarizabilities
	public double H5p;    // H autocorrelation of lag 5 / weighted by atomic polarizabilities
	public double H6p;    // H autocorrelation of lag 6 / weighted by atomic polarizabilities
	public double H7p;    // H autocorrelation of lag 7 / weighted by atomic polarizabilities
	public double H8p;    // H autocorrelation of lag 8 / weighted by atomic polarizabilities
	public double HTp;    // H total index / weighted by atomic polarizabilities
	public double HATS0p; // leverage-weighted autocorrelation of lag 0 / weighted by atomic polarizabilities
	public double HATS1p; // leverage-weighted autocorrelation of lag 1 / weighted by atomic polarizabilities
	public double HATS2p; // leverage-weighted autocorrelation of lag 2 / weighted by atomic polarizabilities
	public double HATS3p; // leverage-weighted autocorrelation of lag 3 / weighted by atomic polarizabilities
	public double HATS4p; // leverage-weighted autocorrelation of lag 4 / weighted by atomic polarizabilities
	public double HATS5p; // leverage-weighted autocorrelation of lag 5 / weighted by atomic polarizabilities
	public double HATS6p; // leverage-weighted autocorrelation of lag 6 / weighted by atomic polarizabilities
	public double HATS7p; // leverage-weighted autocorrelation of lag 7 / weighted by atomic polarizabilities
	public double HATS8p; // leverage-weighted autocorrelation of lag 8 / weighted by atomic polarizabilities
	public double HATSp;  // leverage-weighted total index / weighted by atomic polarizabilities
	public double RCON;   // Randic-type R matrix connectivity
	public double RARS;   // R matrix average row sum
	public double REIG;   // first eigenvalue of the R matrix
	public double R1u;    // R autocorrelation of lag 1 / unweighted
	public double R2u;    // R autocorrelation of lag 2 / unweighted
	public double R3u;    // R autocorrelation of lag 3 / unweighted
	public double R4u;    // R autocorrelation of lag 4 / unweighted
	public double R5u;    // R autocorrelation of lag 5 / unweighted
	public double R6u;    // R autocorrelation of lag 6 / unweighted
	public double R7u;    // R autocorrelation of lag 7 / unweighted
	public double R8u;    // R autocorrelation of lag 8 / unweighted
	public double RTu;    // R total index / unweighted
	public double R1uplus;// R maximal autocorrelation of lag 1 / unweighted
	public double R2uplus;// R maximal autocorrelation of lag 2 / unweighted
	public double R3uplus;// R maximal autocorrelation of lag 3 / unweighted
	public double R4uplus;// R maximal autocorrelation of lag 4 / unweighted
	public double R5uplus;// R maximal autocorrelation of lag 5 / unweighted
	public double R6uplus;// R maximal autocorrelation of lag 6 / unweighted
	public double R7uplus;// R maximal autocorrelation of lag 7 / unweighted
	public double R8uplus;// R maximal autocorrelation of lag 8 / unweighted
	public double RTuplus;// R maximal index / unweighted
	public double R1m;    // R autocorrelation of lag 1 / weighted by atomic masses
	public double R2m;    // R autocorrelation of lag 2 / weighted by atomic masses
	public double R3m;    // R autocorrelation of lag 3 / weighted by atomic masses
	public double R4m;    // R autocorrelation of lag 4 / weighted by atomic masses
	public double R5m;    // R autocorrelation of lag 5 / weighted by atomic masses
	public double R6m;    // R autocorrelation of lag 6 / weighted by atomic masses
	public double R7m;    // R autocorrelation of lag 7 / weighted by atomic masses
	public double R8m;    // R autocorrelation of lag 8 / weighted by atomic masses
	public double RTm;    // R total index / weighted by atomic masses
	public double R1mplus;// R maximal autocorrelation of lag 1 / weighted by atomic masses
	public double R2mplus;// R maximal autocorrelation of lag 2 / weighted by atomic masses
	public double R3mplus;// R maximal autocorrelation of lag 3 / weighted by atomic masses
	public double R4mplus;// R maximal autocorrelation of lag 4 / weighted by atomic masses
	public double R5mplus;// R maximal autocorrelation of lag 5 / weighted by atomic masses
	public double R6mplus;// R maximal autocorrelation of lag 6 / weighted by atomic masses
	public double R7mplus;// R maximal autocorrelation of lag 7 / weighted by atomic masses
	public double R8mplus;// R maximal autocorrelation of lag 8 / weighted by atomic masses
	public double RTmplus;// R maximal index / weighted by atomic masses
	public double R1v;    // R autocorrelation of lag 1 / weighted by atomic van der Waals volumes
	public double R2v;    // R autocorrelation of lag 2 / weighted by atomic van der Waals volumes
	public double R3v;    // R autocorrelation of lag 3 / weighted by atomic van der Waals volumes
	public double R4v;    // R autocorrelation of lag 4 / weighted by atomic van der Waals volumes
	public double R5v;    // R autocorrelation of lag 5 / weighted by atomic van der Waals volumes
	public double R6v;    // R autocorrelation of lag 6 / weighted by atomic van der Waals volumes
	public double R7v;    // R autocorrelation of lag 7 / weighted by atomic van der Waals volumes
	public double R8v;    // R autocorrelation of lag 8 / weighted by atomic van der Waals volumes
	public double RTv;    // R total index / weighted by atomic van der Waals volumes
	public double R1vplus;// R maximal autocorrelation of lag 1 / weighted by atomic van der Waals volumes
	public double R2vplus;// R maximal autocorrelation of lag 2 / weighted by atomic van der Waals volumes
	public double R3vplus;// R maximal autocorrelation of lag 3 / weighted by atomic van der Waals volumes
	public double R4vplus;// R maximal autocorrelation of lag 4 / weighted by atomic van der Waals volumes
	public double R5vplus;// R maximal autocorrelation of lag 5 / weighted by atomic van der Waals volumes
	public double R6vplus;// R maximal autocorrelation of lag 6 / weighted by atomic van der Waals volumes
	public double R7vplus;// R maximal autocorrelation of lag 7 / weighted by atomic van der Waals volumes
	public double R8vplus;// R maximal autocorrelation of lag 8 / weighted by atomic van der Waals volumes
	public double RTvplus;// R maximal index / weighted by atomic van der Waals volumes
	public double R1e;    // R autocorrelation of lag 1 / weighted by atomic Sanderson electronegativities
	public double R2e;    // R autocorrelation of lag 2 / weighted by atomic Sanderson electronegativities
	public double R3e;    // R autocorrelation of lag 3 / weighted by atomic Sanderson electronegativities
	public double R4e;    // R autocorrelation of lag 4 / weighted by atomic Sanderson electronegativities
	public double R5e;    // R autocorrelation of lag 5 / weighted by atomic Sanderson electronegativities
	public double R6e;    // R autocorrelation of lag 6 / weighted by atomic Sanderson electronegativities
	public double R7e;    // R autocorrelation of lag 7 / weighted by atomic Sanderson electronegativities
	public double R8e;    // R autocorrelation of lag 8 / weighted by atomic Sanderson electronegativities
	public double RTe;    // R total index / weighted by atomic Sanderson electronegativities
	public double R1eplus;// R maximal autocorrelation of lag 1 / weighted by atomic Sanderson electronegativities
	public double R2eplus;// R maximal autocorrelation of lag 2 / weighted by atomic Sanderson electronegativities
	public double R3eplus;// R maximal autocorrelation of lag 3 / weighted by atomic Sanderson electronegativities
	public double R4eplus;// R maximal autocorrelation of lag 4 / weighted by atomic Sanderson electronegativities
	public double R5eplus;// R maximal autocorrelation of lag 5 / weighted by atomic Sanderson electronegativities
	public double R6eplus;// R maximal autocorrelation of lag 6 / weighted by atomic Sanderson electronegativities
	public double R7eplus;// R maximal autocorrelation of lag 7 / weighted by atomic Sanderson electronegativities
	public double R8eplus;// R maximal autocorrelation of lag 8 / weighted by atomic Sanderson electronegativities
	public double RTeplus;// R maximal index / weighted by atomic Sanderson electronegativities
	public double R1p;    // R autocorrelation of lag 1 / weighted by atomic polarizabilities
	public double R2p;    // R autocorrelation of lag 2 / weighted by atomic polarizabilities
	public double R3p;    // R autocorrelation of lag 3 / weighted by atomic polarizabilities
	public double R4p;    // R autocorrelation of lag 4 / weighted by atomic polarizabilities
	public double R5p;    // R autocorrelation of lag 5 / weighted by atomic polarizabilities
	public double R6p;    // R autocorrelation of lag 6 / weighted by atomic polarizabilities
	public double R7p;    // R autocorrelation of lag 7 / weighted by atomic polarizabilities
	public double R8p;    // R autocorrelation of lag 8 / weighted by atomic polarizabilities
	public double RTp;    // R total index / weighted by atomic polarizabilities
	public double R1pplus;// R maximal autocorrelation of lag 1 / weighted by atomic polarizabilities
	public double R2pplus;// R maximal autocorrelation of lag 2 / weighted by atomic polarizabilities
	public double R3pplus;// R maximal autocorrelation of lag 3 / weighted by atomic polarizabilities
	public double R4pplus;// R maximal autocorrelation of lag 4 / weighted by atomic polarizabilities
	public double R5pplus;// R maximal autocorrelation of lag 5 / weighted by atomic polarizabilities
	public double R6pplus;// R maximal autocorrelation of lag 6 / weighted by atomic polarizabilities
	public double R7pplus;// R maximal autocorrelation of lag 7 / weighted by atomic polarizabilities
	public double R8pplus;// R maximal autocorrelation of lag 8 / weighted by atomic polarizabilities
	public double RTpplus;// R maximal index / weighted by atomic polarizabilities


	// ***********************************************************************
	//Molecular Orbital (QM) based 3D descriptors

	//Charge based descriptors
	public double qpmax;	// Maximumm positive charge on any atom
	public double qnmax;	//Maximum negative charge on any atom
	public double Qpos;		//Total positive charge
	public double Qneg;		//Total negative charge
	public double Qtot;		//Total absolute atomic charge
	public double Qmean;	//Charge polarization
	public double Q2;		//Total sqared atomic charge
	public double SPP;		//Submolecular polarity parameter - Order 1
	public double SPP2;		//Submolecular polarity parameter - Order 2
	public double Pprime;	//Submolecular polarity parameter/interatomic distance
	public double DP;		//Submolecular polarity parameter/square of interatomic distance
	public double TE1;		//Topological electronic descriptor (between all pairs of atoms)
	public double TE2;		//Topological electronic descriptor (between connected atoms)
	public double PCWTE;	//Partial charge weighted topological electronic index
	public double LDI;		//Local dipole index
	public double ET0;		//Electro-topological index - Order 0
	public double ET1;		//Electro-topological index - Order 1
	public double ET2;		//Electro-topological index - Order 2
	public double ET3;		//Electro-topological index - Order 3

	//Max-min partial charges of atoms
	public double MaxPC_C;	//Maximum partial charge for a carbon atom
	public double MinPC_C;	//Minimum partial charge for a carbon atom
	public double MaxPC_N;	//Maximum partial charge for a nitrogen atom
	public double MinPC_N;	//Minimum partial charge for a nitrogen atom
	public double MaxPC_P;	//Maximum partial charge for a phosphorous atom
	public double MinPC_P;	//Minimum partial charge for a phosphorous atom
	public double MaxPC_O;	//Maximum partial charge for an oxygen atom
	public double MinPC_O;	//Minimum partial charge for an oxygen atom
	public double MaxPC_S;	//Maximum partial charge for a sulfur atom
	public double MinPC_S;	//Minimum partial charge for a sulfur atom
	public double MaxPC_Si;	//Maximum partial charge for a silicon atom
	public double MinPC_Si;	//Minimum partial charge for a silicon atom
	public double MaxPC_Hg;	//Maximum partial charge for a mercury atom
	public double MinPC_Hg;	//Minimum partial charge for a mercury atom
	public double MaxPC_H;	//Maximum partial charge for a hydrogen atom
	public double MinPC_H;	//Minimum partial charge for a hydrogen atom
	public double MaxPC_F;	//Maximum partial charge for a fluorine atom
	public double MinPC_F;	//Minimum partial charge for a fluorine atom
	public double MaxPC_Cl;	//Maximum partial charge for a chlorine atom
	public double MinPC_Cl;	//Minimum partial charge for a chlorine atom
	public double MaxPC_Br;	//Maximum partial charge for a bromine atom
	public double MinPC_Br;	//Minimum partial charge for a bromine atom
	public double MaxPC_I;	//Maximum partial charge for an iodine atom
	public double MinPC_I;	//Minimum partial charge for an iodine atom

	//Electrophilic frontier density (EFD) of atoms
	public double Max_EFD_C;	//maximum EFD of a carbon atom
	public double Min_EFD_C;	//Minimum EFD of a carbon atom
	public double Tot_EFD_C;	//Total EFD of a carbon atom
	public double Av_EFD_C;		//Average EFD of a carbon atom
	public double Max_EFD_N;	//maximum EFD of a nitrogen atom
	public double Min_EFD_N;	//Minimum EFD of a nitrogen atom
	public double Tot_EFD_N;	//Total EFD of a nitrogen atom
	public double Av_EFD_N;		//Average EFD of a nitrogen atom
	public double Max_EFD_P;	//maximum EFD of a phosphorous atom
	public double Min_EFD_P;	//Minimum EFD of a phosphorous atom
	public double Tot_EFD_P;	//Total EFD of a phosphorous atom
	public double Av_EFD_P;		//Average EFD of a phsophorous atom
	public double Max_EFD_O;	//maximum EFD of an oxygen atom
	public double Min_EFD_O;	//Minimum EFD of an oxygen atom
	public double Tot_EFD_O;	//Total EFD of an oxygen atom
	public double Av_EFD_O;		//Average EFD of an oxygen atom
	public double Max_EFD_S;	//maximum EFD of a sulfur atom
	public double Min_EFD_S;	//Minimum EFD of a sulfur atom
	public double Tot_EFD_S;	//Total EFD of a sulfur atom
	public double Av_EFD_S;		//Average EFD of a sulfur atom
	public double Max_EFD_Si;	//maximum EFD of a silicon atom
	public double Min_EFD_Si;	//Minimum EFD of a silicon atom
	public double Tot_EFD_Si;	//Total EFD of a silicon atom
	public double Av_EFD_Si;	//Average EFD of a silicon atom
	public double Max_EFD_Hg;	//maximum EFD of a mercury atom
	public double Min_EFD_Hg;	//Minimum EFD of a mercury atom
	public double Tot_EFD_Hg;	//Total EFD of a mercury atom
	public double Av_EFD_Hg;	//Average EFD of a mercury atom
	public double Max_EFD_H;	//maximum EFD of a hydrogen atom
	public double Min_EFD_H;	//Minimum EFD of a hydrogen atom
	public double Tot_EFD_H;	//Total EFD of a hydrogen atom
	public double Av_EFD_H;		//Average EFD of a hydrogen atom
	public double Max_EFD_F;	//maximum EFD of a fluorine atom
	public double Min_EFD_F;	//Minimum EFD of a fluorine atom
	public double Tot_EFD_F;	//Total EFD of a fluorine atom
	public double Av_EFD_F;		//Average EFD of a fluorine atom
	public double Max_EFD_Cl;	//maximum EFD of a chlorine atom
	public double Min_EFD_Cl;	//Minimum EFD of a chlorine atom
	public double Tot_EFD_Cl;	//Total EFD of a chlorine atom
	public double Av_EFD_Cl;	//Average EFD of a chlorine atom
	public double Max_EFD_Br;	//maximum EFD of a bromine atom
	public double Min_EFD_Br;	//Minimum EFD of a bromine atom
	public double Tot_EFD_Br;	//Total EFD of a bromine atom
	public double Av_EFD_Br;	//Average EFD of a bromine atom
	public double Max_EFD_I;	//maximum EFD of an iodine atom
	public double Min_EFD_I;	//Minimum EFD of an iodine atom
	public double Tot_EFD_I;	//Total EFD of an iodine atom
	public double Av_EFD_I;		//Average EFD of an iodine atom
	public double Max_EFD;		//Maximum EFD of all atoms
	public double Min_EFD;		//Minimum EFD of all atoms
	public double Tot_EFD;		//Total EFD of all atoms
	public double Av_EFD;		//Average EFD of all atoms

	//Nucleophilic frontier density (NFD) of atoms
	public double Max_NFD_C;	//maximum NFD of a carbon atom
	public double Min_NFD_C;	//Minimum NFD of a carbon atom
	public double Tot_NFD_C;	//Total NFD of a carbon atom
	public double Av_NFD_C;		//Average NFD of a carbon atom
	public double Max_NFD_N;	//maximum NFD of a nitrogen atom
	public double Min_NFD_N;	//Minimum NFD of a nitrogen atom
	public double Tot_NFD_N;	//Total NFD of a nitrogen atom
	public double Av_NFD_N;		//Average NFD of a nitrogen atom
	public double Max_NFD_P;	//maximum NFD of a phosphorous atom
	public double Min_NFD_P;	//Minimum NFD of a phosphorous atom
	public double Tot_NFD_P;	//Total NFD of a phosphorous atom
	public double Av_NFD_P;		//Average NFD of a phsophorous atom
	public double Max_NFD_O;	//maximum NFD of an oxygen atom
	public double Min_NFD_O;	//Minimum NFD of an oxygen atom
	public double Tot_NFD_O;	//Total NFD of an oxygen atom
	public double Av_NFD_O;		//Average NFD of an oxygen atom
	public double Max_NFD_S;	//maximum NFD of a sulfur atom
	public double Min_NFD_S;	//Minimum NFD of a sulfur atom
	public double Tot_NFD_S;	//Total NFD of a sulfur atom
	public double Av_NFD_S;		//Average NFD of a sulfur atom
	public double Max_NFD_Si;	//maximum NFD of a silicon atom
	public double Min_NFD_Si;	//Minimum NFD of a silicon atom
	public double Tot_NFD_Si;	//Total NFD of a silicon atom
	public double Av_NFD_Si;	//Average NFD of a silicon atom
	public double Max_NFD_Hg;	//maximum NFD of a mercury atom
	public double Min_NFD_Hg;	//Minimum NFD of a mercury atom
	public double Tot_NFD_Hg;	//Total NFD of a mercury atom
	public double Av_NFD_Hg;	//Average NFD of a mercury atom
	public double Max_NFD_H;	//maximum NFD of a hydrogen atom
	public double Min_NFD_H;	//Minimum NFD of a hydrogen atom
	public double Tot_NFD_H;	//Total NFD of a hydrogen atom
	public double Av_NFD_H;		//Average NFD of a hydrogen atom
	public double Max_NFD_F;	//maximum NFD of a fluorine atom
	public double Min_NFD_F;	//Minimum NFD of a fluorine atom
	public double Tot_NFD_F;	//Total NFD of a fluorine atom
	public double Av_NFD_F;		//Average NFD of a fluorine atom
	public double Max_NFD_Cl;	//maximum NFD of a chlorine atom
	public double Min_NFD_Cl;	//Minimum NFD of a chlorine atom
	public double Tot_NFD_Cl;	//Total NFD of a chlorine atom
	public double Av_NFD_Cl;	//Average NFD of a chlorine atom
	public double Max_NFD_Br;	//maximum NFD of a bromine atom
	public double Min_NFD_Br;	//Minimum NFD of a bromine atom
	public double Tot_NFD_Br;	//Total NFD of a bromine atom
	public double Av_NFD_Br;	//Average NFD of a bromine atom
	public double Max_NFD_I;	//maximum NFD of an iodine atom
	public double Min_NFD_I;	//Minimum NFD of an iodine atom
	public double Tot_NFD_I;	//Total NFD of an iodine atom
	public double Av_NFD_I;		//Average NFD of an iodine atom
	public double Max_NFD;		//Maximum NFD of all atoms
	public double Min_NFD;		//Minimum NFD of all atoms
	public double Tot_NFD;		//Total NFD of all atoms
	public double Av_NFD;		//Average NFD of all atoms

	//Radical Frontier Density (RFD) of atoms
	public double Max_RFD_C;	//maximum RFD of a carbon atom
	public double Min_RFD_C;	//Minimum RFD of a carbon atom
	public double Tot_RFD_C;	//Total RFD of a carbon atom
	public double Av_RFD_C;		//Average RFD of a carbon atom
	public double Max_RFD_N;	//maximum RFD of a nitrogen atom
	public double Min_RFD_N;	//Minimum RFD of a nitrogen atom
	public double Tot_RFD_N;	//Total RFD of a nitrogen atom
	public double Av_RFD_N;		//Average RFD of a nitrogen atom
	public double Max_RFD_P;	//maximum RFD of a phosphorous atom
	public double Min_RFD_P;	//Minimum RFD of a phosphorous atom
	public double Tot_RFD_P;	//Total RFD of a phosphorous atom
	public double Av_RFD_P;		//Average RFD of a phsophorous atom
	public double Max_RFD_O;	//maximum RFD of an oxygen atom
	public double Min_RFD_O;	//Minimum RFD of an oxygen atom
	public double Tot_RFD_O;	//Total RFD of an oxygen atom
	public double Av_RFD_O;		//Average RFD of an oxygen atom
	public double Max_RFD_S;	//maximum RFD of a sulfur atom
	public double Min_RFD_S;	//Minimum RFD of a sulfur atom
	public double Tot_RFD_S;	//Total RFD of a sulfur atom
	public double Av_RFD_S;		//Average RFD of a sulfur atom
	public double Max_RFD_Si;	//maximum RFD of a silicon atom
	public double Min_RFD_Si;	//Minimum RFD of a silicon atom
	public double Tot_RFD_Si;	//Total RFD of a silicon atom
	public double Av_RFD_Si;	//Average RFD of a silicon atom
	public double Max_RFD_Hg;	//maximum RFD of a mercury atom
	public double Min_RFD_Hg;	//Minimum RFD of a mercury atom
	public double Tot_RFD_Hg;	//Total RFD of a mercury atom
	public double Av_RFD_Hg;	//Average RFD of a mercury atom
	public double Max_RFD_H;	//maximum RFD of a hydrogen atom
	public double Min_RFD_H;	//Minimum RFD of a hydrogen atom
	public double Tot_RFD_H;	//Total RFD of a hydrogen atom
	public double Av_RFD_H;		//Average RFD of a hydrogen atom
	public double Max_RFD_F;	//maximum RFD of a fluorine atom
	public double Min_RFD_F;	//Minimum RFD of a fluorine atom
	public double Tot_RFD_F;	//Total RFD of a fluorine atom
	public double Av_RFD_F;		//Average RFD of a fluorine atom
	public double Max_RFD_Cl;	//maximum RFD of a chlorine atom
	public double Min_RFD_Cl;	//Minimum RFD of a chlorine atom
	public double Tot_RFD_Cl;	//Total RFD of a chlorine atom
	public double Av_RFD_Cl;	//Average RFD of a chlorine atom
	public double Max_RFD_Br;	//maximum RFD of a bromine atom
	public double Min_RFD_Br;	//Minimum RFD of a bromine atom
	public double Tot_RFD_Br;	//Total RFD of a bromine atom
	public double Av_RFD_Br;	//Average RFD of a bromine atom
	public double Max_RFD_I;	//maximum RFD of an iodine atom
	public double Min_RFD_I;	//Minimum RFD of an iodine atom
	public double Tot_RFD_I;	//Total RFD of an iodine atom
	public double Av_RFD_I;		//Average RFD of an iodine atom
	public double Max_RFD;		//Maximum RFD of all atoms
	public double Min_RFD;		//Minimum RFD of all atoms
	public double Tot_RFD;		//Total RFD of all atoms
	public double Av_RFD;		//Average RFD of all atoms

	//Electrophilic Frontier Density Index (EFDI) of atoms
	public double Max_EFDI_C;	//maximum EFDI of a carbon atom
	public double Min_EFDI_C;	//Minimum EFDI of a carbon atom
	public double Tot_EFDI_C;	//Total EFDI of a carbon atom
	public double Av_EFDI_C;		//Average EFDI of a carbon atom
	public double Max_EFDI_N;	//maximum EFDI of a nitrogen atom
	public double Min_EFDI_N;	//Minimum EFDI of a nitrogen atom
	public double Tot_EFDI_N;	//Total EFDI of a nitrogen atom
	public double Av_EFDI_N;		//Average EFDI of a nitrogen atom
	public double Max_EFDI_P;	//maximum EFDI of a phosphorous atom
	public double Min_EFDI_P;	//Minimum EFDI of a phosphorous atom
	public double Tot_EFDI_P;	//Total EFDI of a phosphorous atom
	public double Av_EFDI_P;		//Average EFDI of a phsophorous atom
	public double Max_EFDI_O;	//maximum EFDI of an oxygen atom
	public double Min_EFDI_O;	//Minimum EFDI of an oxygen atom
	public double Tot_EFDI_O;	//Total EFDI of an oxygen atom
	public double Av_EFDI_O;		//Average EFDI of an oxygen atom
	public double Max_EFDI_S;	//maximum EFDI of a sulfur atom
	public double Min_EFDI_S;	//Minimum EFDI of a sulfur atom
	public double Tot_EFDI_S;	//Total EFDI of a sulfur atom
	public double Av_EFDI_S;		//Average EFDI of a sulfur atom
	public double Max_EFDI_Si;	//maximum EFDI of a silicon atom
	public double Min_EFDI_Si;	//Minimum EFDI of a silicon atom
	public double Tot_EFDI_Si;	//Total EFDI of a silicon atom
	public double Av_EFDI_Si;	//Average EFDI of a silicon atom
	public double Max_EFDI_Hg;	//maximum EFDI of a mercury atom
	public double Min_EFDI_Hg;	//Minimum EFDI of a mercury atom
	public double Tot_EFDI_Hg;	//Total EFDI of a mercury atom
	public double Av_EFDI_Hg;	//Average EFDI of a mercury atom
	public double Max_EFDI_H;	//maximum EFDI of a hydrogen atom
	public double Min_EFDI_H;	//Minimum EFDI of a hydrogen atom
	public double Tot_EFDI_H;	//Total EFDI of a hydrogen atom
	public double Av_EFDI_H;		//Average EFDI of a hydrogen atom
	public double Max_EFDI_F;	//maximum EFDI of a fluorine atom
	public double Min_EFDI_F;	//Minimum EFDI of a fluorine atom
	public double Tot_EFDI_F;	//Total EFDI of a fluorine atom
	public double Av_EFDI_F;		//Average EFDI of a fluorine atom
	public double Max_EFDI_Cl;	//maximum EFDI of a chlorine atom
	public double Min_EFDI_Cl;	//Minimum EFDI of a chlorine atom
	public double Tot_EFDI_Cl;	//Total EFDI of a chlorine atom
	public double Av_EFDI_Cl;	//Average EFDI of a chlorine atom
	public double Max_EFDI_Br;	//maximum EFDI of a bromine atom
	public double Min_EFDI_Br;	//Minimum EFDI of a bromine atom
	public double Tot_EFDI_Br;	//Total EFDI of a bromine atom
	public double Av_EFDI_Br;	//Average EFDI of a bromine atom
	public double Max_EFDI_I;	//maximum EFDI of an iodine atom
	public double Min_EFDI_I;	//Minimum EFDI of an iodine atom
	public double Tot_EFDI_I;	//Total EFDI of an iodine atom
	public double Av_EFDI_I;		//Average EFDI of an iodine atom
	public double Max_EFDI;		//Maximum EFDI of all atoms
	public double Min_EFDI;		//Minimum EFDI of all atoms
	public double Tot_EFDI;		//Total EFDI of all atoms
	public double Av_EFDI;		//Average EFDI of all atoms

	//Nucleophilic Frontier Density Index (NFDI) of atoms
	public double Max_NFDI_C;	//maximum NFDI of a carbon atom
	public double Min_NFDI_C;	//Minimum NFDI of a carbon atom
	public double Tot_NFDI_C;	//Total NFDI of a carbon atom
	public double Av_NFDI_C;		//Average NFDI of a carbon atom
	public double Max_NFDI_N;	//maximum NFDI of a nitrogen atom
	public double Min_NFDI_N;	//Minimum NFDI of a nitrogen atom
	public double Tot_NFDI_N;	//Total NFDI of a nitrogen atom
	public double Av_NFDI_N;		//Average NFDI of a nitrogen atom
	public double Max_NFDI_P;	//maximum NFDI of a phosphorous atom
	public double Min_NFDI_P;	//Minimum NFDI of a phosphorous atom
	public double Tot_NFDI_P;	//Total NFDI of a phosphorous atom
	public double Av_NFDI_P;		//Average NFDI of a phsophorous atom
	public double Max_NFDI_O;	//maximum NFDI of an oxygen atom
	public double Min_NFDI_O;	//Minimum NFDI of an oxygen atom
	public double Tot_NFDI_O;	//Total NFDI of an oxygen atom
	public double Av_NFDI_O;		//Average NFDI of an oxygen atom
	public double Max_NFDI_S;	//maximum NFDI of a sulfur atom
	public double Min_NFDI_S;	//Minimum NFDI of a sulfur atom
	public double Tot_NFDI_S;	//Total NFDI of a sulfur atom
	public double Av_NFDI_S;		//Average NFDI of a sulfur atom
	public double Max_NFDI_Si;	//maximum NFDI of a silicon atom
	public double Min_NFDI_Si;	//Minimum NFDI of a silicon atom
	public double Tot_NFDI_Si;	//Total NFDI of a silicon atom
	public double Av_NFDI_Si;	//Average NFDI of a silicon atom
	public double Max_NFDI_Hg;	//maximum NFDI of a mercury atom
	public double Min_NFDI_Hg;	//Minimum NFDI of a mercury atom
	public double Tot_NFDI_Hg;	//Total NFDI of a mercury atom
	public double Av_NFDI_Hg;	//Average NFDI of a mercury atom
	public double Max_NFDI_H;	//maximum NFDI of a hydrogen atom
	public double Min_NFDI_H;	//Minimum NFDI of a hydrogen atom
	public double Tot_NFDI_H;	//Total NFDI of a hydrogen atom
	public double Av_NFDI_H;		//Average NFDI of a hydrogen atom
	public double Max_NFDI_F;	//maximum NFDI of a fluorine atom
	public double Min_NFDI_F;	//Minimum NFDI of a fluorine atom
	public double Tot_NFDI_F;	//Total NFDI of a fluorine atom
	public double Av_NFDI_F;		//Average NFDI of a fluorine atom
	public double Max_NFDI_Cl;	//maximum NFDI of a chlorine atom
	public double Min_NFDI_Cl;	//Minimum NFDI of a chlorine atom
	public double Tot_NFDI_Cl;	//Total NFDI of a chlorine atom
	public double Av_NFDI_Cl;	//Average NFDI of a chlorine atom
	public double Max_NFDI_Br;	//maximum NFDI of a bromine atom
	public double Min_NFDI_Br;	//Minimum NFDI of a bromine atom
	public double Tot_NFDI_Br;	//Total NFDI of a bromine atom
	public double Av_NFDI_Br;	//Average NFDI of a bromine atom
	public double Max_NFDI_I;	//maximum NFDI of an iodine atom
	public double Min_NFDI_I;	//Minimum NFDI of an iodine atom
	public double Tot_NFDI_I;	//Total NFDI of an iodine atom
	public double Av_NFDI_I;		//Average NFDI of an iodine atom
	public double Max_NFDI;		//Maximum NFDI of all atoms
	public double Min_NFDI;		//Minimum NFDI of all atoms
	public double Tot_NFDI;		//Total NFDI of all atoms
	public double Av_NFDI;		//Average NFDI of all atoms

	//Radical Frontier Density Index (RFDI) of atoms
	public double Max_RFDI_C;	//maximum RFDI of a carbon atom
	public double Min_RFDI_C;	//Minimum RFDI of a carbon atom
	public double Tot_RFDI_C;	//Total RFDI of a carbon atom
	public double Av_RFDI_C;		//Average RFDI of a carbon atom
	public double Max_RFDI_N;	//maximum RFDI of a nitrogen atom
	public double Min_RFDI_N;	//Minimum RFDI of a nitrogen atom
	public double Tot_RFDI_N;	//Total RFDI of a nitrogen atom
	public double Av_RFDI_N;		//Average RFDI of a nitrogen atom
	public double Max_RFDI_P;	//maximum RFDI of a phosphorous atom
	public double Min_RFDI_P;	//Minimum RFDI of a phosphorous atom
	public double Tot_RFDI_P;	//Total RFDI of a phosphorous atom
	public double Av_RFDI_P;		//Average RFDI of a phsophorous atom
	public double Max_RFDI_O;	//maximum RFDI of an oxygen atom
	public double Min_RFDI_O;	//Minimum RFDI of an oxygen atom
	public double Tot_RFDI_O;	//Total RFDI of an oxygen atom
	public double Av_RFDI_O;		//Average RFDI of an oxygen atom
	public double Max_RFDI_S;	//maximum RFDI of a sulfur atom
	public double Min_RFDI_S;	//Minimum RFDI of a sulfur atom
	public double Tot_RFDI_S;	//Total RFDI of a sulfur atom
	public double Av_RFDI_S;		//Average RFDI of a sulfur atom
	public double Max_RFDI_Si;	//maximum RFDI of a silicon atom
	public double Min_RFDI_Si;	//Minimum RFDI of a silicon atom
	public double Tot_RFDI_Si;	//Total RFDI of a silicon atom
	public double Av_RFDI_Si;	//Average RFDI of a silicon atom
	public double Max_RFDI_Hg;	//maximum RFDI of a mercury atom
	public double Min_RFDI_Hg;	//Minimum RFDI of a mercury atom
	public double Tot_RFDI_Hg;	//Total RFDI of a mercury atom
	public double Av_RFDI_Hg;	//Average RFDI of a mercury atom
	public double Max_RFDI_H;	//maximum RFDI of a hydrogen atom
	public double Min_RFDI_H;	//Minimum RFDI of a hydrogen atom
	public double Tot_RFDI_H;	//Total RFDI of a hydrogen atom
	public double Av_RFDI_H;		//Average RFDI of a hydrogen atom
	public double Max_RFDI_F;	//maximum RFDI of a fluorine atom
	public double Min_RFDI_F;	//Minimum RFDI of a fluorine atom
	public double Tot_RFDI_F;	//Total RFDI of a fluorine atom
	public double Av_RFDI_F;		//Average RFDI of a fluorine atom
	public double Max_RFDI_Cl;	//maximum RFDI of a chlorine atom
	public double Min_RFDI_Cl;	//Minimum RFDI of a chlorine atom
	public double Tot_RFDI_Cl;	//Total RFDI of a chlorine atom
	public double Av_RFDI_Cl;	//Average RFDI of a chlorine atom
	public double Max_RFDI_Br;	//maximum RFDI of a bromine atom
	public double Min_RFDI_Br;	//Minimum RFDI of a bromine atom
	public double Tot_RFDI_Br;	//Total RFDI of a bromine atom
	public double Av_RFDI_Br;	//Average RFDI of a bromine atom
	public double Max_RFDI_I;	//maximum RFDI of an iodine atom
	public double Min_RFDI_I;	//Minimum RFDI of an iodine atom
	public double Tot_RFDI_I;	//Total RFDI of an iodine atom
	public double Av_RFDI_I;		//Average RFDI of an iodine atom
	public double Max_RFDI;		//Maximum RFDI of all atoms
	public double Min_RFDI;		//Minimum RFDI of all atoms
	public double Tot_RFDI;		//Total RFDI of all atoms
	public double Av_RFDI;		//Average RFDI of all atoms

	//Electrophilic Superdelocalizability (ESD) of atoms
	public double Max_ESD_C;	//maximum ESD of a carbon atom
	public double Min_ESD_C;	//Minimum ESD of a carbon atom
	public double Tot_ESD_C;	//Total ESD of a carbon atom
	public double Av_ESD_C;		//Average ESD of a carbon atom
	public double Max_ESD_N;	//maximum ESD of a nitrogen atom
	public double Min_ESD_N;	//Minimum ESD of a nitrogen atom
	public double Tot_ESD_N;	//Total ESD of a nitrogen atom
	public double Av_ESD_N;		//Average ESD of a nitrogen atom
	public double Max_ESD_P;	//maximum ESD of a phosphorous atom
	public double Min_ESD_P;	//Minimum ESD of a phosphorous atom
	public double Tot_ESD_P;	//Total ESD of a phosphorous atom
	public double Av_ESD_P;		//Average ESD of a phsophorous atom
	public double Max_ESD_O;	//maximum ESD of an oxygen atom
	public double Min_ESD_O;	//Minimum ESD of an oxygen atom
	public double Tot_ESD_O;	//Total ESD of an oxygen atom
	public double Av_ESD_O;		//Average ESD of an oxygen atom
	public double Max_ESD_S;	//maximum ESD of a sulfur atom
	public double Min_ESD_S;	//Minimum ESD of a sulfur atom
	public double Tot_ESD_S;	//Total ESD of a sulfur atom
	public double Av_ESD_S;		//Average ESD of a sulfur atom
	public double Max_ESD_Si;	//maximum ESD of a silicon atom
	public double Min_ESD_Si;	//Minimum ESD of a silicon atom
	public double Tot_ESD_Si;	//Total ESD of a silicon atom
	public double Av_ESD_Si;	//Average ESD of a silicon atom
	public double Max_ESD_Hg;	//maximum ESD of a mercury atom
	public double Min_ESD_Hg;	//Minimum ESD of a mercury atom
	public double Tot_ESD_Hg;	//Total ESD of a mercury atom
	public double Av_ESD_Hg;	//Average ESD of a mercury atom
	public double Max_ESD_H;	//maximum ESD of a hydrogen atom
	public double Min_ESD_H;	//Minimum ESD of a hydrogen atom
	public double Tot_ESD_H;	//Total ESD of a hydrogen atom
	public double Av_ESD_H;		//Average ESD of a hydrogen atom
	public double Max_ESD_F;	//maximum ESD of a fluorine atom
	public double Min_ESD_F;	//Minimum ESD of a fluorine atom
	public double Tot_ESD_F;	//Total ESD of a fluorine atom
	public double Av_ESD_F;		//Average ESD of a fluorine atom
	public double Max_ESD_Cl;	//maximum ESD of a chlorine atom
	public double Min_ESD_Cl;	//Minimum ESD of a chlorine atom
	public double Tot_ESD_Cl;	//Total ESD of a chlorine atom
	public double Av_ESD_Cl;	//Average ESD of a chlorine atom
	public double Max_ESD_Br;	//maximum ESD of a bromine atom
	public double Min_ESD_Br;	//Minimum ESD of a bromine atom
	public double Tot_ESD_Br;	//Total ESD of a bromine atom
	public double Av_ESD_Br;	//Average ESD of a bromine atom
	public double Max_ESD_I;	//maximum ESD of an iodine atom
	public double Min_ESD_I;	//Minimum ESD of an iodine atom
	public double Tot_ESD_I;	//Total ESD of an iodine atom
	public double Av_ESD_I;		//Average ESD of an iodine atom
	public double Max_ESD;		//Maximum ESD of all atoms
	public double Min_ESD;		//Minimum ESD of all atoms
	public double Tot_ESD;		//Total ESD of all atoms
	public double Av_ESD;		//Average ESD of all atoms

	//Nucleophilic Superdelocalizability (NSD) of atoms
	public double Max_NSD_C;	//maximum NSD of a carbon atom
	public double Min_NSD_C;	//Minimum NSD of a carbon atom
	public double Tot_NSD_C;	//Total NSD of a carbon atom
	public double Av_NSD_C;		//Average NSD of a carbon atom
	public double Max_NSD_N;	//maximum NSD of a nitrogen atom
	public double Min_NSD_N;	//Minimum NSD of a nitrogen atom
	public double Tot_NSD_N;	//Total NSD of a nitrogen atom
	public double Av_NSD_N;		//Average NSD of a nitrogen atom
	public double Max_NSD_P;	//maximum NSD of a phosphorous atom
	public double Min_NSD_P;	//Minimum NSD of a phosphorous atom
	public double Tot_NSD_P;	//Total NSD of a phosphorous atom
	public double Av_NSD_P;		//Average NSD of a phsophorous atom
	public double Max_NSD_O;	//maximum NSD of an oxygen atom
	public double Min_NSD_O;	//Minimum NSD of an oxygen atom
	public double Tot_NSD_O;	//Total NSD of an oxygen atom
	public double Av_NSD_O;		//Average NSD of an oxygen atom
	public double Max_NSD_S;	//maximum NSD of a sulfur atom
	public double Min_NSD_S;	//Minimum NSD of a sulfur atom
	public double Tot_NSD_S;	//Total NSD of a sulfur atom
	public double Av_NSD_S;		//Average NSD of a sulfur atom
	public double Max_NSD_Si;	//maximum NSD of a silicon atom
	public double Min_NSD_Si;	//Minimum NSD of a silicon atom
	public double Tot_NSD_Si;	//Total NSD of a silicon atom
	public double Av_NSD_Si;	//Average NSD of a silicon atom
	public double Max_NSD_Hg;	//maximum NSD of a mercury atom
	public double Min_NSD_Hg;	//Minimum NSD of a mercury atom
	public double Tot_NSD_Hg;	//Total NSD of a mercury atom
	public double Av_NSD_Hg;	//Average NSD of a mercury atom
	public double Max_NSD_H;	//maximum NSD of a hydrogen atom
	public double Min_NSD_H;	//Minimum NSD of a hydrogen atom
	public double Tot_NSD_H;	//Total NSD of a hydrogen atom
	public double Av_NSD_H;		//Average NSD of a hydrogen atom
	public double Max_NSD_F;	//maximum NSD of a fluorine atom
	public double Min_NSD_F;	//Minimum NSD of a fluorine atom
	public double Tot_NSD_F;	//Total NSD of a fluorine atom
	public double Av_NSD_F;		//Average NSD of a fluorine atom
	public double Max_NSD_Cl;	//maximum NSD of a chlorine atom
	public double Min_NSD_Cl;	//Minimum NSD of a chlorine atom
	public double Tot_NSD_Cl;	//Total NSD of a chlorine atom
	public double Av_NSD_Cl;	//Average NSD of a chlorine atom
	public double Max_NSD_Br;	//maximum NSD of a bromine atom
	public double Min_NSD_Br;	//Minimum NSD of a bromine atom
	public double Tot_NSD_Br;	//Total NSD of a bromine atom
	public double Av_NSD_Br;	//Average NSD of a bromine atom
	public double Max_NSD_I;	//maximum NSD of an iodine atom
	public double Min_NSD_I;	//Minimum NSD of an iodine atom
	public double Tot_NSD_I;	//Total NSD of an iodine atom
	public double Av_NSD_I;		//Average NSD of an iodine atom
	public double Max_NSD;		//Maximum NSD of all atoms
	public double Min_NSD;		//Minimum NSD of all atoms
	public double Tot_NSD;		//Total NSD of all atoms
	public double Av_NSD;		//Average NSD of all atoms

	//Valencies and bond orders of atoms
	public double Min_val_C;	//Minimum valency of a carbon atom
	public double Max_val_C;	//Maximum valency of a carbon atom
	public double Av_val_C;		//Average valency of a carbon atom
	public double Min_BO_C;		//Minimum bond order of a carbon atom
	public double Max_BO_C;		//Maximum bond order of a carbon atom
	public double Av_BO_C;		//Average bond order of a carbon atom
	public double Min_val_N;	//Minimum valency of a nitrogen atom
	public double Max_val_N;	//Maximum valency of a nitrogen atom
	public double Av_val_N;		//Average valency of a nitrogen atom
	public double Min_BO_N;		//Minimum bond order of a nitrogen atom
	public double Max_BO_N;		//Maximum bond order of a nitrogen atom
	public double Av_BO_N;		//Average bond order of a nitrogen atom
	public double Min_val_P;	//Minimum valency of a phosphorous atom
	public double Max_val_P;	//Maximum valency of a phosphorous atom
	public double Av_val_P;		//Average valency of a phsophorous atom
	public double Min_BO_P;		//Minimum bond order of a phosphorous atom
	public double Max_BO_P;		//Maximum bond order of a phsophorous atom
	public double Av_BO_P;		//Average bond order of a phosphorous atom
	public double Min_val_O;	//Minimum valency of an oxygen atom
	public double Max_val_O;	//Maximum valency of an oxygen atom
	public double Av_val_O;		//Average valency of an oxygen atom
	public double Min_BO_O;		//Minimum bond order of an oxygen atom
	public double Max_BO_O;		//Maximum bond order of an oxygen atom
	public double Av_BO_O;		//Average bond order of an oxygen atom
	public double Min_val_S;	//Minimum valency of a sulfur atom
	public double Max_val_S;	//Maximum valency of a sulfur atom
	public double Av_val_S;		//Average valency of a sulfur atom
	public double Min_BO_S;		//Minimum bond order of a sulfur atom
	public double Max_BO_S;		//Maximum bond order of a sulfur atom
	public double Av_BO_S;		//Average bond order of a sulfur atom
	public double Min_val_Si;	//Minimum valency of a silicon atom
	public double Max_val_Si;	//Maximum valency of a silicon atom
	public double Av_val_Si;	//Average valency of a silicon atom
	public double Min_BO_Si;	//Minimum bond order of a silicon atom
	public double Max_BO_Si;	//Maximum bond order of a silicon atom
	public double Av_BO_Si;		//Average bond order of a silicon atom
	public double Min_val_Hg;	//Minimum valency of a mercury atom
	public double Max_val_Hg;	//Maximum valency of a mercury atom
	public double Av_val_Hg;	//Average valency of a mercury atom
	public double Min_BO_Hg;	//Minimum bond order of a mercury atom
	public double Max_BO_Hg;	//Maximum bond order of a mercury atom
	public double Av_BO_Hg;		//Average bond order of a mercury atom
	public double Min_val_H;	//Minimum valency of a hydrogen atom
	public double Max_val_H;	//Maximum valency of a hydrogen atom
	public double Av_val_H;		//Average valency of a hydrogen atom
	public double Min_BO_H;		//Minimum bond order of a hydrogen atom
	public double Max_BO_H;		//Maximum bond order of a hydrogen atom
	public double Av_BO_H;		//Average bond order of a hydrogen atom
	public double Min_val_F;	//Minimum valency of a fluorine atom
	public double Max_val_F;	//Maximum valency of a fluorine atom
	public double Av_val_F;		//Average valency of a fluorine atom
	public double Min_BO_F;		//Minimum bond order of a fluorine atom
	public double Max_BO_F;		//Maximum bond order of a fluorine atom
	public double Av_BO_F;		//Average bond order of a fluorine atom
	public double Min_val_Cl;	//Minimum valency of a chlorine atom
	public double Max_val_Cl;	//Maximum valency of a chlorine atom
	public double Av_val_Cl;	//Average valency of a chlorine atom
	public double Min_BO_Cl;	//Minimum bond order of a chlorine atom
	public double Max_BO_Cl;	//Maximum bond order of a chlorine atom
	public double Av_BO_Cl;		//Average bond order of a chlorine atom
	public double Min_val_Br;	//Minimum valency of a bromine atom
	public double Max_val_Br;	//Maximum valency of a bromine atom
	public double Av_val_Br;	//Average valency of a bromine atom
	public double Min_BO_Br;	//Minimum bond order of a bromine atom
	public double Max_BO_Br;	//Maximum bond order of a bromine atom
	public double Av_BO_Br;		//Average bond order of a bromine atom
	public double Min_val_I;	//Minimum valency of an iodine atom
	public double Max_val_I;	//Maximum valency of an iodine atom
	public double Av_val_I;		//Average valency of an iodine atom
	public double Min_BO_I;		//Minimum bond order of an iodine atom
	public double Max_BO_I;		//Maximum bond order of an iodine atom
	public double Av_BO_I;		//Average bond order of an iodine atom

	//electron-neutron attraction, electron-electron repulsion and atomic state energy for atoms
	public double Min_en_C;		//Minimum e-n attraction energy for a carbon atom
	public double Max_en_C;		//Maximum e-n attraction energy for a carbon atom
	public double Min_ee_C;		//Minimum e-e repulsion energy for a carbon atom
	public double Max_ee_C;		//Maximum e-e repulsion energy for a carbon atom
	public double Min_ase_C;	//Minimum atomic state energy for a carbon atom
	public double Max_ase_C;	//Maximum atomic state energy for a carbon atom
	public double Min_en_N;		//Minimum e-n attraction energy for a nitrogen atom
	public double Max_en_N;		//Maximum e-n attraction energy for a nitrogen atom
	public double Min_ee_N;		//Minimum e-e repulsion energy for a nitrogen atom
	public double Max_ee_N;		//Maximum e-e repulsion energy for a nitrogen atom
	public double Min_ase_N;	//Minimum atomic state energy for a nitrogen atom
	public double Max_ase_N;	//Maximum atomic state energy for a nitrogen atom
	public double Min_en_P;		//Minimum e-n attraction energy for a phosphorous atom
	public double Max_en_P;		//Maximum e-n attraction energy for a phosphorous atom
	public double Min_ee_P;		//Minimum e-e repulsion energy for a phosphorous atom
	public double Max_ee_P;		//Maximum e-e repulsion energy for a phosphorous atom
	public double Min_ase_P;	//Minimum atomic state energy for a phosphorous atom
	public double Max_ase_P;	//Maximum atomic state energy for a phosphorous atom
	public double Min_en_O;		//Minimum e-n attraction energy for an oxygen atom
	public double Max_en_O;		//Maximum e-n attraction energy for an oxygen atom
	public double Min_ee_O;		//Minimum e-e repulsion energy for an oxygen atom
	public double Max_ee_O;		//Maximum e-e repulsion energy for an oxygen atom
	public double Min_ase_O;	//Minimum atomic state energy for an oxygen atom
	public double Max_ase_O;	//Maximum atomic state energy for an oxygen atom
	public double Min_en_S;		//Minimum e-n attraction energy for a sulfur atom
	public double Max_en_S;		//Maximum e-n attraction energy for a sulfur atom
	public double Min_ee_S;		//Minimum e-e repulsion energy for a sulfur atom
	public double Max_ee_S;		//Maximum e-e repulsion energy for a sulfur atom
	public double Min_ase_S;	//Minimum atomic state energy for a sulfur atom
	public double Max_ase_S;	//Maximum atomic state energy for a sulfur atom
	public double Min_en_Si;	//Minimum e-n attraction energy for a silicon atom
	public double Max_en_Si;	//Maximum e-n attraction energy for a silicon atom
	public double Min_ee_Si;	//Minimum e-e repulsion energy for a silicon atom
	public double Max_ee_Si;	//Maximum e-e repulsion energy for a silicon atom
	public double Min_ase_Si;	//Minimum atomic state energy for a silicon atom
	public double Max_ase_Si;	//Maximum atomic state energy for a silicon atom
	public double Min_en_Hg;	//Minimum e-n attraction energy for a mercury atom
	public double Max_en_Hg;	//Maximum e-n attraction energy for a mercury atom
	public double Min_ee_Hg;	//Minimum e-e repulsion energy for a mercury atom
	public double Max_ee_Hg;	//Maximum e-e repulsion energy for a mercury atom
	public double Min_ase_Hg;	//Minimum atomic state energy for a mercury atom
	public double Max_ase_Hg;	//Maximum atomic state energy for a mercury atom
	public double Min_en_H;		//Minimum e-n attraction energy for a hydrogen atom
	public double Max_en_H;		//Maximum e-n attraction energy for a hydrogen atom
	public double Min_ee_H;		//Minimum e-e repulsion energy for a hydrogen atom
	public double Max_ee_H;		//Maximum e-e repulsion energy for a hydrogen atom
	public double Min_ase_H;	//Minimum atomic state energy for a hydrogen atom
	public double Max_ase_H;	//Maximum atomic state energy for a hydrogen atom
	public double Min_en_F;		//Minimum e-n attraction energy for a fluorine atom
	public double Max_en_F;		//Maximum e-n attraction energy for a fluorine atom
	public double Min_ee_F;		//Minimum e-e repulsion energy for a fluorine atom
	public double Max_ee_F;		//Maximum e-e repulsion energy for a fluorine atom
	public double Min_ase_F;	//Minimum atomic state energy for a fluorine atom
	public double Max_ase_F;	//Maximum atomic state energy for a fluorine atom
	public double Min_en_Cl;	//Minimum e-n attraction energy for a chlorine atom
	public double Max_en_Cl;	//Maximum e-n attraction energy for a chlorine atom
	public double Min_ee_Cl;	//Minimum e-e repulsion energy for a chlorine atom
	public double Max_ee_Cl;	//Maximum e-e repulsion energy for a chlorine atom
	public double Min_ase_Cl;	//Minimum atomic state energy for a chlorine atom
	public double Max_ase_Cl;	//Maximum atomic state energy for a chlorine atom
	public double Min_en_Br;	//Minimum e-n attraction energy for a bromine atom
	public double Max_en_Br;	//Maximum e-n attraction energy for a bromine atom
	public double Min_ee_Br;	//Minimum e-e repulsion energy for a bromine atom
	public double Max_ee_Br;	//Maximum e-e repulsion energy for a bromine atom
	public double Min_ase_Br;	//Minimum atomic state energy for a bromine atom
	public double Max_ase_Br;	//Maximum atomic state energy for a bromine atom
	public double Min_en_I;		//Minimum e-n attraction energy for an iodine atom
	public double Max_en_I;		//Maximum e-n attraction energy for an iodine atom
	public double Min_ee_I;		//Minimum e-e repulsion energy for an iodine atom
	public double Max_ee_I;		//Maximum e-e repulsion energy for an iodine atom
	public double Min_ase_I;	//Minimum atomic state energy for an iodine atom
	public double Max_ase_I;	//Maximum atomic state energy for an iodine atom

	//Resonance, e-n, n-n, exchange and totoal interaction energy for bonds
	public double Min_re_C_C;	//Minimum resonance energy for a carbon-carbon bond
	public double Max_re_C_C;	//Maximum resonance energy for a carbon-carbon bond
	public double Min_en_C_C;	//Minimum e-n attraction for a carbon-carbon bond
	public double Max_en_C_C;	//Maximum e-n attraction for a carbon-carbon bond
	public double Min_nn_C_C;	//Minimum n-n repulsion for a carbon-carbon bond
	public double Max_nn_C_C;	//Maximum n-n repulsion for a carbon-carbon bond
	public double Min_eee_C_C;	//Minimum exch. + e-e repulsion energy for a carbon-carbon bond
	public double Max_eee_C_C;	//Maximum exch. + e-e repulsion energy for a carbon-carbon bond
	public double Min_tie_C_C;	//Minimum tot. interaction energy for a carbon-carbon bond
	public double Max_tie_C_C;	//Maximum tot. interaction energy for a carbon-carbon bond
	public double Min_re_C_N;	//Minimum resonance energy for a carbon-nitrogen bond
	public double Max_re_C_N;	//Maximum resonance energy for a carbon-nitrogen bond
	public double Min_en_C_N;	//Minimum e-n attraction for a carbon-nitrogen bond
	public double Max_en_C_N;	//Maximum e-n attraction for a carbon-nitrogen bond
	public double Min_nn_C_N;	//Minimum n-n repulsion for a carbon-nitrogen bond
	public double Max_nn_C_N;	//Maximum n-n repulsion for a carbon-nitrogen bond
	public double Min_eee_C_N;	//Minimum exch. + e-e repulsion energy for a carbon-nitrogen bond
	public double Max_eee_C_N;	//Maximum exch. + e-e repulsion energy for a carbon-nitrogen bond
	public double Min_tie_C_N;	//Minimum tot. interaction energy for a carbon-nitrogen bond
	public double Max_tie_C_N;	//Maximum tot. interaction energy for a carbon-nitrogen bond
	public double Min_re_C_P;	//Minimum resonance energy for a carbon-phosphorous bond
	public double Max_re_C_P;	//Maximum resonance energy for a carbon-phosphorous bond
	public double Min_en_C_P;	//Minimum e-n attraction for a carbon-phosphorous bond
	public double Max_en_C_P;	//Maximum e-n attraction for a carbon-phosphorous bond
	public double Min_nn_C_P;	//Minimum n-n repulsion for a carbon-phosphorous bond
	public double Max_nn_C_P;	//Maximum n-n repulsion for a carbon-phosphorous bond
	public double Min_eee_C_P;	//Minimum exch. + e-e repulsion energy for a carbon-phosphorous bond
	public double Max_eee_C_P;	//Maximum exch. + e-e repulsion energy for a carbon-phosphorous bond
	public double Min_tie_C_P;	//Minimum tot. interaction energy for a carbon-phosphorous bond
	public double Max_tie_C_P;	//Maximum tot. interaction energy for a carbon-phosphorous bond
	public double Min_re_C_O;	//Minimum resonance energy for a carbon-oxygen bond
	public double Max_re_C_O;	//Maximum resonance energy for a carbon-oxygen bond
	public double Min_en_C_O;	//Minimum e-n attraction for a carbon-oxygen bond
	public double Max_en_C_O;	//Maximum e-n attraction for a carbon-oxygen bond
	public double Min_nn_C_O;	//Minimum n-n repulsion for a carbon-oxygen bond
	public double Max_nn_C_O;	//Maximum n-n repulsion for a carbon-oxygen bond
	public double Min_eee_C_O;	//Minimum exch. + e-e repulsion energy for a carbon-oxygen bond
	public double Max_eee_C_O;	//Maximum exch. + e-e repulsion energy for a carbon-oxygen bond
	public double Min_tie_C_O;	//Minimum tot. interaction energy for a carbon-oxygen bond
	public double Max_tie_C_O;	//Maximum tot. interaction energy for a carbon-oxygen bond
	public double Min_re_C_S;	//Minimum resonance energy for a carbon-sulfur bond
	public double Max_re_C_S;	//Maximum resonance energy for a carbon-sulfur bond
	public double Min_en_C_S;	//Minimum e-n attraction for a carbon-sulfur bond
	public double Max_en_C_S;	//Maximum e-n attraction for a carbon-sulfur bond
	public double Min_nn_C_S;	//Minimum n-n repulsion for a carbon-sulfur bond
	public double Max_nn_C_S;	//Maximum n-n repulsion for a carbon-sulfur bond
	public double Min_eee_C_S;	//Minimum exch. + e-e repulsion energy for a carbon-sulfur bond
	public double Max_eee_C_S;	//Maximum exch. + e-e repulsion energy for a carbon-sulfur bond
	public double Min_tie_C_S;	//Minimum tot. interaction energy for a carbon-sulfur bond
	public double Max_tie_C_S;	//Maximum tot. interaction energy for a carbon-sulfur bond
	public double Min_re_C_Si;	//Minimum resonance energy for a carbon-silicon bond
	public double Max_re_C_Si;	//Maximum resonance energy for a carbon-silicon bond
	public double Min_en_C_Si;	//Minimum e-n attraction for a carbon-silicon bond
	public double Max_en_C_Si;	//Maximum e-n attraction for a carbon-silicon bond
	public double Min_nn_C_Si;	//Minimum n-n repulsion for a carbon-silicon bond
	public double Max_nn_C_Si;	//Maximum n-n repulsion for a carbon-silicon bond
	public double Min_eee_C_Si;	//Minimum exch. + e-e repulsion energy for a carbon-silicon bond
	public double Max_eee_C_Si;	//Maximum exch. + e-e repulsion energy for a carbon-silicon bond
	public double Min_tie_C_Si;	//Minimum tot. interaction energy for a carbon-silicon bond
	public double Max_tie_C_Si;	//Maximum tot. interaction energy for a carbon-silicon bond
	public double Min_re_C_Hg;	//Minimum resonance energy for a carbon-mercury bond
	public double Max_re_C_Hg;	//Maximum resonance energy for a carbon-mercury bond
	public double Min_en_C_Hg;	//Minimum e-n attraction for a carbon-mercury bond
	public double Max_en_C_Hg;	//Maximum e-n attraction for a carbon-mercury bond
	public double Min_nn_C_Hg;	//Minimum n-n repulsion for a carbon-mercury bond
	public double Max_nn_C_Hg;	//Maximum n-n repulsion for a carbon-mercury bond
	public double Min_eee_C_Hg;	//Minimum exch. + e-e repulsion energy for a carbon-mercury bond
	public double Max_eee_C_Hg;	//Maximum exch. + e-e repulsion energy for a carbon-mercury bond
	public double Min_tie_C_Hg;	//Minimum tot. interaction energy for a carbon-mercury bond
	public double Max_tie_C_Hg;	//Maximum tot. interaction energy for a carbon-mercury bond
	public double Min_re_C_H;	//Minimum resonance energy for a carbon-hydrogen bond
	public double Max_re_C_H;	//Maximum resonance energy for a carbon-hydrogen bond
	public double Min_en_C_H;	//Minimum e-n attraction for a carbon-hydrogen bond
	public double Max_en_C_H;	//Maximum e-n attraction for a carbon-hydrogen bond
	public double Min_nn_C_H;	//Minimum n-n repulsion for a carbon-hydrogen bond
	public double Max_nn_C_H;	//Maximum n-n repulsion for a carbon-hydrogen bond
	public double Min_eee_C_H;	//Minimum exch. + e-e repulsion energy for a carbon-hydrogen bond
	public double Max_eee_C_H;	//Maximum exch. + e-e repulsion energy for a carbon-hydrogen bond
	public double Min_tie_C_H;	//Minimum tot. interaction energy for a carbon-hydrogen bond
	public double Max_tie_C_H;	//Maximum tot. interaction energy for a carbon-hydrogen bond
	public double Min_re_C_F;	//Minimum resonance energy for a carbon-fluorine bond
	public double Max_re_C_F;	//Maximum resonance energy for a carbon-fluorine bond
	public double Min_en_C_F;	//Minimum e-n attraction for a carbon-fluorine bond
	public double Max_en_C_F;	//Maximum e-n attraction for a carbon-fluorine bond
	public double Min_nn_C_F;	//Minimum n-n repulsion for a carbon-fluorine bond
	public double Max_nn_C_F;	//Maximum n-n repulsion for a carbon-fluorine bond
	public double Min_eee_C_F;	//Minimum exch. + e-e repulsion energy for a carbon-fluorine bond
	public double Max_eee_C_F;	//Maximum exch. + e-e repulsion energy for a carbon-fluorine bond
	public double Min_tie_C_F;	//Minimum tot. interaction energy for a carbon-fluorine bond
	public double Max_tie_C_F;	//Maximum tot. interaction energy for a carbon-fluorine bond
	public double Min_re_C_Cl;	//Minimum resonance energy for a carbon-chlorine bond
	public double Max_re_C_Cl;	//Maximum resonance energy for a carbon-chlorine bond
	public double Min_en_C_Cl;	//Minimum e-n attraction for a carbon-chlorine bond
	public double Max_en_C_Cl;	//Maximum e-n attraction for a carbon-chlorine bond
	public double Min_nn_C_Cl;	//Minimum n-n repulsion for a carbon-chlorine bond
	public double Max_nn_C_Cl;	//Maximum n-n repulsion for a carbon-chlorine bond
	public double Min_eee_C_Cl;	//Minimum exch. + e-e repulsion energy for a carbon-chlorine bond
	public double Max_eee_C_Cl;	//Maximum exch. + e-e repulsion energy for a carbon-chlorine bond
	public double Min_tie_C_Cl;	//Minimum tot. interaction energy for a carbon-chlorine bond
	public double Max_tie_C_Cl;	//Maximum tot. interaction energy for a carbon-chlorine bond
	public double Min_re_C_Br;	//Minimum resonance energy for a carbon-bromine bond
	public double Max_re_C_Br;	//Maximum resonance energy for a carbon-bromine bond
	public double Min_en_C_Br;	//Minimum e-n attraction for a carbon-bromine bond
	public double Max_en_C_Br;	//Maximum e-n attraction for a carbon-bromine bond
	public double Min_nn_C_Br;	//Minimum n-n repulsion for a carbon-bromine bond
	public double Max_nn_C_Br;	//Maximum n-n repulsion for a carbon-bromine bond
	public double Min_eee_C_Br;	//Minimum exch. + e-e repulsion energy for a carbon-bromine bond
	public double Max_eee_C_Br;	//Maximum exch. + e-e repulsion energy for a carbon-bromine bond
	public double Min_tie_C_Br;	//Minimum tot. interaction energy for a carbon-bromine bond
	public double Max_tie_C_Br;	//Maximum tot. interaction energy for a carbon-bromine bond
	public double Min_re_C_I;	//Minimum resonance energy for a carbon-iodine bond
	public double Max_re_C_I;	//Maximum resonance energy for a carbon-iodine bond
	public double Min_en_C_I;	//Minimum e-n attraction for a carbon-iodine bond
	public double Max_en_C_I;	//Maximum e-n attraction for a carbon-iodine bond
	public double Min_nn_C_I;	//Minimum n-n repulsion for a carbon-iodine bond
	public double Max_nn_C_I;	//Maximum n-n repulsion for a carbon-iodine bond
	public double Min_eee_C_I;	//Minimum exch. + e-e repulsion energy for a carbon-iodine bond
	public double Max_eee_C_I;	//Maximum exch. + e-e repulsion energy for a carbon-iodine bond
	public double Min_tie_C_I;	//Minimum tot. interaction energy for a carbon-iodine bond
	public double Max_tie_C_I;	//Maximum tot. interaction energy for a carbon-iodine bond
	public double Min_re_N_N;	//Minimum resonance energy for a nitrogen-nitrogen bond
	public double Max_re_N_N;	//Maximum resonance energy for a nitrogen-nitrogen bond
	public double Min_en_N_N;	//Minimum e-n attraction for a nitrogen-nitrogen bond
	public double Max_en_N_N;	//Maximum e-n attraction for a nitrogen-nitrogen bond
	public double Min_nn_N_N;	//Minimum n-n repulsion for a nitrogen-nitrogen bond
	public double Max_nn_N_N;	//Maximum n-n repulsion for a nitrogen-nitrogen bond
	public double Min_eee_N_N;	//Minimum exch. + e-e repulsion energy for a nitrogen-nitrogen bond
	public double Max_eee_N_N;	//Maximum exch. + e-e repulsion energy for a nitrogen-nitrogen bond
	public double Min_tie_N_N;	//Minimum tot. interaction energy for a nitrogen-nitrogen bond
	public double Max_tie_N_N;	//Maximum tot. interaction energy for a nitrogen-nitrogen bond
	public double Min_re_N_P;	//Minimum resonance energy for a nitrogen-phosphorous bond
	public double Max_re_N_P;	//Maximum resonance energy for a nitrogen-phosphorous bond
	public double Min_en_N_P;	//Minimum e-n attraction for a nitrogen-phosphorous bond
	public double Max_en_N_P;	//Maximum e-n attraction for a nitrogen-phosphorous bond
	public double Min_nn_N_P;	//Minimum n-n repulsion for a nitrogen-phosphorous bond
	public double Max_nn_N_P;	//Maximum n-n repulsion for a nitrogen-phosphorous bond
	public double Min_eee_N_P;	//Minimum exch. + e-e repulsion energy for a nitrogen-phosphorous bond
	public double Max_eee_N_P;	//Maximum exch. + e-e repulsion energy for a nitrogen-phosphorous bond
	public double Min_tie_N_P;	//Minimum tot. interaction energy for a nitrogen-phosphorous bond
	public double Max_tie_N_P;	//Maximum tot. interaction energy for a nitrogen-phosphorous bond
	public double Min_re_N_O;	//Minimum resonance energy for a nitrogen-oxygen bond
	public double Max_re_N_O;	//Maximum resonance energy for a nitrogen-oxygen bond
	public double Min_en_N_O;	//Minimum e-n attraction for a nitrogen-oxygen bond
	public double Max_en_N_O;	//Maximum e-n attraction for a nitrogen-oxygen bond
	public double Min_nn_N_O;	//Minimum n-n repulsion for a nitrogen-oxygen bond
	public double Max_nn_N_O;	//Maximum n-n repulsion for a nitrogen-oxygen bond
	public double Min_eee_N_O;	//Minimum exch. + e-e repulsion energy for a nitrogen-oxygen bond
	public double Max_eee_N_O;	//Maximum exch. + e-e repulsion energy for a nitrogen-oxygen bond
	public double Min_tie_N_O;	//Minimum tot. interaction energy for a nitrogen-oxygen bond
	public double Max_tie_N_O;	//Maximum tot. interaction energy for a nitrogen-oxygen bond
	public double Min_re_N_S;	//Minimum resonance energy for a nitrogen-sulfur bond
	public double Max_re_N_S;	//Maximum resonance energy for a nitrogen-sulfur bond
	public double Min_en_N_S;	//Minimum e-n attraction for a nitrogen-sulfur bond
	public double Max_en_N_S;	//Maximum e-n attraction for a nitrogen-sulfur bond
	public double Min_nn_N_S;	//Minimum n-n repulsion for a nitrogen-sulfur bond
	public double Max_nn_N_S;	//Maximum n-n repulsion for a nitrogen-sulfur bond
	public double Min_eee_N_S;	//Minimum exch. + e-e repulsion energy for a nitrogen-sulfur bond
	public double Max_eee_N_S;	//Maximum exch. + e-e repulsion energy for a nitrogen-sulfur bond
	public double Min_tie_N_S;	//Minimum tot. interaction energy for a nitrogen-sulfur bond
	public double Max_tie_N_S;	//Maximum tot. interaction energy for a nitrogen-sulfur bond
	public double Min_re_N_Si;	//Minimum resonance energy for a nitrogen-silicon bond
	public double Max_re_N_Si;	//Maximum resonance energy for a nitrogen-silicon bond
	public double Min_en_N_Si;	//Minimum e-n attraction for a nitrogen-silicon bond
	public double Max_en_N_Si;	//Maximum e-n attraction for a nitrogen-silicon bond
	public double Min_nn_N_Si;	//Minimum n-n repulsion for a nitrogen-silicon bond
	public double Max_nn_N_Si;	//Maximum n-n repulsion for a nitrogen-silicon bond
	public double Min_eee_N_Si;	//Minimum exch. + e-e repulsion energy for a nitrogen-silicon bond
	public double Max_eee_N_Si;	//Maximum exch. + e-e repulsion energy for a nitrogen-silicon bond
	public double Min_tie_N_Si;	//Minimum tot. interaction energy for a nitrogen-silicon bond
	public double Max_tie_N_Si;	//Maximum tot. interaction energy for a nitrogen-silicon bond
	public double Min_re_N_Hg;	//Minimum resonance energy for a nitrogen-mercury bond
	public double Max_re_N_Hg;	//Maximum resonance energy for a nitrogen-mercury bond
	public double Min_en_N_Hg;	//Minimum e-n attraction for a nitrogen-mercury bond
	public double Max_en_N_Hg;	//Maximum e-n attraction for a nitrogen-mercury bond
	public double Min_nn_N_Hg;	//Minimum n-n repulsion for a nitrogen-mercury bond
	public double Max_nn_N_Hg;	//Maximum n-n repulsion for a nitrogen-mercury bond
	public double Min_eee_N_Hg;	//Minimum exch. + e-e repulsion energy for a nitrogen-mercury bond
	public double Max_eee_N_Hg;	//Maximum exch. + e-e repulsion energy for a nitrogen-mercury bond
	public double Min_tie_N_Hg;	//Minimum tot. interaction energy for a nitrogen-mercury bond
	public double Max_tie_N_Hg;	//Maximum tot. interaction energy for a nitrogen-mercury bond
	public double Min_re_N_H;	//Minimum resonance energy for a nitrogen-hydrogen bond
	public double Max_re_N_H;	//Maximum resonance energy for a nitrogen-hydrogen bond
	public double Min_en_N_H;	//Minimum e-n attraction for a nitrogen-hydrogen bond
	public double Max_en_N_H;	//Maximum e-n attraction for a nitrogen-hydrogen bond
	public double Min_nn_N_H;	//Minimum n-n repulsion for a nitrogen-hydrogen bond
	public double Max_nn_N_H;	//Maximum n-n repulsion for a nitrogen-hydrogen bond
	public double Min_eee_N_H;	//Minimum exch. + e-e repulsion energy for a nitrogen-hydrogen bond
	public double Max_eee_N_H;	//Maximum exch. + e-e repulsion energy for a nitrogen-hydrogen bond
	public double Min_tie_N_H;	//Minimum tot. interaction energy for a nitrogen-hydrogen bond
	public double Max_tie_N_H;	//Maximum tot. interaction energy for a nitrogen-hydrogen bond
	public double Min_re_N_F;	//Minimum resonance energy for a nitrogen-fluorine bond
	public double Max_re_N_F;	//Maximum resonance energy for a nitrogen-fluorine bond
	public double Min_en_N_F;	//Minimum e-n attraction for a nitrogen-fluorine bond
	public double Max_en_N_F;	//Maximum e-n attraction for a nitrogen-fluorine bond
	public double Min_nn_N_F;	//Minimum n-n repulsion for a nitrogen-fluorine bond
	public double Max_nn_N_F;	//Maximum n-n repulsion for a nitrogen-fluorine bond
	public double Min_eee_N_F;	//Minimum exch. + e-e repulsion energy for a nitrogen-fluorine bond
	public double Max_eee_N_F;	//Maximum exch. + e-e repulsion energy for a nitrogen-fluorine bond
	public double Min_tie_N_F;	//Minimum tot. interaction energy for a nitrogen-fluorine bond
	public double Max_tie_N_F;	//Maximum tot. interaction energy for a nitrogen-fluorine bond
	public double Min_re_N_Cl;	//Minimum resonance energy for a nitrogen-chlorine bond
	public double Max_re_N_Cl;	//Maximum resonance energy for a nitrogen-chlorine bond
	public double Min_en_N_Cl;	//Minimum e-n attraction for a nitrogen-chlorine bond
	public double Max_en_N_Cl;	//Maximum e-n attraction for a nitrogen-chlorine bond
	public double Min_nn_N_Cl;	//Minimum n-n repulsion for a nitrogen-chlorine bond
	public double Max_nn_N_Cl;	//Maximum n-n repulsion for a nitrogen-chlorine bond
	public double Min_eee_N_Cl;	//Minimum exch. + e-e repulsion energy for a nitrogen-chlorine bond
	public double Max_eee_N_Cl;	//Maximum exch. + e-e repulsion energy for a nitrogen-chlorine bond
	public double Min_tie_N_Cl;	//Minimum tot. interaction energy for a nitrogen-chlorine bond
	public double Max_tie_N_Cl;	//Maximum tot. interaction energy for a nitrogen-chlorine bond
	public double Min_re_N_Br;	//Minimum resonance energy for a nitrogen-bromine bond
	public double Max_re_N_Br;	//Maximum resonance energy for a nitrogen-bromine bond
	public double Min_en_N_Br;	//Minimum e-n attraction for a nitrogen-bromine bond
	public double Max_en_N_Br;	//Maximum e-n attraction for a nitrogen-bromine bond
	public double Min_nn_N_Br;	//Minimum n-n repulsion for a nitrogen-bromine bond
	public double Max_nn_N_Br;	//Maximum n-n repulsion for a nitrogen-bromine bond
	public double Min_eee_N_Br;	//Minimum exch. + e-e repulsion energy for a nitrogen-bromine bond
	public double Max_eee_N_Br;	//Maximum exch. + e-e repulsion energy for a nitrogen-bromine bond
	public double Min_tie_N_Br;	//Minimum tot. interaction energy for a nitrogen-bromine bond
	public double Max_tie_N_Br;	//Maximum tot. interaction energy for a nitrogen-bromine bond
	public double Min_re_N_I;	//Minimum resonance energy for a nitrogen-iodine bond
	public double Max_re_N_I;	//Maximum resonance energy for a nitrogen-iodine bond
	public double Min_en_N_I;	//Minimum e-n attraction for a nitrogen-iodine bond
	public double Max_en_N_I;	//Maximum e-n attraction for a nitrogen-iodine bond
	public double Min_nn_N_I;	//Minimum n-n repulsion for a nitrogen-iodine bond
	public double Max_nn_N_I;	//Maximum n-n repulsion for a nitrogen-iodine bond
	public double Min_eee_N_I;	//Minimum exch. + e-e repulsion energy for a nitrogen-iodine bond
	public double Max_eee_N_I;	//Maximum exch. + e-e repulsion energy for a nitrogen-iodine bond
	public double Min_tie_N_I;	//Minimum tot. interaction energy for a nitrogen-iodine bond
	public double Max_tie_N_I;	//Maximum tot. interaction energy for a nitrogen-iodine bond
	public double Min_re_P_P;	//Minimum resonance energy for a phosphorous-phosphorous bond
	public double Max_re_P_P;	//Maximum resonance energy for a phosphorous-phosphorous bond
	public double Min_en_P_P;	//Minimum e-n attraction for a phosphorous-phosphorous bond
	public double Max_en_P_P;	//Maximum e-n attraction for a phosphorous-phosphorous bond
	public double Min_nn_P_P;	//Minimum n-n repulsion for a phosphorous-phosphorous bond
	public double Max_nn_P_P;	//Maximum n-n repulsion for a phosphorous-phosphorous bond
	public double Min_eee_P_P;	//Minimum exch. + e-e repulsion energy for a phosphorous-phosphorous bond
	public double Max_eee_P_P;	//Maximum exch. + e-e repulsion energy for a phosphorous-phosphorous bond
	public double Min_tie_P_P;	//Minimum tot. interaction energy for a phosphorous-phosphorous bond
	public double Max_tie_P_P;	//Maximum tot. interaction energy for a phosphorous-phosphorous bond
	public double Min_re_P_O;	//Minimum resonance energy for a phosphorous-oxygen bond
	public double Max_re_P_O;	//Maximum resonance energy for a phosphorous-oxygen bond
	public double Min_en_P_O;	//Minimum e-n attraction for a phosphorous-oxygen bond
	public double Max_en_P_O;	//Maximum e-n attraction for a phosphorous-oxygen bond
	public double Min_nn_P_O;	//Minimum n-n repulsion for a phosphorous-oxygen bond
	public double Max_nn_P_O;	//Maximum n-n repulsion for a phosphorous-oxygen bond
	public double Min_eee_P_O;	//Minimum exch. + e-e repulsion energy for a phosphorous-oxygen bond
	public double Max_eee_P_O;	//Maximum exch. + e-e repulsion energy for a phosphorous-oxygen bond
	public double Min_tie_P_O;	//Minimum tot. interaction energy for a phosphorous-oxygen bond
	public double Max_tie_P_O;	//Maximum tot. interaction energy for a phosphorous-oxygen bond
	public double Min_re_P_S;	//Minimum resonance energy for a phosphorous-sulfur bond
	public double Max_re_P_S;	//Maximum resonance energy for a phosphorous-sulfur bond
	public double Min_en_P_S;	//Minimum e-n attraction for a phosphorous-sulfur bond
	public double Max_en_P_S;	//Maximum e-n attraction for a phosphorous-sulfur bond
	public double Min_nn_P_S;	//Minimum n-n repulsion for a phosphorous-sulfur bond
	public double Max_nn_P_S;	//Maximum n-n repulsion for a phosphorous-sulfur bond
	public double Min_eee_P_S;	//Minimum exch. + e-e repulsion energy for a phosphorous-sulfur bond
	public double Max_eee_P_S;	//Maximum exch. + e-e repulsion energy for a phosphorous-sulfur bond
	public double Min_tie_P_S;	//Minimum tot. interaction energy for a phosphorous-sulfur bond
	public double Max_tie_P_S;	//Maximum tot. interaction energy for a phosphorous-sulfur bond
	public double Min_re_P_Si;	//Minimum resonance energy for a phosphorous-silicon bond
	public double Max_re_P_Si;	//Maximum resonance energy for a phosphorous-silicon bond
	public double Min_en_P_Si;	//Minimum e-n attraction for a phosphorous-silicon bond
	public double Max_en_P_Si;	//Maximum e-n attraction for a phosphorous-silicon bond
	public double Min_nn_P_Si;	//Minimum n-n repulsion for a phosphorous-silicon bond
	public double Max_nn_P_Si;	//Maximum n-n repulsion for a phosphorous-silicon bond
	public double Min_eee_P_Si;	//Minimum exch. + e-e repulsion energy for a phosphorous-silicon bond
	public double Max_eee_P_Si;	//Maximum exch. + e-e repulsion energy for a phosphorous-silicon bond
	public double Min_tie_P_Si;	//Minimum tot. interaction energy for a phosphorous-silicon bond
	public double Max_tie_P_Si;	//Maximum tot. interaction energy for a phosphorous-silicon bond
	public double Min_re_P_Hg;	//Minimum resonance energy for a phosphorous-mercury bond
	public double Max_re_P_Hg;	//Maximum resonance energy for a phosphorous-mercury bond
	public double Min_en_P_Hg;	//Minimum e-n attraction for a phosphorous-mercury bond
	public double Max_en_P_Hg;	//Maximum e-n attraction for a phosphorous-mercury bond
	public double Min_nn_P_Hg;	//Minimum n-n repulsion for a phosphorous-mercury bond
	public double Max_nn_P_Hg;	//Maximum n-n repulsion for a phosphorous-mercury bond
	public double Min_eee_P_Hg;	//Minimum exch. + e-e repulsion energy for a phosphorous-mercury bond
	public double Max_eee_P_Hg;	//Maximum exch. + e-e repulsion energy for a phosphorous-mercury bond
	public double Min_tie_P_Hg;	//Minimum tot. interaction energy for a phosphorous-mercury bond
	public double Max_tie_P_Hg;	//Maximum tot. interaction energy for a phosphorous-mercury bond
	public double Min_re_P_H;	//Minimum resonance energy for a phosphorous-hydrogen bond
	public double Max_re_P_H;	//Maximum resonance energy for a phosphorous-hydrogen bond
	public double Min_en_P_H;	//Minimum e-n attraction for a phosphorous-hydrogen bond
	public double Max_en_P_H;	//Maximum e-n attraction for a phosphorous-hydrogen bond
	public double Min_nn_P_H;	//Minimum n-n repulsion for a phosphorous-hydrogen bond
	public double Max_nn_P_H;	//Maximum n-n repulsion for a phosphorous-hydrogen bond
	public double Min_eee_P_H;	//Minimum exch. + e-e repulsion energy for a phosphorous-hydrogen bond
	public double Max_eee_P_H;	//Maximum exch. + e-e repulsion energy for a phosphorous-hydrogen bond
	public double Min_tie_P_H;	//Minimum tot. interaction energy for a phosphorous-hydrogen bond
	public double Max_tie_P_H;	//Maximum tot. interaction energy for a phosphorous-hydrogen bond
	public double Min_re_P_F;	//Minimum resonance energy for a phosphorous-fluorine bond
	public double Max_re_P_F;	//Maximum resonance energy for a phosphorous-fluorine bond
	public double Min_en_P_F;	//Minimum e-n attraction for a phosphorous-fluorine bond
	public double Max_en_P_F;	//Maximum e-n attraction for a phosphorous-fluorine bond
	public double Min_nn_P_F;	//Minimum n-n repulsion for a phosphorous-fluorine bond
	public double Max_nn_P_F;	//Maximum n-n repulsion for a phosphorous-fluorine bond
	public double Min_eee_P_F;	//Minimum exch. + e-e repulsion energy for a phosphorous-fluorine bond
	public double Max_eee_P_F;	//Maximum exch. + e-e repulsion energy for a phosphorous-fluorine bond
	public double Min_tie_P_F;	//Minimum tot. interaction energy for a phosphorous-fluorine bond
	public double Max_tie_P_F;	//Maximum tot. interaction energy for a phosphorous-fluorine bond
	public double Min_re_P_Cl;	//Minimum resonance energy for a phosphorous-chlorine bond
	public double Max_re_P_Cl;	//Maximum resonance energy for a phosphorous-chlorine bond
	public double Min_en_P_Cl;	//Minimum e-n attraction for a phosphorous-chlorine bond
	public double Max_en_P_Cl;	//Maximum e-n attraction for a phosphorous-chlorine bond
	public double Min_nn_P_Cl;	//Minimum n-n repulsion for a phosphorous-chlorine bond
	public double Max_nn_P_Cl;	//Maximum n-n repulsion for a phosphorous-chlorine bond
	public double Min_eee_P_Cl;	//Minimum exch. + e-e repulsion energy for a phosphorous-chlorine bond
	public double Max_eee_P_Cl;	//Maximum exch. + e-e repulsion energy for a phosphorous-chlorine bond
	public double Min_tie_P_Cl;	//Minimum tot. interaction energy for a phosphorous-chlorine bond
	public double Max_tie_P_Cl;	//Maximum tot. interaction energy for a phosphorous-chlorine bond
	public double Min_re_P_Br;	//Minimum resonance energy for a phosphorous-bromine bond
	public double Max_re_P_Br;	//Maximum resonance energy for a phosphorous-bromine bond
	public double Min_en_P_Br;	//Minimum e-n attraction for a phosphorous-bromine bond
	public double Max_en_P_Br;	//Maximum e-n attraction for a phosphorous-bromine bond
	public double Min_nn_P_Br;	//Minimum n-n repulsion for a phosphorous-bromine bond
	public double Max_nn_P_Br;	//Maximum n-n repulsion for a phosphorous-bromine bond
	public double Min_eee_P_Br;	//Minimum exch. + e-e repulsion energy for a phosphorous-bromine bond
	public double Max_eee_P_Br;	//Maximum exch. + e-e repulsion energy for a phosphorous-bromine bond
	public double Min_tie_P_Br;	//Minimum tot. interaction energy for a phosphorous-bromine bond
	public double Max_tie_P_Br;	//Maximum tot. interaction energy for a phosphorous-bromine bond
	public double Min_re_P_I;	//Minimum resonance energy for a phosphorous-iodine bond
	public double Max_re_P_I;	//Maximum resonance energy for a phosphorous-iodine bond
	public double Min_en_P_I;	//Minimum e-n attraction for a phosphorous-iodine bond
	public double Max_en_P_I;	//Maximum e-n attraction for a phosphorous-iodine bond
	public double Min_nn_P_I;	//Minimum n-n repulsion for a phosphorous-iodine bond
	public double Max_nn_P_I;	//Maximum n-n repulsion for a phosphorous-iodine bond
	public double Min_eee_P_I;	//Minimum exch. + e-e repulsion energy for a phosphorous-iodine bond
	public double Max_eee_P_I;	//Maximum exch. + e-e repulsion energy for a phosphorous-iodine bond
	public double Min_tie_P_I;	//Minimum tot. interaction energy for a phosphorous-iodine bond
	public double Max_tie_P_I;	//Maximum tot. interaction energy for a phosphorous-iodine bond
	public double Min_re_O_O;	//Minimum resonance energy for an oxygen-oxygen bond
	public double Max_re_O_O;	//Maximum resonance energy for an oxygen-oxygen bond
	public double Min_en_O_O;	//Minimum e-n attraction for an oxygen-oxygen bond
	public double Max_en_O_O;	//Maximum e-n attraction for an oxygen-oxygen bond
	public double Min_nn_O_O;	//Minimum n-n repulsion for an oxygen-oxygen bond
	public double Max_nn_O_O;	//Maximum n-n repulsion for an oxygen-oxygen bond
	public double Min_eee_O_O;	//Minimum exch. + e-e repulsion energy for an oxygen-oxygen bond
	public double Max_eee_O_O;	//Maximum exch. + e-e repulsion energy for an oxygen-oxygen bond
	public double Min_tie_O_O;	//Minimum tot. interaction energy for an oxygen-oxygen bond
	public double Max_tie_O_O;	//Maximum tot. interaction energy for an oxygen-oxygen bond
	public double Min_re_O_S;	//Minimum resonance energy for an oxygen-sulfur bond
	public double Max_re_O_S;	//Maximum resonance energy for an oxygen-sulfur bond
	public double Min_en_O_S;	//Minimum e-n attraction for an oxygen-sulfur bond
	public double Max_en_O_S;	//Maximum e-n attraction for an oxygen-sulfur bond
	public double Min_nn_O_S;	//Minimum n-n repulsion for an oxygen-sulfur bond
	public double Max_nn_O_S;	//Maximum n-n repulsion for an oxygen-sulfur bond
	public double Min_eee_O_S;	//Minimum exch. + e-e repulsion energy for an oxygen-sulfur bond
	public double Max_eee_O_S;	//Maximum exch. + e-e repulsion energy for an oxygen-sulfur bond
	public double Min_tie_O_S;	//Minimum tot. interaction energy for an oxygen-sulfur bond
	public double Max_tie_O_S;	//Maximum tot. interaction energy for an oxygen-sulfur bond
	public double Min_re_O_Si;	//Minimum resonance energy for an oxygen-silicon bond
	public double Max_re_O_Si;	//Maximum resonance energy for an oxygen-silicon bond
	public double Min_en_O_Si;	//Minimum e-n attraction for an oxygen-silicon bond
	public double Max_en_O_Si;	//Maximum e-n attraction for an oxygen-silicon bond
	public double Min_nn_O_Si;	//Minimum n-n repulsion for an oxygen-silicon bond
	public double Max_nn_O_Si;	//Maximum n-n repulsion for an oxygen-silicon bond
	public double Min_eee_O_Si;	//Minimum exch. + e-e repulsion energy for an oxygen-silicon bond
	public double Max_eee_O_Si;	//Maximum exch. + e-e repulsion energy for an oxygen-silicon bond
	public double Min_tie_O_Si;	//Minimum tot. interaction energy for an oxygen-silicon bond
	public double Max_tie_O_Si;	//Maximum tot. interaction energy for an oxygen-silicon bond
	public double Min_re_O_Hg;	//Minimum resonance energy for an oxygen-mercury bond
	public double Max_re_O_Hg;	//Maximum resonance energy for an oxygen-mercury bond
	public double Min_en_O_Hg;	//Minimum e-n attraction for an oxygen-mercury bond
	public double Max_en_O_Hg;	//Maximum e-n attraction for an oxygen-mercury bond
	public double Min_nn_O_Hg;	//Minimum n-n repulsion for an oxygen-mercury bond
	public double Max_nn_O_Hg;	//Maximum n-n repulsion for an oxygen-mercury bond
	public double Min_eee_O_Hg;	//Minimum exch. + e-e repulsion energy for an oxygen-mercury bond
	public double Max_eee_O_Hg;	//Maximum exch. + e-e repulsion energy for an oxygen-mercury bond
	public double Min_tie_O_Hg;	//Minimum tot. interaction energy for an oxygen-mercury bond
	public double Max_tie_O_Hg;	//Maximum tot. interaction energy for an oxygen-mercury bond
	public double Min_re_O_H;	//Minimum resonance energy for an oxygen-hydrogen bond
	public double Max_re_O_H;	//Maximum resonance energy for an oxygen-hydrogen bond
	public double Min_en_O_H;	//Minimum e-n attraction for an oxygen-hydrogen bond
	public double Max_en_O_H;	//Maximum e-n attraction for an oxygen-hydrogen bond
	public double Min_nn_O_H;	//Minimum n-n repulsion for an oxygen-hydrogen bond
	public double Max_nn_O_H;	//Maximum n-n repulsion for an oxygen-hydrogen bond
	public double Min_eee_O_H;	//Minimum exch. + e-e repulsion energy for an oxygen-hydrogen bond
	public double Max_eee_O_H;	//Maximum exch. + e-e repulsion energy for an oxygen-hydrogen bond
	public double Min_tie_O_H;	//Minimum tot. interaction energy for an oxygen-hydrogen bond
	public double Max_tie_O_H;	//Maximum tot. interaction energy for an oxygen-hydrogen bond
	public double Min_re_O_F;	//Minimum resonance energy for an oxygen-fluorine bond
	public double Max_re_O_F;	//Maximum resonance energy for an oxygen-fluorine bond
	public double Min_en_O_F;	//Minimum e-n attraction for an oxygen-fluorine bond
	public double Max_en_O_F;	//Maximum e-n attraction for an oxygen-fluorine bond
	public double Min_nn_O_F;	//Minimum n-n repulsion for an oxygen-fluorine bond
	public double Max_nn_O_F;	//Maximum n-n repulsion for an oxygen-fluorine bond
	public double Min_eee_O_F;	//Minimum exch. + e-e repulsion energy for an oxygen-fluorine bond
	public double Max_eee_O_F;	//Maximum exch. + e-e repulsion energy for an oxygen-fluorine bond
	public double Min_tie_O_F;	//Minimum tot. interaction energy for an oxygen-fluorine bond
	public double Max_tie_O_F;	//Maximum tot. interaction energy for an oxygen-fluorine bond
	public double Min_re_O_Cl;	//Minimum resonance energy for an oxygen-chlorine bond
	public double Max_re_O_Cl;	//Maximum resonance energy for an oxygen-chlorine bond
	public double Min_en_O_Cl;	//Minimum e-n attraction for an oxygen-chlorine bond
	public double Max_en_O_Cl;	//Maximum e-n attraction for an oxygen-chlorine bond
	public double Min_nn_O_Cl;	//Minimum n-n repulsion for an oxygen-chlorine bond
	public double Max_nn_O_Cl;	//Maximum n-n repulsion for an oxygen-chlorine bond
	public double Min_eee_O_Cl;	//Minimum exch. + e-e repulsion energy for an oxygen-chlorine bond
	public double Max_eee_O_Cl;	//Maximum exch. + e-e repulsion energy for an oxygen-chlorine bond
	public double Min_tie_O_Cl;	//Minimum tot. interaction energy for an oxygen-chlorine bond
	public double Max_tie_O_Cl;	//Maximum tot. interaction energy for an oxygen-chlorine bond
	public double Min_re_O_Br;	//Minimum resonance energy for an oxygen-bromine bond
	public double Max_re_O_Br;	//Maximum resonance energy for an oxygen-bromine bond
	public double Min_en_O_Br;	//Minimum e-n attraction for an oxygen-bromine bond
	public double Max_en_O_Br;	//Maximum e-n attraction for an oxygen-bromine bond
	public double Min_nn_O_Br;	//Minimum n-n repulsion for an oxygen-bromine bond
	public double Max_nn_O_Br;	//Maximum n-n repulsion for an oxygen-bromine bond
	public double Min_eee_O_Br;	//Minimum exch. + e-e repulsion energy for an oxygen-bromine bond
	public double Max_eee_O_Br;	//Maximum exch. + e-e repulsion energy for an oxygen-bromine bond
	public double Min_tie_O_Br;	//Minimum tot. interaction energy for an oxygen-bromine bond
	public double Max_tie_O_Br;	//Maximum tot. interaction energy for an oxygen-bromine bond
	public double Min_re_O_I;	//Minimum resonance energy for an oxygen-iodine bond
	public double Max_re_O_I;	//Maximum resonance energy for an oxygen-iodine bond
	public double Min_en_O_I;	//Minimum e-n attraction for an oxygen-iodine bond
	public double Max_en_O_I;	//Maximum e-n attraction for an oxygen-iodine bond
	public double Min_nn_O_I;	//Minimum n-n repulsion for an oxygen-iodine bond
	public double Max_nn_O_I;	//Maximum n-n repulsion for an oxygen-iodine bond
	public double Min_eee_O_I;	//Minimum exch. + e-e repulsion energy for an oxygen-iodine bond
	public double Max_eee_O_I;	//Maximum exch. + e-e repulsion energy for an oxygen-iodine bond
	public double Min_tie_O_I;	//Minimum tot. interaction energy for an oxygen-iodine bond
	public double Max_tie_O_I;	//Maximum tot. interaction energy for an oxygen-iodine bond
	public double Min_re_S_S;	//Minimum resonance energy for a sulfur-sulfur bond
	public double Max_re_S_S;	//Maximum resonance energy for a sulfur-sulfur bond
	public double Min_en_S_S;	//Minimum e-n attraction for a sulfur-sulfur bond
	public double Max_en_S_S;	//Maximum e-n attraction for a sulfur-sulfur bond
	public double Min_nn_S_S;	//Minimum n-n repulsion for a sulfur-sulfur bond
	public double Max_nn_S_S;	//Maximum n-n repulsion for a sulfur-sulfur bond
	public double Min_eee_S_S;	//Minimum exch. + e-e repulsion energy for a sulfur-sulfur bond
	public double Max_eee_S_S;	//Maximum exch. + e-e repulsion energy for a sulfur-sulfur bond
	public double Min_tie_S_S;	//Minimum tot. interaction energy for a sulfur-sulfur bond
	public double Max_tie_S_S;	//Maximum tot. interaction energy for a sulfur-sulfur bond
	public double Min_re_S_Si;	//Minimum resonance energy for a sulfur-silicon bond
	public double Max_re_S_Si;	//Maximum resonance energy for a sulfur-silicon bond
	public double Min_en_S_Si;	//Minimum e-n attraction for a sulfur-silicon bond
	public double Max_en_S_Si;	//Maximum e-n attraction for a sulfur-silicon bond
	public double Min_nn_S_Si;	//Minimum n-n repulsion for a sulfur-silicon bond
	public double Max_nn_S_Si;	//Maximum n-n repulsion for a sulfur-silicon bond
	public double Min_eee_S_Si;	//Minimum exch. + e-e repulsion energy for a sulfur-silicon bond
	public double Max_eee_S_Si;	//Maximum exch. + e-e repulsion energy for a sulfur-silicon bond
	public double Min_tie_S_Si;	//Minimum tot. interaction energy for a sulfur-silicon bond
	public double Max_tie_S_Si;	//Maximum tot. interaction energy for a sulfur-silicon bond
	public double Min_re_S_Hg;	//Minimum resonance energy for a sulfur-silicon bond
	public double Max_re_S_Hg;	//Maximum resonance energy for a sulfur-silicon bond
	public double Min_en_S_Hg;	//Minimum e-n attraction for a sulfur-silicon bond
	public double Max_en_S_Hg;	//Maximum e-n attraction for a sulfur-silicon bond
	public double Min_nn_S_Hg;	//Minimum n-n repulsion for a sulfur-silicon bond
	public double Max_nn_S_Hg;	//Maximum n-n repulsion for a sulfur-silicon bond
	public double Min_eee_S_Hg;	//Minimum exch. + e-e repulsion energy for a sulfur-silicon bond
	public double Max_eee_S_Hg;	//Maximum exch. + e-e repulsion energy for a sulfur-silicon bond
	public double Min_tie_S_Hg;	//Minimum tot. interaction energy for a sulfur-silicon bond
	public double Max_tie_S_Hg;	//Maximum tot. interaction energy for a sulfur-silicon bond
	public double Min_re_S_H;	//Minimum resonance energy for a sulfur-hydrogen bond
	public double Max_re_S_H;	//Maximum resonance energy for a sulfur-hydrogen bond
	public double Min_en_S_H;	//Minimum e-n attraction for a sulfur-hydrogen bond
	public double Max_en_S_H;	//Maximum e-n attraction for a sulfur-hydrogen bond
	public double Min_nn_S_H;	//Minimum n-n repulsion for a sulfur-hydrogen bond
	public double Max_nn_S_H;	//Maximum n-n repulsion for a sulfur-hydrogen bond
	public double Min_eee_S_H;	//Minimum exch. + e-e repulsion energy for a sulfur-hydrogen bond
	public double Max_eee_S_H;	//Maximum exch. + e-e repulsion energy for a sulfur-hydrogen bond
	public double Min_tie_S_H;	//Minimum tot. interaction energy for a sulfur-hydrogen bond
	public double Max_tie_S_H;	//Maximum tot. interaction energy for a sulfur-hydrogen bond
	public double Min_re_S_F;	//Minimum resonance energy for a sulfur-fluorine bond
	public double Max_re_S_F;	//Maximum resonance energy for a sulfur-fluorine bond
	public double Min_en_S_F;	//Minimum e-n attraction for a sulfur-fluorine bond
	public double Max_en_S_F;	//Maximum e-n attraction for a sulfur-fluorine bond
	public double Min_nn_S_F;	//Minimum n-n repulsion for a sulfur-fluorine bond
	public double Max_nn_S_F;	//Maximum n-n repulsion for a sulfur-fluorine bond
	public double Min_eee_S_F;	//Minimum exch. + e-e repulsion energy for a sulfur-fluorine bond
	public double Max_eee_S_F;	//Maximum exch. + e-e repulsion energy for a sulfur-fluorine bond
	public double Min_tie_S_F;	//Minimum tot. interaction energy for a sulfur-fluorine bond
	public double Max_tie_S_F;	//Maximum tot. interaction energy for a sulfur-fluorine bond
	public double Min_re_S_Cl;	//Minimum resonance energy for a sulfur-chlorine bond
	public double Max_re_S_Cl;	//Maximum resonance energy for a sulfur-chlorine bond
	public double Min_en_S_Cl;	//Minimum e-n attraction for a sulfur-chlorine bond
	public double Max_en_S_Cl;	//Maximum e-n attraction for a sulfur-chlorine bond
	public double Min_nn_S_Cl;	//Minimum n-n repulsion for a sulfur-chlorine bond
	public double Max_nn_S_Cl;	//Maximum n-n repulsion for a sulfur-chlorine bond
	public double Min_eee_S_Cl;	//Minimum exch. + e-e repulsion energy for a sulfur-chlorine bond
	public double Max_eee_S_Cl;	//Maximum exch. + e-e repulsion energy for a sulfur-chlorine bond
	public double Min_tie_S_Cl;	//Minimum tot. interaction energy for a sulfur-chlorine bond
	public double Max_tie_S_Cl;	//Maximum tot. interaction energy for a sulfur-chlorine bond
	public double Min_re_S_Br;	//Minimum resonance energy for a sulfur-bromine bond
	public double Max_re_S_Br;	//Maximum resonance energy for a sulfur-bromine bond
	public double Min_en_S_Br;	//Minimum e-n attraction for a sulfur-bromine bond
	public double Max_en_S_Br;	//Maximum e-n attraction for a sulfur-bromine bond
	public double Min_nn_S_Br;	//Minimum n-n repulsion for a sulfur-bromine bond
	public double Max_nn_S_Br;	//Maximum n-n repulsion for a sulfur-bromine bond
	public double Min_eee_S_Br;	//Minimum exch. + e-e repulsion energy for a sulfur-bromine bond
	public double Max_eee_S_Br;	//Maximum exch. + e-e repulsion energy for a sulfur-bromine bond
	public double Min_tie_S_Br;	//Minimum tot. interaction energy for a sulfur-bromine bond
	public double Max_tie_S_Br;	//Maximum tot. interaction energy for a sulfur-bromine bond
	public double Min_re_S_I;	//Minimum resonance energy for a sulfur-iodine bond
	public double Max_re_S_I;	//Maximum resonance energy for a sulfur-iodine bond
	public double Min_en_S_I;	//Minimum e-n attraction for a sulfur-iodine bond
	public double Max_en_S_I;	//Maximum e-n attraction for a sulfur-iodine bond
	public double Min_nn_S_I;	//Minimum n-n repulsion for a sulfur-iodine bond
	public double Max_nn_S_I;	//Maximum n-n repulsion for a sulfur-iodine bond
	public double Min_eee_S_I;	//Minimum exch. + e-e repulsion energy for a sulfur-iodine bond
	public double Max_eee_S_I;	//Maximum exch. + e-e repulsion energy for a sulfur-iodine bond
	public double Min_tie_S_I;	//Minimum tot. interaction energy for a sulfur-iodine bond
	public double Max_tie_S_I;	//Maximum tot. interaction energy for a sulfur-iodine bond
	public double Min_re_Si_Si;	//Minimum resonance energy for a silicon-silicon bond
	public double Max_re_Si_Si;	//Maximum resonance energy for a silicon-silicon bond
	public double Min_en_Si_Si;	//Minimum e-n attraction for a silicon-silicon bond
	public double Max_en_Si_Si;	//Maximum e-n attraction for a silicon-silicon bond
	public double Min_nn_Si_Si;	//Minimum n-n repulsion for a silicon-silicon bond
	public double Max_nn_Si_Si;	//Maximum n-n repulsion for a silicon-silicon bond
	public double Min_eee_Si_Si;	//Minimum exch. + e-e repulsion energy for a silicon-silicon bond
	public double Max_eee_Si_Si;	//Maximum exch. + e-e repulsion energy for a silicon-silicon bond
	public double Min_tie_Si_Si;	//Minimum tot. interaction energy for a silicon-silicon bond
	public double Max_tie_Si_Si;	//Maximum tot. interaction energy for a silicon-silicon bond
	public double Min_re_Si_Hg;	//Minimum resonance energy for a silicon-mercury bond
	public double Max_re_Si_Hg;	//Maximum resonance energy for a silicon-mercury bond
	public double Min_en_Si_Hg;	//Minimum e-n attraction for a silicon-mercury bond
	public double Max_en_Si_Hg;	//Maximum e-n attraction for a silicon-mercury bond
	public double Min_nn_Si_Hg;	//Minimum n-n repulsion for a silicon-mercury bond
	public double Max_nn_Si_Hg;	//Maximum n-n repulsion for a silicon-mercury bond
	public double Min_eee_Si_Hg;	//Minimum exch. + e-e repulsion energy for a silicon-mercury bond
	public double Max_eee_Si_Hg;	//Maximum exch. + e-e repulsion energy for a silicon-mercury bond
	public double Min_tie_Si_Hg;	//Minimum tot. interaction energy for a silicon-mercury bond
	public double Max_tie_Si_Hg;	//Maximum tot. interaction energy for a silicon-mercury bond
	public double Min_re_Si_H;	//Minimum resonance energy for a silicon-hydrogen bond
	public double Max_re_Si_H;	//Maximum resonance energy for a silicon-hydrogen bond
	public double Min_en_Si_H;	//Minimum e-n attraction for a silicon-hydrogen bond
	public double Max_en_Si_H;	//Maximum e-n attraction for a silicon-hydrogen bond
	public double Min_nn_Si_H;	//Minimum n-n repulsion for a silicon-hydrogen bond
	public double Max_nn_Si_H;	//Maximum n-n repulsion for a silicon-hydrogen bond
	public double Min_eee_Si_H;	//Minimum exch. + e-e repulsion energy for a silicon-hydrogen bond
	public double Max_eee_Si_H;	//Maximum exch. + e-e repulsion energy for a silicon-hydrogen bond
	public double Min_tie_Si_H;	//Minimum tot. interaction energy for a silicon-hydrogen bond
	public double Max_tie_Si_H;	//Maximum tot. interaction energy for a silicon-hydrogen bond
	public double Min_re_Si_F;	//Minimum resonance energy for a silicon-fluorine bond
	public double Max_re_Si_F;	//Maximum resonance energy for a silicon-fluorine bond
	public double Min_en_Si_F;	//Minimum e-n attraction for a silicon-fluorine bond
	public double Max_en_Si_F;	//Maximum e-n attraction for a silicon-fluorine bond
	public double Min_nn_Si_F;	//Minimum n-n repulsion for a silicon-fluorine bond
	public double Max_nn_Si_F;	//Maximum n-n repulsion for a silicon-fluorine bond
	public double Min_eee_Si_F;	//Minimum exch. + e-e repulsion energy for a silicon-fluorine bond
	public double Max_eee_Si_F;	//Maximum exch. + e-e repulsion energy for a silicon-fluorine bond
	public double Min_tie_Si_F;	//Minimum tot. interaction energy for a silicon-fluorine bond
	public double Max_tie_Si_F;	//Maximum tot. interaction energy for a silicon-fluorine bond
	public double Min_re_Si_Cl;	//Minimum resonance energy for a silicon-chlorine bond
	public double Max_re_Si_Cl;	//Maximum resonance energy for a silicon-chlorine bond
	public double Min_en_Si_Cl;	//Minimum e-n attraction for a silicon-chlorine bond
	public double Max_en_Si_Cl;	//Maximum e-n attraction for a silicon-chlorine bond
	public double Min_nn_Si_Cl;	//Minimum n-n repulsion for a silicon-chlorine bond
	public double Max_nn_Si_Cl;	//Maximum n-n repulsion for a silicon-chlorine bond
	public double Min_eee_Si_Cl;	//Minimum exch. + e-e repulsion energy for a silicon-chlorine bond
	public double Max_eee_Si_Cl;	//Maximum exch. + e-e repulsion energy for a silicon-chlorine bond
	public double Min_tie_Si_Cl;	//Minimum tot. interaction energy for a silicon-chlorine bond
	public double Max_tie_Si_Cl;	//Maximum tot. interaction energy for a silicon-chlorine bond
	public double Min_re_Si_Br;	//Minimum resonance energy for a silicon-bromine bond
	public double Max_re_Si_Br;	//Maximum resonance energy for a silicon-bromine bond
	public double Min_en_Si_Br;	//Minimum e-n attraction for a silicon-bromine bond
	public double Max_en_Si_Br;	//Maximum e-n attraction for a silicon-bromine bond
	public double Min_nn_Si_Br;	//Minimum n-n repulsion for a silicon-bromine bond
	public double Max_nn_Si_Br;	//Maximum n-n repulsion for a silicon-bromine bond
	public double Min_eee_Si_Br;	//Minimum exch. + e-e repulsion energy for a silicon-bromine bond
	public double Max_eee_Si_Br;	//Maximum exch. + e-e repulsion energy for a silicon-bromine bond
	public double Min_tie_Si_Br;	//Minimum tot. interaction energy for a silicon-bromine bond
	public double Max_tie_Si_Br;	//Maximum tot. interaction energy for a silicon-bromine bond
	public double Min_re_Si_I;	//Minimum resonance energy for a silicon-iodine bond
	public double Max_re_Si_I;	//Maximum resonance energy for a silicon-iodine bond
	public double Min_en_Si_I;	//Minimum e-n attraction for a silicon-iodine bond
	public double Max_en_Si_I;	//Maximum e-n attraction for a silicon-iodine bond
	public double Min_nn_Si_I;	//Minimum n-n repulsion for a silicon-iodine bond
	public double Max_nn_Si_I;	//Maximum n-n repulsion for a silicon-iodine bond
	public double Min_eee_Si_I;	//Minimum exch. + e-e repulsion energy for a silicon-iodine bond
	public double Max_eee_Si_I;	//Maximum exch. + e-e repulsion energy for a silicon-iodine bond
	public double Min_tie_Si_I;	//Minimum tot. interaction energy for a silicon-iodine bond
	public double Max_tie_Si_I;	//Maximum tot. interaction energy for a silicon-iodine bond
	public double Min_re_Hg_Hg;	//Minimum resonance energy for a mercury-mercury bond
	public double Max_re_Hg_Hg;	//Maximum resonance energy for a mercury-mercury bond
	public double Min_en_Hg_Hg;	//Minimum e-n attraction for a mercury-mercury bond
	public double Max_en_Hg_Hg;	//Maximum e-n attraction for a mercury-mercury bond
	public double Min_nn_Hg_Hg;	//Minimum n-n repulsion for a mercury-mercury bond
	public double Max_nn_Hg_Hg;	//Maximum n-n repulsion for a mercury-mercury bond
	public double Min_eee_Hg_Hg;	//Minimum exch. + e-e repulsion energy for a mercury-mercury bond
	public double Max_eee_Hg_Hg;	//Maximum exch. + e-e repulsion energy for a mercury-mercury bond
	public double Min_tie_Hg_Hg;	//Minimum tot. interaction energy for a mercury-mercury bond
	public double Max_tie_Hg_Hg;	//Maximum tot. interaction energy for a mercury-mercury bond
	public double Min_re_Hg_H;	//Minimum resonance energy for a mercury-hydrogen bond
	public double Max_re_Hg_H;	//Maximum resonance energy for a mercury-hydrogen bond
	public double Min_en_Hg_H;	//Minimum e-n attraction for a mercury-hydrogen bond
	public double Max_en_Hg_H;	//Maximum e-n attraction for a mercury-hydrogen bond
	public double Min_nn_Hg_H;	//Minimum n-n repulsion for a mercury-hydrogen bond
	public double Max_nn_Hg_H;	//Maximum n-n repulsion for a mercury-hydrogen bond
	public double Min_eee_Hg_H;	//Minimum exch. + e-e repulsion energy for a mercury-hydrogen bond
	public double Max_eee_Hg_H;	//Maximum exch. + e-e repulsion energy for a mercury-hydrogen bond
	public double Min_tie_Hg_H;	//Minimum tot. interaction energy for a mercury-hydrogen bond
	public double Max_tie_Hg_H;	//Maximum tot. interaction energy for a mercury-hydrogen bond
	public double Min_re_Hg_F;	//Minimum resonance energy for a mercury-fluorine bond
	public double Max_re_Hg_F;	//Maximum resonance energy for a mercury-fluorine bond
	public double Min_en_Hg_F;	//Minimum e-n attraction for a mercury-fluorine bond
	public double Max_en_Hg_F;	//Maximum e-n attraction for a mercury-fluorine bond
	public double Min_nn_Hg_F;	//Minimum n-n repulsion for a mercury-fluorine bond
	public double Max_nn_Hg_F;	//Maximum n-n repulsion for a mercury-fluorine bond
	public double Min_eee_Hg_F;	//Minimum exch. + e-e repulsion energy for a mercury-fluorine bond
	public double Max_eee_Hg_F;	//Maximum exch. + e-e repulsion energy for a mercury-fluorine bond
	public double Min_tie_Hg_F;	//Minimum tot. interaction energy for a mercury-fluorine bond
	public double Max_tie_Hg_F;	//Maximum tot. interaction energy for a mercury-fluorine bond
	public double Min_re_Hg_Cl;	//Minimum resonance energy for a mercury-chlorine bond
	public double Max_re_Hg_Cl;	//Maximum resonance energy for a mercury-chlorine bond
	public double Min_en_Hg_Cl;	//Minimum e-n attraction for a mercury-chlorine bond
	public double Max_en_Hg_Cl;	//Maximum e-n attraction for a mercury-chlorine bond
	public double Min_nn_Hg_Cl;	//Minimum n-n repulsion for a mercury-chlorine bond
	public double Max_nn_Hg_Cl;	//Maximum n-n repulsion for a mercury-chlorine bond
	public double Min_eee_Hg_Cl;	//Minimum exch. + e-e repulsion energy for a mercury-chlorine bond
	public double Max_eee_Hg_Cl;	//Maximum exch. + e-e repulsion energy for a mercury-chlorine bond
	public double Min_tie_Hg_Cl;	//Minimum tot. interaction energy for a mercury-chlorine bond
	public double Max_tie_Hg_Cl;	//Maximum tot. interaction energy for a mercury-chlorine bond
	public double Min_re_Hg_Br;	//Minimum resonance energy for a mercury-bromine bond
	public double Max_re_Hg_Br;	//Maximum resonance energy for a mercury-bromine bond
	public double Min_en_Hg_Br;	//Minimum e-n attraction for a mercury-bromine bond
	public double Max_en_Hg_Br;	//Maximum e-n attraction for a mercury-bromine bond
	public double Min_nn_Hg_Br;	//Minimum n-n repulsion for a mercury-bromine bond
	public double Max_nn_Hg_Br;	//Maximum n-n repulsion for a mercury-bromine bond
	public double Min_eee_Hg_Br;	//Minimum exch. + e-e repulsion energy for a mercury-bromine bond
	public double Max_eee_Hg_Br;	//Maximum exch. + e-e repulsion energy for a mercury-bromine bond
	public double Min_tie_Hg_Br;	//Minimum tot. interaction energy for a mercury-bromine bond
	public double Max_tie_Hg_Br;	//Maximum tot. interaction energy for a mercury-bromine bond
	public double Min_re_Hg_I;	//Minimum resonance energy for a mercury-iodine bond
	public double Max_re_Hg_I;	//Maximum resonance energy for a mercury-iodine bond
	public double Min_en_Hg_I;	//Minimum e-n attraction for a mercury-iodine bond
	public double Max_en_Hg_I;	//Maximum e-n attraction for a mercury-iodine bond
	public double Min_nn_Hg_I;	//Minimum n-n repulsion for a mercury-iodine bond
	public double Max_nn_Hg_I;	//Maximum n-n repulsion for a mercury-iodine bond
	public double Min_eee_Hg_I;	//Minimum exch. + e-e repulsion energy for a mercury-iodine bond
	public double Max_eee_Hg_I;	//Maximum exch. + e-e repulsion energy for a mercury-iodine bond
	public double Min_tie_Hg_I;	//Minimum tot. interaction energy for a mercury-iodine bond
	public double Max_tie_Hg_I;	//Maximum tot. interaction energy for a mercury-iodine bond
	public double Min_re_H_H;	//Minimum resonance energy for a hydrogen-hydrogen bond
	public double Max_re_H_H;	//Maximum resonance energy for a hydrogen-hydrogen bond
	public double Min_en_H_H;	//Minimum e-n attraction for a hydrogen-hydrogen bond
	public double Max_en_H_H;	//Maximum e-n attraction for a hydrogen-hydrogen bond
	public double Min_nn_H_H;	//Minimum n-n repulsion for a hydrogen-hydrogen bond
	public double Max_nn_H_H;	//Maximum n-n repulsion for a hydrogen-hydrogen bond
	public double Min_eee_H_H;	//Minimum exch. + e-e repulsion energy for a hydrogen-hydrogen bond
	public double Max_eee_H_H;	//Maximum exch. + e-e repulsion energy for a hydrogen-hydrogen bond
	public double Min_tie_H_H;	//Minimum tot. interaction energy for a hydrogen-hydrogen bond
	public double Max_tie_H_H;	//Maximum tot. interaction energy for a hydrogen-hydrogen bond
	public double Min_re_H_F;	//Minimum resonance energy for a hydrogen-fluorine bond
	public double Max_re_H_F;	//Maximum resonance energy for a hydrogen-fluorine bond
	public double Min_en_H_F;	//Minimum e-n attraction for a hydrogen-fluorine bond
	public double Max_en_H_F;	//Maximum e-n attraction for a hydrogen-fluorine bond
	public double Min_nn_H_F;	//Minimum n-n repulsion for a hydrogen-fluorine bond
	public double Max_nn_H_F;	//Maximum n-n repulsion for a hydrogen-fluorine bond
	public double Min_eee_H_F;	//Minimum exch. + e-e repulsion energy for a hydrogen-fluorine bond
	public double Max_eee_H_F;	//Maximum exch. + e-e repulsion energy for a hydrogen-fluorine bond
	public double Min_tie_H_F;	//Minimum tot. interaction energy for a hydrogen-fluorine bond
	public double Max_tie_H_F;	//Maximum tot. interaction energy for a hydrogen-fluorine bond
	public double Min_re_H_Cl;	//Minimum resonance energy for a hydrogen-chlorine bond
	public double Max_re_H_Cl;	//Maximum resonance energy for a hydrogen-chlorine bond
	public double Min_en_H_Cl;	//Minimum e-n attraction for a hydrogen-chlorine bond
	public double Max_en_H_Cl;	//Maximum e-n attraction for a hydrogen-chlorine bond
	public double Min_nn_H_Cl;	//Minimum n-n repulsion for a hydrogen-chlorine bond
	public double Max_nn_H_Cl;	//Maximum n-n repulsion for a hydrogen-chlorine bond
	public double Min_eee_H_Cl;	//Minimum exch. + e-e repulsion energy for a hydrogen-chlorine bond
	public double Max_eee_H_Cl;	//Maximum exch. + e-e repulsion energy for a hydrogen-chlorine bond
	public double Min_tie_H_Cl;	//Minimum tot. interaction energy for a hydrogen-chlorine bond
	public double Max_tie_H_Cl;	//Maximum tot. interaction energy for a hydrogen-chlorine bond
	public double Min_re_H_Br;	//Minimum resonance energy for a hydrogen-bromine bond
	public double Max_re_H_Br;	//Maximum resonance energy for a hydrogen-bromine bond
	public double Min_en_H_Br;	//Minimum e-n attraction for a hydrogen-bromine bond
	public double Max_en_H_Br;	//Maximum e-n attraction for a hydrogen-bromine bond
	public double Min_nn_H_Br;	//Minimum n-n repulsion for a hydrogen-bromine bond
	public double Max_nn_H_Br;	//Maximum n-n repulsion for a hydrogen-bromine bond
	public double Min_eee_H_Br;	//Minimum exch. + e-e repulsion energy for a hydrogen-bromine bond
	public double Max_eee_H_Br;	//Maximum exch. + e-e repulsion energy for a hydrogen-bromine bond
	public double Min_tie_H_Br;	//Minimum tot. interaction energy for a hydrogen-bromine bond
	public double Max_tie_H_Br;	//Maximum tot. interaction energy for a hydrogen-bromine bond
	public double Min_re_H_I;	//Minimum resonance energy for a hydrogen-iodine bond
	public double Max_re_H_I;	//Maximum resonance energy for a hydrogen-iodine bond
	public double Min_en_H_I;	//Minimum e-n attraction for a hydrogen-iodine bond
	public double Max_en_H_I;	//Maximum e-n attraction for a hydrogen-iodine bond
	public double Min_nn_H_I;	//Minimum n-n repulsion for a hydrogen-iodine bond
	public double Max_nn_H_I;	//Maximum n-n repulsion for a hydrogen-iodine bond
	public double Min_eee_H_I;	//Minimum exch. + e-e repulsion energy for a hydrogen-iodine bond
	public double Max_eee_H_I;	//Maximum exch. + e-e repulsion energy for a hydrogen-iodine bond
	public double Min_tie_H_I;	//Minimum tot. interaction energy for a hydrogen-iodine bond
	public double Max_tie_H_I;	//Maximum tot. interaction energy for a hydrogen-iodine bond
	public double Min_re_F_F;	//Minimum resonance energy for a fluorine-fluorine bond
	public double Max_re_F_F;	//Maximum resonance energy for a fluorine-fluorine bond
	public double Min_en_F_F;	//Minimum e-n attraction for a fluorine-fluorine bond
	public double Max_en_F_F;	//Maximum e-n attraction for a fluorine-fluorine bond
	public double Min_nn_F_F;	//Minimum n-n repulsion for a fluorine-fluorine bond
	public double Max_nn_F_F;	//Maximum n-n repulsion for a fluorine-fluorine bond
	public double Min_eee_F_F;	//Minimum exch. + e-e repulsion energy for a fluorine-fluorine bond
	public double Max_eee_F_F;	//Maximum exch. + e-e repulsion energy for a fluorine-fluorine bond
	public double Min_tie_F_F;	//Minimum tot. interaction energy for a fluorine-fluorine bond
	public double Max_tie_F_F;	//Maximum tot. interaction energy for a fluorine-fluorine bond
	public double Min_re_F_Cl;	//Minimum resonance energy for a fluorine-chlorine bond
	public double Max_re_F_Cl;	//Maximum resonance energy for a fluorine-chlorine bond
	public double Min_en_F_Cl;	//Minimum e-n attraction for a fluorine-chlorine bond
	public double Max_en_F_Cl;	//Maximum e-n attraction for a fluorine-chlorine bond
	public double Min_nn_F_Cl;	//Minimum n-n repulsion for a fluorine-chlorine bond
	public double Max_nn_F_Cl;	//Maximum n-n repulsion for a fluorine-chlorine bond
	public double Min_eee_F_Cl;	//Minimum exch. + e-e repulsion energy for a fluorine-chlorine bond
	public double Max_eee_F_Cl;	//Maximum exch. + e-e repulsion energy for a fluorine-chlorine bond
	public double Min_tie_F_Cl;	//Minimum tot. interaction energy for a fluorine-chlorine bond
	public double Max_tie_F_Cl;	//Maximum tot. interaction energy for a fluorine-chlorine bond
	public double Min_re_F_Br;	//Minimum resonance energy for a fluorine-bromine bond
	public double Max_re_F_Br;	//Maximum resonance energy for a fluorine-bromine bond
	public double Min_en_F_Br;	//Minimum e-n attraction for a fluorine-bromine bond
	public double Max_en_F_Br;	//Maximum e-n attraction for a fluorine-bromine bond
	public double Min_nn_F_Br;	//Minimum n-n repulsion for a fluorine-bromine bond
	public double Max_nn_F_Br;	//Maximum n-n repulsion for a fluorine-bromine bond
	public double Min_eee_F_Br;	//Minimum exch. + e-e repulsion energy for a fluorine-bromine bond
	public double Max_eee_F_Br;	//Maximum exch. + e-e repulsion energy for a fluorine-bromine bond
	public double Min_tie_F_Br;	//Minimum tot. interaction energy for a fluorine-bromine bond
	public double Max_tie_F_Br;	//Maximum tot. interaction energy for a fluorine-bromine bond
	public double Min_re_F_I;	//Minimum resonance energy for a fluorine-iodine bond
	public double Max_re_F_I;	//Maximum resonance energy for a fluorine-iodine bond
	public double Min_en_F_I;	//Minimum e-n attraction for a fluorine-iodine bond
	public double Max_en_F_I;	//Maximum e-n attraction for a fluorine-iodine bond
	public double Min_nn_F_I;	//Minimum n-n repulsion for a fluorine-iodine bond
	public double Max_nn_F_I;	//Maximum n-n repulsion for a fluorine-iodine bond
	public double Min_eee_F_I;	//Minimum exch. + e-e repulsion energy for a fluorine-iodine bond
	public double Max_eee_F_I;	//Maximum exch. + e-e repulsion energy for a fluorine-iodine bond
	public double Min_tie_F_I;	//Minimum tot. interaction energy for a fluorine-iodine bond
	public double Max_tie_F_I;	//Maximum tot. interaction energy for a fluorine-iodine bond
	public double Min_re_Cl_Cl;	//Minimum resonance energy for a chlorine-chlorine bond
	public double Max_re_Cl_Cl;	//Maximum resonance energy for a chlorine-chlorine bond
	public double Min_en_Cl_Cl;	//Minimum e-n attraction for a chlorine-chlorine bond
	public double Max_en_Cl_Cl;	//Maximum e-n attraction for a chlorine-chlorine bond
	public double Min_nn_Cl_Cl;	//Minimum n-n repulsion for a chlorine-chlorine bond
	public double Max_nn_Cl_Cl;	//Maximum n-n repulsion for a chlorine-chlorine bond
	public double Min_eee_Cl_Cl;	//Minimum exch. + e-e repulsion energy for a chlorine-chlorine bond
	public double Max_eee_Cl_Cl;	//Maximum exch. + e-e repulsion energy for a chlorine-chlorine bond
	public double Min_tie_Cl_Cl;	//Minimum tot. interaction energy for a chlorine-chlorine bond
	public double Max_tie_Cl_Cl;	//Maximum tot. interaction energy for a chlorine-chlorine bond
	public double Min_re_Cl_Br;	//Minimum resonance energy for a chlorine-bromine bond
	public double Max_re_Cl_Br;	//Maximum resonance energy for a chlorine-bromine bond
	public double Min_en_Cl_Br;	//Minimum e-n attraction for a chlorine-bromine bond
	public double Max_en_Cl_Br;	//Maximum e-n attraction for a chlorine-bromine bond
	public double Min_nn_Cl_Br;	//Minimum n-n repulsion for a chlorine-bromine bond
	public double Max_nn_Cl_Br;	//Maximum n-n repulsion for a chlorine-bromine bond
	public double Min_eee_Cl_Br;	//Minimum exch. + e-e repulsion energy for a chlorine-bromine bond
	public double Max_eee_Cl_Br;	//Maximum exch. + e-e repulsion energy for a chlorine-bromine bond
	public double Min_tie_Cl_Br;	//Minimum tot. interaction energy for a chlorine-bromine bond
	public double Max_tie_Cl_Br;	//Maximum tot. interaction energy for a chlorine-bromine bond
	public double Min_re_Cl_I;	//Minimum resonance energy for a chlorine-iodine bond
	public double Max_re_Cl_I;	//Maximum resonance energy for a chlorine-iodine bond
	public double Min_en_Cl_I;	//Minimum e-n attraction for a chlorine-iodine bond
	public double Max_en_Cl_I;	//Maximum e-n attraction for a chlorine-iodine bond
	public double Min_nn_Cl_I;	//Minimum n-n repulsion for a chlorine-iodine bond
	public double Max_nn_Cl_I;	//Maximum n-n repulsion for a chlorine-iodine bond
	public double Min_eee_Cl_I;	//Minimum exch. + e-e repulsion energy for a chlorine-iodine bond
	public double Max_eee_Cl_I;	//Maximum exch. + e-e repulsion energy for a chlorine-iodine bond
	public double Min_tie_Cl_I;	//Minimum tot. interaction energy for a chlorine-iodine bond
	public double Max_tie_Cl_I;	//Maximum tot. interaction energy for a chlorine-iodine bond
	public double Min_re_Br_Br;	//Minimum resonance energy for a bromine-bromine bond
	public double Max_re_Br_Br;	//Maximum resonance energy for a bromine-bromine bond
	public double Min_en_Br_Br;	//Minimum e-n attraction for a bromine-bromine bond
	public double Max_en_Br_Br;	//Maximum e-n attraction for a bromine-bromine bond
	public double Min_nn_Br_Br;	//Minimum n-n repulsion for a bromine-bromine bond
	public double Max_nn_Br_Br;	//Maximum n-n repulsion for a bromine-bromine bond
	public double Min_eee_Br_Br;	//Minimum exch. + e-e repulsion energy for a bromine-bromine bond
	public double Max_eee_Br_Br;	//Maximum exch. + e-e repulsion energy for a bromine-bromine bond
	public double Min_tie_Br_Br;	//Minimum tot. interaction energy for a bromine-bromine bond
	public double Max_tie_Br_Br;	//Maximum tot. interaction energy for a bromine-bromine bond
	public double Min_re_Br_I;	//Minimum resonance energy for a bromine-iodine bond
	public double Max_re_Br_I;	//Maximum resonance energy for a bromine-iodine bond
	public double Min_en_Br_I;	//Minimum e-n attraction for a bromine-iodine bond
	public double Max_en_Br_I;	//Maximum e-n attraction for a bromine-iodine bond
	public double Min_nn_Br_I;	//Minimum n-n repulsion for a bromine-iodine bond
	public double Max_nn_Br_I;	//Maximum n-n repulsion for a bromine-iodine bond
	public double Min_eee_Br_I;	//Minimum exch. + e-e repulsion energy for a bromine-iodine bond
	public double Max_eee_Br_I;	//Maximum exch. + e-e repulsion energy for a bromine-iodine bond
	public double Min_tie_Br_I;	//Minimum tot. interaction energy for a bromine-iodine bond
	public double Max_tie_Br_I;	//Maximum tot. interaction energy for a bromine-iodine bond
	public double Min_re_I_I;	//Minimum resonance energy for an iodine-iodine bond
	public double Max_re_I_I;	//Maximum resonance energy for an iodine-iodine bond
	public double Min_en_I_I;	//Minimum e-n attraction for an iodine-iodine bond
	public double Max_en_I_I;	//Maximum e-n attraction for an iodine-iodine bond
	public double Min_nn_I_I;	//Minimum n-n repulsion for an iodine-iodine bond
	public double Max_nn_I_I;	//Maximum n-n repulsion for an iodine-iodine bond
	public double Min_eee_I_I;	//Minimum exch. + e-e repulsion energy for an iodine-iodine bond
	public double Max_eee_I_I;	//Maximum exch. + e-e repulsion energy for an iodine-iodine bond
	public double Min_tie_I_I;	//Minimum tot. interaction energy for an iodine-iodine bond
	public double Max_tie_I_I;	//Maximum tot. interaction energy for an iodine-iodine bond

	//1-center and 2-center energies
	public double Tot_1c_en;	//Tot. molecular 1-center e-n attraction energy
	public double Av_1c_en;		//Average molecular 1-center e-n attraction energy
	public double Tot_1c_ee;	//Total molecular 1-center e-e repulsion energy
	public double Av_1c_ee;		//Average molecular 1-center e-e repulsion energy
	public double Tot_2c_res;	//Total molecular 2-center resonance energy
	public double Av_2c_res;	//Average molecular 2-center resonance energy
	public double Tot_2c_exch;	//Total molecular 2-center exchange energy
	public double Av_2c_exch;	//Average molecular 2-center exchange energy
	public double Tot_2c_en;	//Total molecular 2-center e-n attraction energy
	public double Av_2c_en;		//Average molecular 2-center e-n attraction energy
	public double Tot_2c_ee;	//Total molecular 2-center e-e repulsion energy
	public double Av_2c_ee;		//Average molecular 2-center e-e repulsion energy

	//Charged Partial surface area (CPSA) descriptors
	public double pnsa1;		//Partial negative surface area
	public double ppsa1;		//Partial positive surface area
	public double pnsa2;		//Total charge weighted negative surface area
	public double ppsa2;		//Total charge weighted positive surface area
	public double pnsa3;		//Atomic charge weighted negative surface area
	public double ppsa3;		//Atomic charge weighted positive surface area
	public double dpsa1;		//Difference in charged partial surface area
	public double dpsa2;		//Difference in total charge weighted surface area
	public double dpsa3;		//Difference in atomic charge weighted surface area
	public double fnsa1;		//Fractional partial negative surface area
	public double fnsa2;		//Fractional total charge weighted negative surface area
	public double fnsa3;		//Fractional atomic charge weighted negative surface area
	public double fpsa1;		//Fractional partial positive surface area
	public double fpsa2;		//Fractional total charge weighted positive surface area
	public double fpsa3;		//Fractional atomic charge weighted positive surface area
	public double wnsa1;		//Surface weighted partial negative surface area
	public double wnsa2;		//Surface weighted total charge weighted negative surface area
	public double wnsa3;		//Surface weighted atomic charge weighted negative surface area
	public double wpsa1;		//Surface weighted partial positive surface area
	public double wpsa2;		//Surface weighted total charge weighted positive surface area
	public double wpsa3;		//Surface weighted atomic charge weighted positive surface area
	public double rncg;			//Relative negative charge
	public double rpcg;			//Relative positive charge
	public double rncs;			//Relative negative charged surface area
	public double rpcs;			//Relative positive charged surface area
	public double tasa;			//Total hydrophobic surface area
	public double tpsa;			//Total polar surface area
	public double rasa;			//Relative hydrophobic surface area
	public double rpsa;			//Relative polar surface area

	//Hydrogen bond Charged Partial Surface Area (HB-CPSA) descriptors
	public double nhba;			//Number of hydrogen bond acceptors
	public double nhbd;			//Number of hydrogen bond donors
	public double rhta;			//Ratio of number of donor groups to acceptor groups
	public double hdsa;			//Sum of surface areas of hydrogens that can be donated
	public double rsah;			//Average surface area of hydrogens that can be donated
	public double fhdsa;		//Fraction of total molecular surface area associated with hydrogens
	public double hasa;			//Sum of surface areas of all hydrogen bond acceptor atoms
	public double hasa2;		//Variant of hasa as defined by Katritzsky et al.
	public double rsaa;			//Average surface area of hydrogen bond acceptor groups
	public double fhasa;		//Fraction of total molecular surface area associated with hydrogen bonds
	public double fhasa2;		//Variant of hasa2 as defined by Katritzsky et al.
	public double hdca;			//Sum of charged surface areas of hydrogens that can be donated
	public double csa2h;		//Charged surface areas of hydrogen atoms
	public double csa2cl;		//Charged surface areas of chlorine atoms
	public double fhdca;		//Charged surface area of hydrogens that can be donated relative to TMSA
	public double hdca2;		//Variant of fhdca as defined by Katritzsky et al.
	public double hdsa2;		//Variant of hdca2 as defined by Katritzsky et al.
	public double haca;			//Sum of charged surface areas of hydrogen bond acceptors
	public double scaa2;		//Average hydrogen bond acceptor charged surface area
	public double fhaca;		//Charged surface area of hydrogen bond acceptors relative to TMSA
	public double hbsa;			//Sum of surface areas of both hydrogen bond acceptors and donors
	public double fhbsa;		//Surface areas of hydrogen bond donors and acceptors relative to TMSA
	public double hbca;			//Sum of charged surface areas of both hydrogen bond acceptors and donors
	public double fhbca;		//Charged surface areas of both hydrogen bond acceptors and donors relative to TMSA
	public double chaa1;		//Sum of partial charges on hydrogen bond acceptor atoms
	public double chaa2;		//Average value of chaa1
	public double acgd;			//Average difference in charge between all pairs of hydrogen bond donors
	public double hrpcg;		//Relative positive charge (RPCG) restricted to hydrogen bond donors
	public double hrncg;		//Relative negative charge (RNCG) restricted to hydrogen bond acceptors
	public double hrpcs;		//Relative positive charge surface area (RPCS) restricted to hydrogen bond donors
	public double hrncs;		//Relative negative charge surface area (RNCS) restricted to hydrogen bond acceptors
	public double chgd;			//Max. difference in charge between a hydrogen bond donor and acceptor

	//Surface area and volume descriptors
	public double vdWSA;		//van der Waals surface area
	public double vdWV;			//van der Waals volume
	public double SASA;			//Solvent accessible surface area
	public double SAV;			//Solvent accessible volume
	public double cesa;			//Contact re-entrant surface area
	public double cev;			//Contact re-entrant volume
	public double ovality;		//Ovality

	//Frontier molecular orbital descriptors
	public double HOMO_2;		//Energy of the HOMO-2 orbital
	public double HOMO_1;		//Energy of the HOMO-1 orbital
	public double HOMO;			//Energy of the HOMO orbital
	public double LUMO;			//Energy of the LUMO orbital
	public double LUMO_1;		//Energy of the LUMO+1 orbital
	public double LUMO_2;		//Energy of the LUMO+2 orbital
	public double aoep_HOMO_2;	//Atomic orbital electronic population of HOMO-2 orbital
	public double aoep_HOMO_1;	//Atomic orbital electronic population of HOMO-1 orbital
	public double aoep_HOMO;	//Atomic orbital electronic population of HOMO orbital
	public double aoep_LUMO;	//Atomic orbital electronic population of LUMO orbital
	public double aoep_LUMO_1;	//Atomic orbital electronic population of LUMO+1 orbital
	public double aoep_LUMO_2;	//Atomic orbital electronic population of LUMO+2 orbital
	public double hardness;		//Hardness index
	public double HOMO_LUMO;	//HOMO - LUMO energy gap
	public double HOMLUM;		//HOMO/LUMO energy fraction
	public double elect;		//Atom electronegativity
	public double ip;			//Ionization potential

	//Other descriptors
//	public double mw;			//Molecular weight
	public double Dipole_charge;	//Total point charge component of the dipole moment
	public double Dipole_hybrid;	//Total hybrid component of the dipole moment
	public double Dipole_X;			//Dipole moment in the X axis
	public double Dipole_Y;			//Dipole moment in the Y axis
	public double Dipole_Z;			//Dipole moment in the Z axis
	public double Dipole;			//Dipole moment
	public double PMI_A;			//Principal Moment of Inertia A
	public double PMI_B;			//Principal Moment of Inertia B
	public double PMI_C;			//Principal Moment of Inertia C
	public double apolar;			//Alpha polarizability
	public double bpolar;			//Beta polarizability
	public double gpolar;			//Gamma polarizability
	public double eccentricity;		//Molecular eccentricity
	public double asphericity;		//Asphericity
	public double sphericity;		//Sphericity


	// ***********************************************************************
	// convenience string arrays for accessing fields using arrays:

	public String[] varlist2d = { "strChi", "strKappa", "strES", "strES_acnt","strHES",
			"strESMaxMin", "str2D","strIC", "strMDEC", "strBurden", "strTopo",
			"strCon", "str2DAuto","strWalk","strMP" }; // array of list names

//	public String [] vlist2d = {"strIC"};

//	String[] varlist3d = {"strWHIM","strGetaway","strMOMI","strRandicShape", "strQM"};
	String[] varlist3d = {"strWHIM","strGetaway","strMOMI","strRandicShape"};

	public String[] strxp = { "x0", "x1", "x2", "xp3", "xp4", "xp5", "xp6",
			"xp7", "xp8", "xp9", "xp10" };

	public String[] strxvp = { "xv0", "xv1", "xv2", "xvp3", "xvp4", "xvp5",
			"xvp6", "xvp7", "xvp8", "xvp9", "xvp10" };

	public String[] strChi = { "x0", "x1", "x2", "xp3", "xp4", "xp5", "xp6",
			"xp7", "xp8", "xp9", "xp10", "xc3", "xc4", "xpc4", "xch3", "xch4",
			"xch5", "xch6", "xch7", "xch8", "xch9", "xch10", "knotp", "xv0",
			"xv1", "xv2", "xvp3", "xvp4", "xvp5", "xvp6", "xvp7", "xvp8",
			"xvp9", "xvp10", "xvc3", "xvc4", "xvpc4", "xvch3", "xvch4",
			"xvch5", "xvch6", "xvch7", "xvch8", "xvch9", "xvch10", "knotpv" };

//	public String [] strKappa={"k0", "k1", "k2", "k3","ka0", "ka1", "ka2", "ka3","phia"};
	public String [] strKappa={"k1", "k2", "k3","ka1", "ka2", "ka3","phia"};


	public String[] strES={"SsCH3","SdCH2","SssCH2","StCH","SdsCH","SaaCH","SsssCH","SddC","StsC","SdssC","SsaaC","SaaaC","SssssC","SsNH2","SdNH","SssNH","SaaNH","StN","SdsN","SaaN","SsssN","SaaaN","SsaaN","SsNH3p","SssNH2p","SsssNHp","SddNp","StsNp","SdNm","SdssNp","SssssNp","SsaaNp","SsOH","SdO","SssO","SaaO","SsOm","SsF","SsCl","SsBr","SsI","SsSiH3","SssSiH2","SsssSiH","SssssSi","SsssP","SdsssP","SddsP","SsssssP","SdssPH","SsPH2","SssPH","SsSH","SdS","SssS","SaaS","SdssS","SddssS","SssssssS","SsSnH3","SssSnH2","SsssSnH","SssssSn","SsPbH3","SssPbH2","SsssPbH","SssssPb","SsssAs","SdsssAs","SddsAs","SsssssAs","SsAsH2","SssAsH","SsHgp","SssHg","SsHg","StCm"};

//	public String[] strES_acnt={"SsCH3_acnt","SdCH2_acnt","SssCH2_acnt","StCH_acnt","SdsCH_acnt","SaaCH_acnt","SsssCH_acnt","SddC_acnt","StsC_acnt","SdssC_acnt","SsaaC_acnt","SaaaC_acnt","SssssC_acnt","SsNH2_acnt","SdNH_acnt","SssNH_acnt","SaaNH_acnt","StN_acnt","SdsN_acnt","SaaN_acnt","SsssN_acnt","SaaaN_acnt","SsaaN_acnt","SsNH3p_acnt","SssNH2p_acnt","SsssNHp_acnt","SddNp_acnt","StsNp_acnt","SdNm_acnt","SdssNp_acnt","SssssNp_acnt","SsaaNp_acnt","SsOH_acnt","SdO_acnt","SssO_acnt","SaaO_acnt","SsOm_acnt","SsF_acnt","SsCl_acnt","SsBr_acnt","SsI_acnt","SsSiH3_acnt","SssSiH2_acnt","SsssSiH_acnt","SssssSi_acnt","SsssP_acnt","SdsssP_acnt","SddsP_acnt","SsssssP_acnt","SdssPH_acnt","SsPH2_acnt","SssPH_acnt","SsSH_acnt","SdS_acnt","SssS_acnt","SaaS_acnt","SdssS_acnt","SddssS_acnt","SssssssS_acnt","SsSnH3_acnt","SssSnH2_acnt","SsssSnH_acnt","SssssSn_acnt","SsPbH3_acnt","SssPbH2_acnt","SsssPbH_acnt","SssssPb_acnt","SsssAs_acnt","SdsssAs_acnt","SddsAs_acnt","SsssssAs_acnt","SsAsH2_acnt","SssAsH_acnt","SsHgp_acnt","SssHg_acnt","SsHg_acnt","StCm_acnt"};
	public String[] strES_acnt=null;//now generated in constructor


	public String[] strHES = { "SHsOH", "SHother", "SHdNH", "SHsSH", "SHsNH2",
			"SHssNH", "SHtCH", "SHCHnX" };

	public String [] strESMaxMin={"Hmax","Gmax","Hmin","Gmin","Hmaxpos","Hminneg"};



	public String[] str2D = { "numHBa", "numHBd", "numwHBd",
			"SHHBa", "SHHBd", "SHwHBd", "Qs", "Qv", "Qsv" };


	public String[] strIC = { "ic","si","ib","I","maxic","ssi","R","eim","iadje", "iadjem", "iadjm", "iadjmm", "ivdem",
			"ivdmm", "ieadje", "ieadjem", "ieadjm", "ieadjmm", "ide", "idem",
			 "idm", "idmm", "iddem", "iddmm", "iede", "iedem", "We",
			"iedm", "iedmm", "tvc", "icyce", "icycem", "icycm", "icycmm",
			"tti", "ttvi" };
	public String[] strMP={"XLOGP","XLOGP2","ALOGP","ALOGP2","AMR","Hy","Ui"};


	public String[] strMDEC = { "MDEC11", "MDEC12", "MDEC13", "MDEC14",
			"MDEC22", "MDEC23", "MDEC24", "MDEC33", "MDEC34", "MDEC44",
			"MDEO11", "MDEO12", "MDEO22", "MDEN11", "MDEN12", "MDEN13",
			"MDEN22", "MDEN23", "MDEN33" };

	public String[] strBurden = { "BEHm1", "BEHm2", "BEHm3", "BEHm4", "BEHm5",
			"BEHm6", "BEHm7", "BEHm8", "BELm1", "BELm2", "BELm3", "BELm4",
			"BELm5", "BELm6", "BELm7", "BELm8", "BEHv1", "BEHv2", "BEHv3",
			"BEHv4", "BEHv5", "BEHv6", "BEHv7", "BEHv8", "BELv1", "BELv2",
			"BELv3", "BELv4", "BELv5", "BELv6", "BELv7", "BELv8", "BEHe1",
			"BEHe2", "BEHe3", "BEHe4", "BEHe5", "BEHe6", "BEHe7", "BEHe8",
			"BELe1", "BELe2", "BELe3", "BELe4", "BELe5", "BELe6", "BELe7",
			"BELe8", "BEHp1", "BEHp2", "BEHp3", "BEHp4", "BEHp5", "BEHp6",
			"BEHp7", "BEHp8", "BELp1", "BELp2", "BELp3", "BELp4", "BELp5",
			"BELp6", "BELp7", "BELp8" };


	public String[] strGetaway = { "HGM", "ITH", "ISH", "HIC", "H0u", "H1u",
			"H2u", "H3u", "H4u", "H5u", "H6u", "H7u", "H8u", "HTu", "HATS0u",
			"HATS1u", "HATS2u", "HATS3u", "HATS4u", "HATS5u", "HATS6u",
			"HATS7u", "HATS8u", "HATSu", "H0m", "H1m", "H2m", "H3m", "H4m",
			"H5m", "H6m", "H7m", "H8m", "HTm", "HATS0m", "HATS1m", "HATS2m",
			"HATS3m", "HATS4m", "HATS5m", "HATS6m", "HATS7m", "HATS8m",
			"HATSm", "H0v", "H1v", "H2v", "H3v", "H4v", "H5v", "H6v", "H7v",
			"H8v", "HTv", "HATS0v", "HATS1v", "HATS2v", "HATS3v", "HATS4v",
			"HATS5v", "HATS6v", "HATS7v", "HATS8v", "HATSv", "H0e", "H1e",
			"H2e", "H3e", "H4e", "H5e", "H6e", "H7e", "H8e", "HTe", "HATS0e",
			"HATS1e", "HATS2e", "HATS3e", "HATS4e", "HATS5e", "HATS6e",
			"HATS7e", "HATS8e", "HATSe", "H0p", "H1p", "H2p", "H3p", "H4p",
			"H5p", "H6p", "H7p", "H8p", "HTp", "HATS0p", "HATS1p", "HATS2p",
			"HATS3p", "HATS4p", "HATS5p", "HATS6p", "HATS7p", "HATS8p",
			"HATSp", "RCON", "RARS", "REIG", "R1u", "R2u", "R3u", "R4u", "R5u",
			"R6u", "R7u", "R8u", "RTu", "R1uplus", "R2uplus", "R3uplus",
			"R4uplus", "R5uplus", "R6uplus", "R7uplus", "R8uplus", "RTuplus",
			"R1m", "R2m", "R3m", "R4m", "R5m", "R6m", "R7m", "R8m", "RTm",
			"R1mplus", "R2mplus", "R3mplus", "R4mplus", "R5mplus", "R6mplus",
			"R7mplus", "R8mplus", "RTmplus", "R1v", "R2v", "R3v", "R4v", "R5v",
			"R6v", "R7v", "R8v", "RTv", "R1vplus", "R2vplus", "R3vplus",
			"R4vplus", "R5vplus", "R6vplus", "R7vplus", "R8vplus", "RTvplus",
			"R1e", "R2e", "R3e", "R4e", "R5e", "R6e", "R7e", "R8e", "RTe",
			"R1eplus", "R2eplus", "R3eplus", "R4eplus", "R5eplus", "R6eplus",
			"R7eplus", "R8eplus", "RTeplus", "R1p", "R2p", "R3p", "R4p", "R5p",
			"R6p", "R7p", "R8p", "RTp", "R1pplus", "R2pplus", "R3pplus",
			"R4pplus", "R5pplus", "R6pplus", "R7pplus", "R8pplus", "RTpplus" };


	public String[] strTopo = { "ZM1", "ZM1V", "ZM2", "ZM2V", "J", "Jt", "BAC",
			"Lop", "ICR", "TIE", "MAXDN", "MAXDP", "DELS", "W", "WA" };


	public String[] str2DAuto = { "ATS1m", "ATS2m", "ATS3m", "ATS4m", "ATS5m",
			"ATS6m", "ATS7m", "ATS8m", "ATS1v", "ATS2v", "ATS3v", "ATS4v",
			"ATS5v", "ATS6v", "ATS7v", "ATS8v", "ATS1e", "ATS2e", "ATS3e",
			"ATS4e", "ATS5e", "ATS6e", "ATS7e", "ATS8e", "ATS1p", "ATS2p",
			"ATS3p", "ATS4p", "ATS5p", "ATS6p", "ATS7p", "ATS8p", "MATS1m",
			"MATS2m", "MATS3m", "MATS4m", "MATS5m", "MATS6m", "MATS7m",
			"MATS8m", "MATS1v", "MATS2v", "MATS3v", "MATS4v", "MATS5v",
			"MATS6v", "MATS7v", "MATS8v", "MATS1e", "MATS2e", "MATS3e",
			"MATS4e", "MATS5e", "MATS6e", "MATS7e", "MATS8e", "MATS1p",
			"MATS2p", "MATS3p", "MATS4p", "MATS5p", "MATS6p", "MATS7p",
			"MATS8p", "GATS1m", "GATS2m", "GATS3m", "GATS4m", "GATS5m",
			"GATS6m", "GATS7m", "GATS8m", "GATS1v", "GATS2v", "GATS3v",
			"GATS4v", "GATS5v", "GATS6v", "GATS7v", "GATS8v", "GATS1e",
			"GATS2e", "GATS3e", "GATS4e", "GATS5e", "GATS6e", "GATS7e",
			"GATS8e", "GATS1p", "GATS2p", "GATS3p", "GATS4p", "GATS5p",
			"GATS6p", "GATS7p", "GATS8p" };




	public String[] strCon = { "MW", "AMW", "Sv", "Se", "Sp", "Ss", "Mv", "Me",
			"Mp", "Ms", "nAT", "nSK", "nBT", "nBO", "nBM", "SCBO", "ARR",
			"nCIC", "nCIR", "nDB", "nTB", "nAB", "nH", "nC",
			"nN", "nO", "nP", "nS", "nF", "nCL", "nBR", "nI", "nB",
			"nX", "nR03", "nR04", "nR05", "nR06", "nR07", "nR08", "nR09",
			"nR10", "nR11", "nR12", "nBnz" };

	public String[] strWHIM = { "L1u", "L2u", "L3u", "P1u", "P2u", "G1u",
			"G2u", "G3u", "E1u", "E2u", "E3u", "L1m", "L2m", "L3m", "P1m",
			"P2m", "G1m", "G2m", "G3m", "E1m", "E2m", "E3m", "L1v", "L2v",
			"L3v", "P1v", "P2v", "G1v", "G2v", "G3v", "E1v", "E2v", "E3v",
			"L1e", "L2e", "L3e", "P1e", "P2e", "G1e", "G2e", "G3e", "E1e",
			"E2e", "E3e", "L1p", "L2p", "L3p", "P1p", "P2p", "G1p", "G2p",
			"G3p", "E1p", "E2p", "E3p", "L1s", "L2s", "L3s", "P1s", "P2s",
			"G1s", "G2s", "G3s", "E1s", "E2s", "E3s", "Tu", "Tm", "Tv", "Te",
			"Tp", "Ts", "Au", "Am", "Av", "Ae", "Ap", "As", "Vu", "Vm", "Vv",
			"Ve", "Vp", "Vs", "Gu", "Gm", "Gs", "Ku", "Km", "Kv", "Ke", "Kp",
			"Ks", "Du", "Dm", "Dv", "De", "Dp", "Ds" };


	public String[] strWalk = { "MWC01", "MWC02", "MWC03", "MWC04", "MWC05",
			"MWC06", "MWC07", "MWC08", "MWC09", "MWC10", "SRW01", "SRW02",
			"SRW03", "SRW04", "SRW05", "SRW06", "SRW07", "SRW08", "SRW09",
			"SRW10", "MPC01", "MPC02", "MPC03", "MPC04", "MPC05", "MPC06",
			"MPC07", "MPC08", "MPC09", "MPC10", "piPC01", "piPC02", "piPC03",
			"piPC04", "piPC05", "piPC06", "piPC07", "piPC08", "piPC09",
			"piPC10", "piID", "TWC", "TPC", "CID", "CID2", "BID" };


	public String [] strMOMI={"MOMI1","MOMI2","MOMI3","MOMI4","MOMI5","MOMI6","MOMI7"};

	public String[] strRandicShape = { "DP01", "DP02", "DP03", "DP04", "DP05",
			"DP06", "DP07", "DP08", "DP09", "DP10", "DP11", "DP12", "DP13",
			"DP14", "DP15", "DP16", "DP17", "DP18", "DP19", "DP20", "SP01",
			"SP02", "SP03", "SP04", "SP05", "SP06", "SP07", "SP08", "SP09",
			"SP10", "SP11", "SP12", "SP13", "SP14", "SP15", "SP16", "SP17",
			"SP18", "SP19", "SP20","SHP2" };


	public String [] strQM = {"qpmax", "qnmax", "Qpos", "Qneg", "Qtot", "Qmean", "Q2",
			"SPP", "SPP2", "Pprime", "DP", "TE1", "TE2", "PCWTE", "LDI", "ET0", "ET1",
			"ET2", "ET3", "MaxPC_C", "MinPC_C", "MaxPC_N", "MinPC_N", "MaxPC_P",
			"MinPC_P", "MaxPC_O", "MinPC_O", "MaxPC_S", "MinPC_S", "MaxPC_Si", "MinPC_Si",
			"MaxPC_Hg", "MinPC_Hg", "MaxPC_H", "MinPC_H", "MaxPC_F", "MinPC_F", "MaxPC_Cl",
			"MinPC_Cl", "MaxPC_Br", "MinPC_Br", "MaxPC_I", "MinPC_I", "Max_EFD_C",
			"Min_EFD_C", "Tot_EFD_C", "Av_EFD_C", "Max_EFD_N", "Min_EFD_N", "Tot_EFD_N",
			"Av_EFD_N", "Max_EFD_P", "Min_EFD_P", "Tot_EFD_P", "Av_EFD_P", "Max_EFD_O",
			"Min_EFD_O", "Tot_EFD_O", "Av_EFD_O", "Max_EFD_S", "Min_EFD_S", "Tot_EFD_S",
			"Av_EFD_S", "Max_EFD_Si", "Min_EFD_Si", "Tot_EFD_Si", "Av_EFD_Si", "Max_EFD_Hg",
			"Min_EFD_Hg", "Tot_EFD_Hg", "Av_EFD_Hg", "Max_EFD_H", "Min_EFD_H", "Tot_EFD_H",
			"Av_EFD_H", "Max_EFD_F", "Min_EFD_F", "Tot_EFD_F", "Av_EFD_F", "Max_EFD_Cl",
			"Min_EFD_Cl", "Tot_EFD_Cl", "Av_EFD_Cl", "Max_EFD_Br", "Min_EFD_Br", "Tot_EFD_Br",
			"Av_EFD_Br", "Max_EFD_I", "Min_EFD_I", "Tot_EFD_I", "Av_EFD_I", "Max_EFD",
			"Min_EFD", "Tot_EFD", "Av_EFD", "Max_NFD_C", "Min_NFD_C", "Tot_NFD_C", "Av_NFD_C",
			"Max_NFD_N", "Min_NFD_N", "Tot_NFD_N", "Av_NFD_N", "Max_NFD_P", "Min_NFD_P",
			"Tot_NFD_P", "Av_NFD_P", "Max_NFD_O", "Min_NFD_O", "Tot_NFD_O", "Av_NFD_O",
			"Max_NFD_S", "Min_NFD_S", "Tot_NFD_S", "Av_NFD_S", "Max_NFD_Si", "Min_NFD_Si",
			"Tot_NFD_Si", "Av_NFD_Si", "Max_NFD_Hg", "Min_NFD_Hg", "Tot_NFD_Hg", "Av_NFD_Hg",
			"Max_NFD_H", "Min_NFD_H", "Tot_NFD_H", "Av_NFD_H", "Max_NFD_F", "Min_NFD_F",
			"Tot_NFD_F", "Av_NFD_F", "Max_NFD_Cl", "Min_NFD_Cl", "Tot_NFD_Cl", "Av_NFD_Cl",
			"Max_NFD_Br", "Min_NFD_Br", "Tot_NFD_Br", "Av_NFD_Br", "Max_NFD_I", "Min_NFD_I",
			"Tot_NFD_I", "Av_NFD_I", "Max_NFD", "Min_NFD", "Tot_NFD", "Av_NFD", "Max_RFD_C",
			"Min_RFD_C", "Tot_RFD_C", "Av_RFD_C", "Max_RFD_N", "Min_RFD_N", "Tot_RFD_N",
			"Av_RFD_N", "Max_RFD_P", "Min_RFD_P", "Tot_RFD_P", "Av_RFD_P", "Max_RFD_O",
			"Min_RFD_O", "Tot_RFD_O", "Av_RFD_O", "Max_RFD_S", "Min_RFD_S", "Tot_RFD_S",
			"Av_RFD_S", "Max_RFD_Si", "Min_RFD_Si", "Tot_RFD_Si", "Av_RFD_Si", "Max_RFD_Hg",
			"Min_RFD_Hg", "Tot_RFD_Hg", "Av_RFD_Hg", "Max_RFD_H", "Min_RFD_H", "Tot_RFD_H",
			"Av_RFD_H", "Max_RFD_F", "Min_RFD_F", "Tot_RFD_F", "Av_RFD_F", "Max_RFD_Cl",
			"Min_RFD_Cl", "Tot_RFD_Cl", "Av_RFD_Cl", "Max_RFD_Br", "Min_RFD_Br", "Tot_RFD_Br",
			"Av_RFD_Br", "Max_RFD_I", "Min_RFD_I", "Tot_RFD_I", "Av_RFD_I", "Max_RFD",
			"Min_RFD", "Tot_RFD", "Av_RFD", "Max_EFDI_C", "Min_EFDI_C", "Tot_EFDI_C",
			"Av_EFDI_C", "Max_EFDI_N", "Min_EFDI_N", "Tot_EFDI_N", "Av_EFDI_N", "Max_EFDI_P",
			"Min_EFDI_P", "Tot_EFDI_P", "Av_EFDI_P", "Max_EFDI_O", "Min_EFDI_O", "Tot_EFDI_O",
			"Av_EFDI_O", "Max_EFDI_S", "Min_EFDI_S", "Tot_EFDI_S", "Av_EFDI_S", "Max_EFDI_Si",
			"Min_EFDI_Si", "Tot_EFDI_Si", "Av_EFDI_Si", "Max_EFDI_Hg", "Min_EFDI_Hg",
			"Tot_EFDI_Hg", "Av_EFDI_Hg", "Max_EFDI_H", "Min_EFDI_H", "Tot_EFDI_H", "Av_EFDI_H",
			"Max_EFDI_F", "Min_EFDI_F", "Tot_EFDI_F", "Av_EFDI_F", "Max_EFDI_Cl",
			"Min_EFDI_Cl", "Tot_EFDI_Cl", "Av_EFDI_Cl", "Max_EFDI_Br", "Min_EFDI_Br",
			"Tot_EFDI_Br", "Av_EFDI_Br", "Max_EFDI_I", "Min_EFDI_I", "Tot_EFDI_I", "Av_EFDI_I",
			"Max_EFDI", "Min_EFDI", "Tot_EFDI", "Av_EFDI", "Max_NFDI_C", "Min_NFDI_C",
			"Tot_NFDI_C", "Av_NFDI_C", "Max_NFDI_N", "Min_NFDI_N", "Tot_NFDI_N", "Av_NFDI_N",
			"Max_NFDI_P", "Min_NFDI_P", "Tot_NFDI_P", "Av_NFDI_P", "Max_NFDI_O", "Min_NFDI_O",
			"Tot_NFDI_O", "Av_NFDI_O", "Max_NFDI_S", "Min_NFDI_S", "Tot_NFDI_S", "Av_NFDI_S",
			"Max_NFDI_Si", "Min_NFDI_Si", "Tot_NFDI_Si", "Av_NFDI_Si", "Max_NFDI_Hg",
			"Min_NFDI_Hg", "Tot_NFDI_Hg", "Av_NFDI_Hg", "Max_NFDI_H", "Min_NFDI_H",
			"Tot_NFDI_H", "Av_NFDI_H", "Max_NFDI_F", "Min_NFDI_F", "Tot_NFDI_F", "Av_NFDI_F",
			"Max_NFDI_Cl", "Min_NFDI_Cl", "Tot_NFDI_Cl", "Av_NFDI_Cl", "Max_NFDI_Br",
			"Min_NFDI_Br", "Tot_NFDI_Br", "Av_NFDI_Br", "Max_NFDI_I", "Min_NFDI_I",
			"Tot_NFDI_I", "Av_NFDI_I", "Max_NFDI", "Min_NFDI", "Tot_NFDI", "Av_NFDI",
			"Max_RFDI_C", "Min_RFDI_C", "Tot_RFDI_C", "Av_RFDI_C", "Max_RFDI_N", "Min_RFDI_N",
			"Tot_RFDI_N", "Av_RFDI_N", "Max_RFDI_P", "Min_RFDI_P", "Tot_RFDI_P", "Av_RFDI_P",
			"Max_RFDI_O", "Min_RFDI_O", "Tot_RFDI_O", "Av_RFDI_O", "Max_RFDI_S", "Min_RFDI_S",
			"Tot_RFDI_S", "Av_RFDI_S", "Max_RFDI_Si", "Min_RFDI_Si", "Tot_RFDI_Si",
			"Av_RFDI_Si", "Max_RFDI_Hg", "Min_RFDI_Hg", "Tot_RFDI_Hg", "Av_RFDI_Hg",
			"Max_RFDI_H", "Min_RFDI_H", "Tot_RFDI_H", "Av_RFDI_H", "Max_RFDI_F", "Min_RFDI_F",
			"Tot_RFDI_F", "Av_RFDI_F", "Max_RFDI_Cl", "Min_RFDI_Cl", "Tot_RFDI_Cl", "Av_RFDI_Cl",
			"Max_RFDI_Br", "Min_RFDI_Br", "Tot_RFDI_Br", "Av_RFDI_Br", "Max_RFDI_I",
			"Min_RFDI_I", "Tot_RFDI_I", "Av_RFDI_I", "Max_RFDI", "Min_RFDI", "Tot_RFDI",
			"Av_RFDI", "Max_ESD_C", "Min_ESD_C", "Tot_ESD_C", "Av_ESD_C", "Max_ESD_N",
			"Min_ESD_N", "Tot_ESD_N", "Av_ESD_N", "Max_ESD_P", "Min_ESD_P", "Tot_ESD_P",
			"Av_ESD_P", "Max_ESD_O", "Min_ESD_O", "Tot_ESD_O", "Av_ESD_O", "Max_ESD_S",
			"Min_ESD_S", "Tot_ESD_S", "Av_ESD_S", "Max_ESD_Si", "Min_ESD_Si", "Tot_ESD_Si",
			"Av_ESD_Si", "Max_ESD_Hg", "Min_ESD_Hg", "Tot_ESD_Hg", "Av_ESD_Hg", "Max_ESD_H",
			"Min_ESD_H", "Tot_ESD_H", "Av_ESD_H", "Max_ESD_F", "Min_ESD_F", "Tot_ESD_F",
			"Av_ESD_F", "Max_ESD_Cl", "Min_ESD_Cl", "Tot_ESD_Cl", "Av_ESD_Cl", "Max_ESD_Br",
			"Min_ESD_Br", "Tot_ESD_Br", "Av_ESD_Br", "Max_ESD_I", "Min_ESD_I", "Tot_ESD_I",
			"Av_ESD_I", "Max_ESD", "Min_ESD", "Tot_ESD", "Av_ESD", "Max_NSD_C", "Min_NSD_C",
			"Tot_NSD_C", "Av_NSD_C", "Max_NSD_N", "Min_NSD_N", "Tot_NSD_N", "Av_NSD_N",
			"Max_NSD_P", "Min_NSD_P", "Tot_NSD_P", "Av_NSD_P", "Max_NSD_O", "Min_NSD_O",
			"Tot_NSD_O", "Av_NSD_O", "Max_NSD_S", "Min_NSD_S", "Tot_NSD_S", "Av_NSD_S",
			"Max_NSD_Si", "Min_NSD_Si", "Tot_NSD_Si", "Av_NSD_Si", "Max_NSD_Hg", "Min_NSD_Hg",
			"Tot_NSD_Hg", "Av_NSD_Hg", "Max_NSD_H", "Min_NSD_H", "Tot_NSD_H", "Av_NSD_H",
			"Max_NSD_F", "Min_NSD_F", "Tot_NSD_F", "Av_NSD_F", "Max_NSD_Cl", "Min_NSD_Cl",
			"Tot_NSD_Cl", "Av_NSD_Cl", "Max_NSD_Br", "Min_NSD_Br", "Tot_NSD_Br", "Av_NSD_Br",
			"Max_NSD_I", "Min_NSD_I", "Tot_NSD_I", "Av_NSD_I", "Max_NSD", "Min_NSD",
			"Tot_NSD", "Av_NSD", "Min_val_C", "Max_val_C", "Av_val_C", "Min_BO_C", "Max_BO_C",
			"Av_BO_C", "Min_val_N", "Max_val_N", "Av_val_N", "Min_BO_N", "Max_BO_N",
			"Av_BO_N", "Min_val_P", "Max_val_P", "Av_val_P", "Min_BO_P", "Max_BO_P",
			"Av_BO_P", "Min_val_O", "Max_val_O", "Av_val_O", "Min_BO_O", "Max_BO_O",
			"Av_BO_O", "Min_val_S", "Max_val_S", "Av_val_S", "Min_BO_S", "Max_BO_S",
			"Av_BO_S", "Min_val_Si", "Max_val_Si", "Av_val_Si", "Min_BO_Si", "Max_BO_Si",
			"Av_BO_Si", "Min_val_Hg", "Max_val_Hg", "Av_val_Hg", "Min_BO_Hg", "Max_BO_Hg",
			"Av_BO_Hg", "Min_val_H", "Max_val_H", "Av_val_H", "Min_BO_H", "Max_BO_H",
			"Av_BO_H", "Min_val_F", "Max_val_F", "Av_val_F", "Min_BO_F", "Max_BO_F",
			"Av_BO_F", "Min_val_Cl", "Max_val_Cl", "Av_val_Cl", "Min_BO_Cl", "Max_BO_Cl",
			"Av_BO_Cl", "Min_val_Br", "Max_val_Br", "Av_val_Br", "Min_BO_Br", "Max_BO_Br",
			"Av_BO_Br", "Min_val_I", "Max_val_I", "Av_val_I", "Min_BO_I", "Max_BO_I",
			"Av_BO_I", "Min_en_C", "Max_en_C", "Min_ee_C", "Max_ee_C", "Min_ase_C",
			"Max_ase_C", "Min_en_N", "Max_en_N", "Min_ee_N", "Max_ee_N", "Min_ase_N",
			"Max_ase_N", "Min_en_P", "Max_en_P", "Min_ee_P", "Max_ee_P", "Min_ase_P",
			"Max_ase_P", "Min_en_O", "Max_en_O", "Min_ee_O", "Max_ee_O", "Min_ase_O",
			"Max_ase_O", "Min_en_S", "Max_en_S", "Min_ee_S", "Max_ee_S", "Min_ase_S",
			"Max_ase_S", "Min_en_Si", "Max_en_Si", "Min_ee_Si", "Max_ee_Si", "Min_ase_Si",
			"Max_ase_Si", "Min_en_Hg", "Max_en_Hg", "Min_ee_Hg", "Max_ee_Hg", "Min_ase_Hg",
			"Max_ase_Hg", "Min_en_H", "Max_en_H", "Min_ee_H", "Max_ee_H", "Min_ase_H",
			"Max_ase_H", "Min_en_F", "Max_en_F", "Min_ee_F", "Max_ee_F", "Min_ase_F",
			"Max_ase_F", "Min_en_Cl", "Max_en_Cl", "Min_ee_Cl", "Max_ee_Cl", "Min_ase_Cl",
			"Max_ase_Cl", "Min_en_Br", "Max_en_Br", "Min_ee_Br", "Max_ee_Br", "Min_ase_Br",
			"Max_ase_Br", "Min_en_I", "Max_en_I", "Min_ee_I", "Max_ee_I", "Min_ase_I",
			"Max_ase_I", "Min_re_C_C", "Max_re_C_C", "Min_en_C_C", "Max_en_C_C", "Min_nn_C_C",
			"Max_nn_C_C", "Min_eee_C_C", "Max_eee_C_C", "Min_tie_C_C", "Max_tie_C_C",
			"Min_re_C_N", "Max_re_C_N", "Min_en_C_N", "Max_en_C_N", "Min_nn_C_N",
			"Max_nn_C_N", "Min_eee_C_N", "Max_eee_C_N", "Min_tie_C_N", "Max_tie_C_N",
			"Min_re_C_P", "Max_re_C_P", "Min_en_C_P", "Max_en_C_P", "Min_nn_C_P", "Max_nn_C_P",
			"Min_eee_C_P", "Max_eee_C_P", "Min_tie_C_P", "Max_tie_C_P", "Min_re_C_O",
			"Max_re_C_O", "Min_en_C_O", "Max_en_C_O", "Min_nn_C_O", "Max_nn_C_O",
			"Min_eee_C_O", "Max_eee_C_O", "Min_tie_C_O", "Max_tie_C_O", "Min_re_C_S",
			"Max_re_C_S", "Min_en_C_S", "Max_en_C_S", "Min_nn_C_S", "Max_nn_C_S",
			"Min_eee_C_S", "Max_eee_C_S", "Min_tie_C_S", "Max_tie_C_S", "Min_re_C_Si",
			"Max_re_C_Si", "Min_en_C_Si", "Max_en_C_Si", "Min_nn_C_Si", "Max_nn_C_Si",
			"Min_eee_C_Si", "Max_eee_C_Si", "Min_tie_C_Si", "Max_tie_C_Si", "Min_re_C_Hg",
			"Max_re_C_Hg", "Min_en_C_Hg", "Max_en_C_Hg", "Min_nn_C_Hg", "Max_nn_C_Hg",
			"Min_eee_C_Hg", "Max_eee_C_Hg", "Min_tie_C_Hg", "Max_tie_C_Hg", "Min_re_C_H",
			"Max_re_C_H", "Min_en_C_H", "Max_en_C_H", "Min_nn_C_H", "Max_nn_C_H",
			"Min_eee_C_H", "Max_eee_C_H", "Min_tie_C_H", "Max_tie_C_H", "Min_re_C_F",
			"Max_re_C_F", "Min_en_C_F", "Max_en_C_F", "Min_nn_C_F", "Max_nn_C_F",
			"Min_eee_C_F", "Max_eee_C_F", "Min_tie_C_F", "Max_tie_C_F", "Min_re_C_Cl",
			"Max_re_C_Cl", "Min_en_C_Cl", "Max_en_C_Cl", "Min_nn_C_Cl", "Max_nn_C_Cl",
			"Min_eee_C_Cl", "Max_eee_C_Cl", "Min_tie_C_Cl", "Max_tie_C_Cl", "Min_re_C_Br",
			"Max_re_C_Br", "Min_en_C_Br", "Max_en_C_Br", "Min_nn_C_Br", "Max_nn_C_Br",
			"Min_eee_C_Br", "Max_eee_C_Br", "Min_tie_C_Br", "Max_tie_C_Br", "Min_re_C_I",
			"Max_re_C_I", "Min_en_C_I", "Max_en_C_I", "Min_nn_C_I", "Max_nn_C_I",
			"Min_eee_C_I", "Max_eee_C_I", "Min_tie_C_I", "Max_tie_C_I", "Min_re_N_N",
			"Max_re_N_N", "Min_en_N_N", "Max_en_N_N", "Min_nn_N_N", "Max_nn_N_N",
			"Min_eee_N_N", "Max_eee_N_N", "Min_tie_N_N", "Max_tie_N_N", "Min_re_N_P",
			"Max_re_N_P", "Min_en_N_P", "Max_en_N_P", "Min_nn_N_P", "Max_nn_N_P",
			"Min_eee_N_P", "Max_eee_N_P", "Min_tie_N_P", "Max_tie_N_P", "Min_re_N_O",
			"Max_re_N_O", "Min_en_N_O", "Max_en_N_O", "Min_nn_N_O", "Max_nn_N_O",
			"Min_eee_N_O", "Max_eee_N_O", "Min_tie_N_O", "Max_tie_N_O", "Min_re_N_S",
			"Max_re_N_S", "Min_en_N_S", "Max_en_N_S", "Min_nn_N_S", "Max_nn_N_S",
			"Min_eee_N_S", "Max_eee_N_S", "Min_tie_N_S", "Max_tie_N_S", "Min_re_N_Si",
			"Max_re_N_Si", "Min_en_N_Si", "Max_en_N_Si", "Min_nn_N_Si", "Max_nn_N_Si",
			"Min_eee_N_Si", "Max_eee_N_Si", "Min_tie_N_Si", "Max_tie_N_Si", "Min_re_N_Hg",
			"Max_re_N_Hg", "Min_en_N_Hg", "Max_en_N_Hg", "Min_nn_N_Hg", "Max_nn_N_Hg",
			"Min_eee_N_Hg", "Max_eee_N_Hg", "Min_tie_N_Hg", "Max_tie_N_Hg", "Min_re_N_H",
			"Max_re_N_H", "Min_en_N_H", "Max_en_N_H", "Min_nn_N_H", "Max_nn_N_H",
			"Min_eee_N_H", "Max_eee_N_H", "Min_tie_N_H", "Max_tie_N_H", "Min_re_N_F",
			"Max_re_N_F", "Min_en_N_F", "Max_en_N_F", "Min_nn_N_F", "Max_nn_N_F",
			"Min_eee_N_F", "Max_eee_N_F", "Min_tie_N_F", "Max_tie_N_F", "Min_re_N_Cl",
			"Max_re_N_Cl", "Min_en_N_Cl", "Max_en_N_Cl", "Min_nn_N_Cl", "Max_nn_N_Cl",
			"Min_eee_N_Cl", "Max_eee_N_Cl", "Min_tie_N_Cl", "Max_tie_N_Cl", "Min_re_N_Br",
			"Max_re_N_Br", "Min_en_N_Br", "Max_en_N_Br", "Min_nn_N_Br", "Max_nn_N_Br",
			"Min_eee_N_Br", "Max_eee_N_Br", "Min_tie_N_Br", "Max_tie_N_Br", "Min_re_N_I",
			"Max_re_N_I", "Min_en_N_I", "Max_en_N_I", "Min_nn_N_I", "Max_nn_N_I",
			"Min_eee_N_I", "Max_eee_N_I", "Min_tie_N_I", "Max_tie_N_I", "Min_re_P_P",
			"Max_re_P_P", "Min_en_P_P", "Max_en_P_P", "Min_nn_P_P", "Max_nn_P_P",
			"Min_eee_P_P", "Max_eee_P_P", "Min_tie_P_P", "Max_tie_P_P", "Min_re_P_O",
			"Max_re_P_O", "Min_en_P_O", "Max_en_P_O", "Min_nn_P_O", "Max_nn_P_O",
			"Min_eee_P_O", "Max_eee_P_O", "Min_tie_P_O", "Max_tie_P_O", "Min_re_P_S",
			"Max_re_P_S", "Min_en_P_S", "Max_en_P_S", "Min_nn_P_S", "Max_nn_P_S",
			"Min_eee_P_S", "Max_eee_P_S", "Min_tie_P_S", "Max_tie_P_S", "Min_re_P_Si",
			"Max_re_P_Si", "Min_en_P_Si", "Max_en_P_Si", "Min_nn_P_Si", "Max_nn_P_Si",
			"Min_eee_P_Si", "Max_eee_P_Si", "Min_tie_P_Si", "Max_tie_P_Si", "Min_re_P_Hg",
			"Max_re_P_Hg", "Min_en_P_Hg", "Max_en_P_Hg", "Min_nn_P_Hg", "Max_nn_P_Hg",
			"Min_eee_P_Hg", "Max_eee_P_Hg", "Min_tie_P_Hg", "Max_tie_P_Hg", "Min_re_P_H",
			"Max_re_P_H", "Min_en_P_H", "Max_en_P_H", "Min_nn_P_H", "Max_nn_P_H",
			"Min_eee_P_H", "Max_eee_P_H", "Min_tie_P_H", "Max_tie_P_H", "Min_re_P_F",
			"Max_re_P_F", "Min_en_P_F", "Max_en_P_F", "Min_nn_P_F", "Max_nn_P_F",
			"Min_eee_P_F", "Max_eee_P_F", "Min_tie_P_F", "Max_tie_P_F", "Min_re_P_Cl",
			"Max_re_P_Cl", "Min_en_P_Cl", "Max_en_P_Cl", "Min_nn_P_Cl", "Max_nn_P_Cl",
			"Min_eee_P_Cl", "Max_eee_P_Cl", "Min_tie_P_Cl", "Max_tie_P_Cl", "Min_re_P_Br",
			"Max_re_P_Br", "Min_en_P_Br", "Max_en_P_Br", "Min_nn_P_Br", "Max_nn_P_Br",
			"Min_eee_P_Br", "Max_eee_P_Br", "Min_tie_P_Br", "Max_tie_P_Br", "Min_re_P_I",
			"Max_re_P_I", "Min_en_P_I", "Max_en_P_I", "Min_nn_P_I", "Max_nn_P_I",
			"Min_eee_P_I", "Max_eee_P_I", "Min_tie_P_I", "Max_tie_P_I", "Min_re_O_O",
			"Max_re_O_O", "Min_en_O_O", "Max_en_O_O", "Min_nn_O_O", "Max_nn_O_O",
			"Min_eee_O_O", "Max_eee_O_O", "Min_tie_O_O", "Max_tie_O_O", "Min_re_O_S",
			"Max_re_O_S", "Min_en_O_S", "Max_en_O_S", "Min_nn_O_S", "Max_nn_O_S",
			"Min_eee_O_S", "Max_eee_O_S", "Min_tie_O_S", "Max_tie_O_S", "Min_re_O_Si",
			"Max_re_O_Si", "Min_en_O_Si", "Max_en_O_Si", "Min_nn_O_Si", "Max_nn_O_Si",
			"Min_eee_O_Si", "Max_eee_O_Si", "Min_tie_O_Si", "Max_tie_O_Si", "Min_re_O_Hg",
			"Max_re_O_Hg", "Min_en_O_Hg", "Max_en_O_Hg", "Min_nn_O_Hg", "Max_nn_O_Hg",
			"Min_eee_O_Hg", "Max_eee_O_Hg", "Min_tie_O_Hg", "Max_tie_O_Hg", "Min_re_O_H",
			"Max_re_O_H", "Min_en_O_H", "Max_en_O_H", "Min_nn_O_H", "Max_nn_O_H",
			"Min_eee_O_H", "Max_eee_O_H", "Min_tie_O_H", "Max_tie_O_H", "Min_re_O_F",
			"Max_re_O_F", "Min_en_O_F", "Max_en_O_F", "Min_nn_O_F", "Max_nn_O_F",
			"Min_eee_O_F", "Max_eee_O_F", "Min_tie_O_F", "Max_tie_O_F", "Min_re_O_Cl",
			"Max_re_O_Cl", "Min_en_O_Cl", "Max_en_O_Cl", "Min_nn_O_Cl", "Max_nn_O_Cl",
			"Min_eee_O_Cl", "Max_eee_O_Cl", "Min_tie_O_Cl", "Max_tie_O_Cl", "Min_re_O_Br",
			"Max_re_O_Br", "Min_en_O_Br", "Max_en_O_Br", "Min_nn_O_Br", "Max_nn_O_Br",
			"Min_eee_O_Br", "Max_eee_O_Br", "Min_tie_O_Br", "Max_tie_O_Br", "Min_re_O_I",
			"Max_re_O_I", "Min_en_O_I", "Max_en_O_I", "Min_nn_O_I", "Max_nn_O_I",
			"Min_eee_O_I", "Max_eee_O_I", "Min_tie_O_I", "Max_tie_O_I", "Min_re_S_S",
			"Max_re_S_S", "Min_en_S_S", "Max_en_S_S", "Min_nn_S_S", "Max_nn_S_S",
			"Min_eee_S_S", "Max_eee_S_S", "Min_tie_S_S", "Max_tie_S_S", "Min_re_S_Si",
			"Max_re_S_Si", "Min_en_S_Si", "Max_en_S_Si", "Min_nn_S_Si", "Max_nn_S_Si",
			"Min_eee_S_Si", "Max_eee_S_Si", "Min_tie_S_Si", "Max_tie_S_Si", "Min_re_S_Hg",
			"Max_re_S_Hg", "Min_en_S_Hg", "Max_en_S_Hg", "Min_nn_S_Hg", "Max_nn_S_Hg",
			"Min_eee_S_Hg", "Max_eee_S_Hg", "Min_tie_S_Hg", "Max_tie_S_Hg", "Min_re_S_H",
			"Max_re_S_H", "Min_en_S_H", "Max_en_S_H", "Min_nn_S_H", "Max_nn_S_H",
			"Min_eee_S_H", "Max_eee_S_H", "Min_tie_S_H", "Max_tie_S_H", "Min_re_S_F",
			"Max_re_S_F", "Min_en_S_F", "Max_en_S_F", "Min_nn_S_F", "Max_nn_S_F",
			"Min_eee_S_F", "Max_eee_S_F", "Min_tie_S_F", "Max_tie_S_F", "Min_re_S_Cl",
			"Max_re_S_Cl", "Min_en_S_Cl", "Max_en_S_Cl", "Min_nn_S_Cl", "Max_nn_S_Cl",
			"Min_eee_S_Cl", "Max_eee_S_Cl", "Min_tie_S_Cl", "Max_tie_S_Cl", "Min_re_S_Br",
			"Max_re_S_Br", "Min_en_S_Br", "Max_en_S_Br", "Min_nn_S_Br", "Max_nn_S_Br",
			"Min_eee_S_Br", "Max_eee_S_Br", "Min_tie_S_Br", "Max_tie_S_Br", "Min_re_S_I",
			"Max_re_S_I", "Min_en_S_I", "Max_en_S_I", "Min_nn_S_I", "Max_nn_S_I",
			"Min_eee_S_I", "Max_eee_S_I", "Min_tie_S_I", "Max_tie_S_I", "Min_re_Si_Si",
			"Max_re_Si_Si", "Min_en_Si_Si", "Max_en_Si_Si", "Min_nn_Si_Si", "Max_nn_Si_Si",
			"Min_eee_Si_Si", "Max_eee_Si_Si", "Min_tie_Si_Si", "Max_tie_Si_Si", "Min_re_Si_Hg",
			"Max_re_Si_Hg", "Min_en_Si_Hg", "Max_en_Si_Hg", "Min_nn_Si_Hg", "Max_nn_Si_Hg",
			"Min_eee_Si_Hg", "Max_eee_Si_Hg", "Min_tie_Si_Hg", "Max_tie_Si_Hg", "Min_re_Si_H",
			"Max_re_Si_H", "Min_en_Si_H", "Max_en_Si_H", "Min_nn_Si_H", "Max_nn_Si_H",
			"Min_eee_Si_H", "Max_eee_Si_H", "Min_tie_Si_H", "Max_tie_Si_H", "Min_re_Si_F",
			"Max_re_Si_F", "Min_en_Si_F", "Max_en_Si_F", "Min_nn_Si_F", "Max_nn_Si_F",
			"Min_eee_Si_F", "Max_eee_Si_F", "Min_tie_Si_F", "Max_tie_Si_F", "Min_re_Si_Cl",
			"Max_re_Si_Cl", "Min_en_Si_Cl", "Max_en_Si_Cl", "Min_nn_Si_Cl", "Max_nn_Si_Cl",
			"Min_eee_Si_Cl", "Max_eee_Si_Cl", "Min_tie_Si_Cl", "Max_tie_Si_Cl",
			"Min_re_Si_Br", "Max_re_Si_Br", "Min_en_Si_Br", "Max_en_Si_Br", "Min_nn_Si_Br",
			"Max_nn_Si_Br", "Min_eee_Si_Br", "Max_eee_Si_Br", "Min_tie_Si_Br", "Max_tie_Si_Br",
			"Min_re_Si_I", "Max_re_Si_I", "Min_en_Si_I", "Max_en_Si_I", "Min_nn_Si_I",
			"Max_nn_Si_I", "Min_eee_Si_I", "Max_eee_Si_I", "Min_tie_Si_I", "Max_tie_Si_I",
			"Min_re_Hg_Hg", "Max_re_Hg_Hg", "Min_en_Hg_Hg", "Max_en_Hg_Hg", "Min_nn_Hg_Hg",
			"Max_nn_Hg_Hg", "Min_eee_Hg_Hg", "Max_eee_Hg_Hg", "Min_tie_Hg_Hg", "Max_tie_Hg_Hg",
			"Min_re_Hg_H", "Max_re_Hg_H", "Min_en_Hg_H", "Max_en_Hg_H", "Min_nn_Hg_H",
			"Max_nn_Hg_H", "Min_eee_Hg_H", "Max_eee_Hg_H", "Min_tie_Hg_H", "Max_tie_Hg_H",
			"Min_re_Hg_F", "Max_re_Hg_F", "Min_en_Hg_F", "Max_en_Hg_F", "Min_nn_Hg_F",
			"Max_nn_Hg_F", "Min_eee_Hg_F", "Max_eee_Hg_F", "Min_tie_Hg_F", "Max_tie_Hg_F",
			"Min_re_Hg_Cl", "Max_re_Hg_Cl", "Min_en_Hg_Cl", "Max_en_Hg_Cl", "Min_nn_Hg_Cl",
			"Max_nn_Hg_Cl", "Min_eee_Hg_Cl", "Max_eee_Hg_Cl", "Min_tie_Hg_Cl", "Max_tie_Hg_Cl",
			"Min_re_Hg_Br", "Max_re_Hg_Br", "Min_en_Hg_Br", "Max_en_Hg_Br", "Min_nn_Hg_Br",
			"Max_nn_Hg_Br", "Min_eee_Hg_Br", "Max_eee_Hg_Br", "Min_tie_Hg_Br", "Max_tie_Hg_Br",
			"Min_re_Hg_I", "Max_re_Hg_I", "Min_en_Hg_I", "Max_en_Hg_I", "Min_nn_Hg_I",
			"Max_nn_Hg_I", "Min_eee_Hg_I", "Max_eee_Hg_I", "Min_tie_Hg_I", "Max_tie_Hg_I",
			"Min_re_H_H", "Max_re_H_H", "Min_en_H_H", "Max_en_H_H", "Min_nn_H_H", "Max_nn_H_H",
			"Min_eee_H_H", "Max_eee_H_H", "Min_tie_H_H", "Max_tie_H_H", "Min_re_H_F",
			"Max_re_H_F", "Min_en_H_F", "Max_en_H_F", "Min_nn_H_F", "Max_nn_H_F", "Min_eee_H_F",
			"Max_eee_H_F", "Min_tie_H_F", "Max_tie_H_F", "Min_re_H_Cl", "Max_re_H_Cl",
			"Min_en_H_Cl", "Max_en_H_Cl", "Min_nn_H_Cl", "Max_nn_H_Cl", "Min_eee_H_Cl",
			"Max_eee_H_Cl", "Min_tie_H_Cl", "Max_tie_H_Cl", "Min_re_H_Br", "Max_re_H_Br",
			"Min_en_H_Br", "Max_en_H_Br", "Min_nn_H_Br", "Max_nn_H_Br", "Min_eee_H_Br",
			"Max_eee_H_Br", "Min_tie_H_Br", "Max_tie_H_Br", "Min_re_H_I", "Max_re_H_I",
			"Min_en_H_I", "Max_en_H_I", "Min_nn_H_I", "Max_nn_H_I", "Min_eee_H_I",
			"Max_eee_H_I", "Min_tie_H_I", "Max_tie_H_I", "Min_re_F_F", "Max_re_F_F",
			"Min_en_F_F", "Max_en_F_F", "Min_nn_F_F", "Max_nn_F_F", "Min_eee_F_F", "Max_eee_F_F",
			"Min_tie_F_F", "Max_tie_F_F", "Min_re_F_Cl", "Max_re_F_Cl", "Min_en_F_Cl",
			"Max_en_F_Cl", "Min_nn_F_Cl", "Max_nn_F_Cl", "Min_eee_F_Cl", "Max_eee_F_Cl",
			"Min_tie_F_Cl", "Max_tie_F_Cl", "Min_re_F_Br", "Max_re_F_Br", "Min_en_F_Br",
			"Max_en_F_Br", "Min_nn_F_Br", "Max_nn_F_Br", "Min_eee_F_Br", "Max_eee_F_Br",
			"Min_tie_F_Br", "Max_tie_F_Br", "Min_re_F_I", "Max_re_F_I", "Min_en_F_I",
			"Max_en_F_I", "Min_nn_F_I", "Max_nn_F_I", "Min_eee_F_I", "Max_eee_F_I",
			"Min_tie_F_I", "Max_tie_F_I", "Min_re_Cl_Cl", "Max_re_Cl_Cl", "Min_en_Cl_Cl",
			"Max_en_Cl_Cl", "Min_nn_Cl_Cl", "Max_nn_Cl_Cl", "Min_eee_Cl_Cl", "Max_eee_Cl_Cl",
			"Min_tie_Cl_Cl", "Max_tie_Cl_Cl", "Min_re_Cl_Br", "Max_re_Cl_Br", "Min_en_Cl_Br",
			"Max_en_Cl_Br", "Min_nn_Cl_Br", "Max_nn_Cl_Br", "Min_eee_Cl_Br", "Max_eee_Cl_Br",
			"Min_tie_Cl_Br", "Max_tie_Cl_Br", "Min_re_Cl_I", "Max_re_Cl_I", "Min_en_Cl_I",
			"Max_en_Cl_I", "Min_nn_Cl_I", "Max_nn_Cl_I", "Min_eee_Cl_I", "Max_eee_Cl_I",
			"Min_tie_Cl_I", "Max_tie_Cl_I", "Min_re_Br_Br", "Max_re_Br_Br", "Min_en_Br_Br",
			"Max_en_Br_Br", "Min_nn_Br_Br", "Max_nn_Br_Br", "Min_eee_Br_Br", "Max_eee_Br_Br",
			"Min_tie_Br_Br", "Max_tie_Br_Br", "Min_re_Br_I", "Max_re_Br_I", "Min_en_Br_I",
			"Max_en_Br_I", "Min_nn_Br_I", "Max_nn_Br_I", "Min_eee_Br_I", "Max_eee_Br_I",
			"Min_tie_Br_I", "Max_tie_Br_I", "Min_re_I_I", "Max_re_I_I", "Min_en_I_I",
			"Max_en_I_I", "Min_nn_I_I", "Max_nn_I_I", "Min_eee_I_I", "Max_eee_I_I",
			"Min_tie_I_I", "Max_tie_I_I", "Tot_1c_en", "Av_1c_en", "Tot_1c_ee", "Av_1c_ee",
			"Tot_2c_res", "Av_2c_res", "Tot_2c_exch", "Av_2c_exch", "Tot_2c_en", "Av_2c_en",
			"Tot_2c_ee", "Av_2c_ee", "pnsa1", "ppsa1", "pnsa2", "ppsa2", "pnsa3", "ppsa3",
			"dpsa1", "dpsa2", "dpsa3", "fnsa1", "fnsa2", "fnsa3", "fpsa1", "fpsa2", "fpsa3",
			"wnsa1", "wnsa2", "wnsa3", "wpsa1", "wpsa2", "wpsa3", "rncg", "rpcg", "rncs",
			"rpcs", "tasa", "tpsa", "rasa", "rpsa", "nhba", "nhbd", "rhta", "hdsa", "rsah",
			"fhdsa", "hasa", "hasa2", "rsaa", "fhasa", "fhasa2", "hdca", "csa2h", "csa2cl",
			"fhdca", "hdca2", "hdsa2", "haca", "scaa2", "fhaca", "hbsa", "fhbsa", "hbca",
			"fhbca", "chaa1", "chaa2", "acgd", "hrpcg", "hrncg", "hrpcs", "hrncs", "chgd",
			"vdWSA", "vdWV", "SASA", "SAV", "cesa", "cev", "ovality", "HOMO_2", "HOMO_1",
			"HOMO", "LUMO", "LUMO_1", "LUMO_2", "aoep_HOMO_2", "aoep_HOMO_1", "aoep_HOMO",
			"aoep_LUMO", "aoep_LUMO_1", "aoep_LUMO_2", "hardness", "HOMO_LUMO", "HOMLUM",
			"elect", "ip", "Dipole_charge", "Dipole_hybrid", "Dipole_X", "Dipole_Y",
			"Dipole_Z", "Dipole", "PMI_A", "PMI_B", "PMI_C", "apolar", "bpolar", "gpolar",
			"eccentricity", "asphericity", "sphericity"};


	public Hashtable<String,String> htVarDef;


	public void Initialize() {
		try {
			Object obj=this;
			Field fields []=this.getClass().getFields();

			for (int i=0;i<=fields.length-1;i++) {
				String type = fields[i].getType().toString();

				if (type.equals("int")) {
					fields[i].setInt(obj, -9999);
				} else if (type.equals("double")) {
					fields[i].setDouble(obj, -9999);
				}

			}

		} catch (Exception e) {

		}

	}


	/**
	 * Writes out descriptors to readable file for looking at all descriptors in
	 * a single vertical table
	 * 
	 * @param m
	 * @param dd
	 */
	public Vector<String> getDescriptorValues(String Delimiter) {
		Vector<String>vals=new Vector<>();
		
		try {
			
			DecimalFormat myD8 = new DecimalFormat("0.########");
			String[] varlist = varlist2d;

			
			for (int i = 0; i <= varlist.length - 1; i++) {
				Field myField = this.getClass().getField(varlist[i]);
				String[] names = (String[]) myField.get(this);
				for (int j = 0; j <= names.length - 1; j++) {
					Field myField2 = this.getClass().getField(names[j]);

					String val = myD8.format(myField2.getDouble(this));

					vals.add(names[j] + Delimiter + val);
				}
			}

			for (int i = 0; i < strFragments.length; i++) {
				String strVar = strFragments[i];
				double Val = (Double) FragmentList.get(strVar);
				vals.add(strVar + Delimiter + Val);
			}
			

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return vals;
	}
	
	/**
	 * Store as LinkedHashMap so values are in the same order as we put them in
	 * Also great for exporting later to csv
	 * 	 * 
	 * @return
	 */
	public LinkedHashMap<String,String> convertToLinkedHashMap() {

		LinkedHashMap<String,String> lhm=new LinkedHashMap<>();
		
		try {
						
			if (Error.contentEquals("OK"))
				lhm.put("Error","");
			else 
				lhm.put("Error",Error);
			
			if (!Error.contentEquals("OK")) {
				return lhm;
			}
						
			DecimalFormat myD8 = new DecimalFormat("0.########");
			DecimalFormat myD = new DecimalFormat("0");

			String[] varlist = varlist2d;

			for (int i = 0; i <= varlist.length - 1; i++) {
				Field myField = this.getClass().getField(varlist[i]);
				String[] names = (String[]) myField.get(this);
				for (int j = 0; j <= names.length - 1; j++) {
					Field myField2 = this.getClass().getField(names[j]);
					String val = myD8.format(myField2.getDouble(this));
					lhm.put(names[j],val);
				}
			}

			for (int i = 0; i < strFragments.length; i++) {
				String strVar = strFragments[i];
				String val = myD.format(FragmentList.get(strVar));
				lhm.put(strVar,val);
			}


		} catch (Exception e) {
			e.printStackTrace();
		}
		return lhm;
	}
	
	
	public Vector<String>getDescriptorValues() {

		Vector <String>values=new Vector<>();
		
		try {
			
			if (!Error.contentEquals("OK")) {
				values.add("Error: "+Error.replace("\n", "").replace("\r", ""));
				return values;
			}
							
			DecimalFormat myD8 = new DecimalFormat("0.########");
			DecimalFormat myD = new DecimalFormat("0");

			String[] varlist = varlist2d;

			for (int i = 0; i <= varlist.length - 1; i++) {
				Field myField = this.getClass().getField(varlist[i]);
				String[] names = (String[]) myField.get(this);
				for (int j = 0; j <= names.length - 1; j++) {
					Field myField2 = this.getClass().getField(names[j]);
					String val = myD8.format(myField2.getDouble(this));

					values.add(val);
				}
			}

			for (int i = 0; i < strFragments.length; i++) {
				String strVar = strFragments[i];
				String val = myD.format(FragmentList.get(strVar));
				values.add(val);
			}


		} catch (Exception e) {
			e.printStackTrace();
		}
		return values;
	}
	
	void load3DDefs() {

		try {

			String file = "3d-descriptor-list.csv";
			InputStream ins = this.getClass().getClassLoader()
					.getResourceAsStream(file);
			InputStreamReader isr = new InputStreamReader(ins);

			// create ArrayLists to store the descriptor number and name of
			// descriptor
			ArrayList<Integer> num = new ArrayList<Integer>(); // List to store
																// the
																// descriptor
																// number
			ArrayList<String> dname = new ArrayList<String>(); // List to store
																// the
																// descriptor
																// name
			ArrayList<String> desc = new ArrayList<String>();// List to store
																// descriptor
																// definitions

			// Read the descriptor list file, and populate the arrays
			// if (file.exists()) {
			String line;
			StringTokenizer token;
			BufferedReader in = new BufferedReader(isr);
			line = in.readLine(); // contains the heading
			line = in.readLine(); // should contain the number, name and
									// description separated by commas
			String junk; // variable to store the tokens temporarily
			int number;

			while (in.ready() && line != null) {
				token = new StringTokenizer(line, "\t,;");
				number = (new Integer(token.nextToken())).intValue();
				num.add(new Integer(number)); // Adds descriptor number to the
												// array
				junk = token.nextToken();

				dname.add(new String(junk)); // Adds descriptor name to the
												// array
				junk = token.nextToken();


				desc.add(new String(junk)); //adds description to array
				line = in.readLine(); // Read the next line
			} // while in.ready()
			in.close();
			isr.close();
			ins.close();
			// }//if file.exists()


			for (int i = 0; i < dname.size(); i++) {
				this.htVarDef.put((String) dname.get(i), (String) desc.get(i));
			}
		} catch (Exception e) {
			e.printStackTrace();

		}
	}


	void loadDefs(String filename) {

		try {

			InputStream ins=this.getClass().getClassLoader().getResourceAsStream(filename);

			InputStreamReader isr=new InputStreamReader(ins);
			BufferedReader br=new BufferedReader(isr);

			String Line="12345";
			while (true) {
				Line=br.readLine();
				if (!(Line instanceof String)) break;

				LinkedList<String> ll=Utilities.Parse(Line,"\t");

				String strvar=(String)ll.get(0);
				String strdef=(String)ll.get(1);
				htVarDef.put(strvar,strdef);
//				System.out.println(strvar+"\t"+strdef);


			}
			br.close();
			isr.close();
			ins.close();
			//read the file that contains the list of 3D descriptors


		} catch (Exception e) {
			e.printStackTrace();
		}


	}

	public DescriptorData() {
		this.htVarDef=new Hashtable();

		if (Load2dDescriptions) this.loadDefs("variable definitions.txt");
		if (this.LoadMy3dDescriptions) this.loadDefs("variable definitions-My3d.txt");
		if (this.LoadRaghu3dDescriptions) this.load3DDefs();

		//generate Estate count array:
		strES_acnt=new String [strES.length];
		for (int i=0;i<strES.length;i++) strES_acnt[i]=strES[i]+"_acnt";

	}
	
	@Override
	public String toString() {
		Vector<String>fieldNames=getDescriptorNames();
		Vector<String>fieldValues=getDescriptorValues();
		
		if (fieldNames.size()!=fieldValues.size()) return "Error: # of fieldNames doesnt match # of fieldValues";
		
		String result=toStringDescriptorNames()+"\n"+toStringDescriptorValues();
			
		return result;
	}
	
	
	public String toStringDescriptorValues() {
		
		Vector<String>fieldValues=getDescriptorValues();
		String result="";
				
		for (int i=0;i<fieldValues.size();i++) {
			result+=fieldValues.get(i);
			if (i<fieldValues.size()-1) result+="\t";			
		}
		
		
		return result;
	}
	
	public static String toStringDescriptorNames() {
		
		Vector<String>fieldNames=getDescriptorNames();
		String result="";
				
		for (int i=0;i<fieldNames.size();i++) {
			result+=fieldNames.get(i);
			if (i<fieldNames.size()-1) result+="\t";			
		}
		
		
		return result;
	}

	
	public String toStringTranspose() {
		Vector<String>fieldNames=getDescriptorNames();
		Vector<String>fieldValues=getDescriptorValues();
		
		if (fieldValues.size()==1) return fieldValues.get(0);
		
		if (fieldNames.size()!=fieldValues.size()) return "Error: # of fieldNames doesnt match # of fieldValues";
		
		String result="";
		
		for (int i=0;i<fieldNames.size();i++) {
			result+=fieldNames.get(i)+"\t"+fieldValues.get(i)+"\r\n";
		}
				
		return result;
	}
	

	public static Vector<String> getDescriptorNames() {

		Vector<String>vecNames=new Vector<>();
						
		DescriptorData dd = new DescriptorData();

		String[] varlist = dd.varlist2d;

		try {
			for (int i = 0; i <= varlist.length - 1; i++) {
				Field myField = dd.getClass().getField(varlist[i]);
				String[] names = (String[]) myField.get(dd);

				for (int j = 0; j < names.length; j++) {
					vecNames.add(names[j]);
				}
			}

			for (int i = 0; i < strFragments.length; i++) {
				vecNames.add(strFragments[i]);
			}

		} catch (Exception ex) {
			System.out.println(ex.getMessage());
		}

//		System.out.println(vecNames.size());

		return vecNames;
	}
	
} // end class




