package gov.epa.TEST.Descriptors.DescriptorCalculations;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorUtilities.AtomicProperties;

//import ToxPredictor.Utilities.CDKUtilities;

import java.lang.reflect.*;
import java.util.List;
import java.util.Vector;
import java.io.*;

//TODO: what should sub return if have missing fragment?
// Just report sum for other fragments? 
//Or report as -9999 and then dont use descriptor if have this 
// value for any chemicals in cluster?

// this class calculates ALOGP (Ghose-Crippen LogKow) and the Ghose_Crippen molar refractivity


//TODO- method order changes when you use Method[] methods = this.getClass().getDeclaredMethods();- need to write code so methods are called in specific order!
//TODO for 606-34-8, ALOGP is sometimes	1.8614 and sometimes 2.6422- need to set order that methods are called since getDeclaredMethods is random!!!

/**
 * @author TMARTI02
 */
public class ALOGP {

	IAtomContainer m;

	DescriptorData dd;

	IRingSet rs;

	String[] Fragment; // estate fragments for each atom

	AtomicProperties ap; // needed to retrieve electronegativities

	public int[] frags = new int[121]; // counts of each type of fragment in
										// the molecule

	public int[] alogpfrag; // alogp fragments for each atom (used to see which
							// atoms have missing fragments)

	double[] fragval = new double[121];// coefficients for alogp model

	double[] refracval = new double[121]; // coefficients for refractivity
											// model

	String[] fragdesc = new String[121];

	public String UnassignedAtoms = "";

	public boolean debug = false;
//	 boolean debug=false;
	 public double ALogPwomissing;

	public ALOGP() {
		// fragments for ALOGP from Ghose et al., 1998
		fragval[1] = -1.5603;
		fragval[2] = -1.012;
		fragval[3] = -0.6681;
		fragval[4] = -0.3698;
		fragval[5] = -1.788;
		fragval[6] = -1.2486;
		fragval[7] = -1.0305;
		fragval[8] = -0.6805;
		fragval[9] = -0.3858;
		fragval[10] = 0.7555;
		fragval[11] = -0.2849;
		fragval[12] = 0.02;
		fragval[13] = 0.7894;
		fragval[14] = 1.6422;
		fragval[15] = -0.7866;
		fragval[16] = -0.3962;
		fragval[17] = 0.0383;
		fragval[18] = -0.8051;
		fragval[19] = -0.2129;
		fragval[20] = 0.2432;
		fragval[21] = 0.4697;
		fragval[22] = 0.2952;
		fragval[23] = 0;
		fragval[24] = -0.3251;
		fragval[25] = 0.1492;
		fragval[26] = 0.1539;
		fragval[27] = 0.0005;
		fragval[28] = 0.2361;
		fragval[29] = 0.3514;
		fragval[30] = 0.1814;
		fragval[31] = 0.0901;
		fragval[32] = 0.5142;
		fragval[33] = -0.3723;
		fragval[34] = 0.2813;
		fragval[35] = 0.1191;
		fragval[36] = -0.132;
		fragval[37] = -0.0244;
		fragval[38] = -0.2405;
		fragval[39] = -0.0909;
		fragval[40] = -0.1002;
		fragval[41] = 0.4182;
		fragval[42] = -0.2147;
		fragval[43] = -0.0009;
		fragval[44] = 0.1388;
		fragval[45] = 0;
		fragval[46] = 0.7341;
		fragval[47] = 0.6301;
		fragval[48] = 0.518;
		fragval[49] = -0.0371;
		fragval[50] = -0.1036;
		fragval[51] = 0.5234;
		fragval[52] = 0.6666;
		fragval[53] = 0.5372;
		fragval[54] = 0.6338;
		fragval[55] = 0.362;
		fragval[56] = -0.3567;
		fragval[57] = -0.0127;
		fragval[58] = -0.0233;
		fragval[59] = -0.1541;
		fragval[60] = 0.0324;
		fragval[61] = 1.052;
		fragval[62] = -0.7941;
		fragval[63] = 0.4165;
		fragval[64] = 0.6601;
		fragval[65] = 0;
		fragval[66] = -0.5427;
		fragval[67] = -0.3168;
		fragval[68] = 0.0132;
		fragval[69] = -0.3883;
		fragval[70] = -0.0389;
		fragval[71] = 0.1087;
		fragval[72] = -0.5113;
		fragval[73] = 0.1259;
		fragval[74] = 0.1349;
		fragval[75] = -0.1624;
		fragval[76] = -2.0585;
		fragval[77] = -1.915;
		fragval[78] = 0.4208;
		fragval[79] = -1.4439;
		fragval[80] = 0;
		fragval[81] = 0.4797;
		fragval[82] = 0.2358;
		fragval[83] = 0.1029;
		fragval[84] = 0.3566;
		fragval[85] = 0.1988;
		fragval[86] = 0.7443;
		fragval[87] = 0.5337;
		fragval[88] = 0.2996;
		fragval[89] = 0.8155;
		fragval[90] = 0.4856;
		fragval[91] = 0.8888;
		fragval[92] = 0.7452;
		fragval[93] = 0.5034;
		fragval[94] = 0.8995;
		fragval[95] = 0.5946;
		fragval[96] = 1.4201;
		fragval[97] = 1.1472;
		fragval[98] = 0;
		fragval[99] = 0.7293;
		fragval[100] = 0.7173;
		fragval[101] = 0;
		fragval[102] = -2.6737;
		fragval[103] = -2.4178;
		fragval[104] = -3.1121;
		fragval[105] = 0;
		fragval[106] = 0.6146;
		fragval[107] = 0.5906;
		fragval[108] = 0.8758;
		fragval[109] = -0.4979;
		fragval[110] = -0.3786;
		fragval[111] = 1.5188;
		fragval[112] = 1.0255;
		fragval[113] = 0;
		fragval[114] = 0;
		fragval[115] = 0;
		fragval[116] = -0.9359;
		fragval[117] = -0.1726;
		fragval[118] = -0.7966;
		fragval[119] = 0.6705;
		fragval[120] = -0.4801;

		fragdesc[1] = "CH3R / CH4";
		fragdesc[2] = "CH2R2";
		fragdesc[3] = "CHR3";
		fragdesc[4] = "CR4";
		fragdesc[5] = "CH3X";
		fragdesc[6] = "CH2RX";
		fragdesc[7] = "CH2X2";
		fragdesc[8] = "CHR2X";
		fragdesc[9] = "CHRX2";
		fragdesc[10] = "CHX3";
		fragdesc[11] = "CR3X";
		fragdesc[12] = "CR2X2";
		fragdesc[13] = "CRX3";
		fragdesc[14] = "CX4";
		fragdesc[15] = " =CH2";
		fragdesc[16] = " =CHR";
		fragdesc[17] = " =CR2";
		fragdesc[18] = " =CHX";
		fragdesc[19] = " =CRX";
		fragdesc[20] = " =CX2";
		fragdesc[21] = "#CH";
		fragdesc[22] = "#CR / R=C=R";
		fragdesc[23] = "#CX";
		fragdesc[24] = "R--CH--R";
		fragdesc[25] = "R--CR--R";
		fragdesc[26] = "R--CX--R";
		fragdesc[27] = "R--CH--X";
		fragdesc[28] = "R--CR--X";
		fragdesc[29] = "R--CX--X";
		fragdesc[30] = "X--CH--X";
		fragdesc[31] = "X--CR--X";
		fragdesc[32] = "X--CX--X";
		fragdesc[33] = "R--CH..X";
		fragdesc[34] = "R--CR..X";
		fragdesc[35] = "R--CX..X";
		fragdesc[36] = "Al-CH=X";
		fragdesc[37] = "Ar-CH=X";
		fragdesc[38] = "Al-C(=X)-Al";
		fragdesc[39] = "Ar-C(=X)-R";
		fragdesc[40] = "R-C(=X)-X / R-C#X / X=C=X";
		fragdesc[41] = "X-C(=X)-X";
		fragdesc[42] = "X--CH..X";
		fragdesc[43] = "X--CR..X";
		fragdesc[44] = "X--CX..X";
		fragdesc[45] = "undefined";
		fragdesc[46] = "H attached to C0(sp3) no X attached to next C";
		fragdesc[47] = "H attached to C1(sp3) / C0(sp2)";
		fragdesc[48] = "H attached to C2(sp3) / C1(sp2) / C0(sp)";
		fragdesc[49] = "H attached to C3(sp3) / C2(sp2) / C3(sp2) / C3(sp)";
		fragdesc[50] = "H attached to heteroatom";
		fragdesc[51] = "H attached to alpha-C";
		fragdesc[52] = "H attached to C0(sp3) with 1X attached to next C";
		fragdesc[53] = "H attached to C0(sp3) with 2X attached to next C";
		fragdesc[54] = "H attached to C0(sp3) with 3X attached to next C";
		fragdesc[55] = "H attached to C0(sp3) with 4X attached to next C";
		fragdesc[56] = "alcohol";
		fragdesc[57] = "phenol / enol / carboxyl OH";
		fragdesc[58] = "=O";
		fragdesc[59] = "Al-O-Al";
		fragdesc[60] = "Al-O-Ar / Ar-O-Ar / R..O..R / R-O-C=X";
		fragdesc[61] = "O--";
		fragdesc[62] = "O- (negatively charged)";
		fragdesc[63] = "R-O-O-R";
		fragdesc[64] = "Any-Se-Any";
		fragdesc[65] = "=Se";
		fragdesc[66] = "Al-NH2";
		fragdesc[67] = "Al2-NH";
		fragdesc[68] = "Al3-N";
		fragdesc[69] = "Ar-NH2 / X-NH2";
		fragdesc[70] = "Ar-NH-Al";
		fragdesc[71] = "Ar-NAl2";
		fragdesc[72] = "RCO-N< / >N-X=X";
		fragdesc[73] = "Ar2NH / Ar3N / Ar2N-Al / R..N..R";
		fragdesc[74] = "R#N / R=N-";
		fragdesc[75] = "R--N--R / R--N--X";
		fragdesc[76] = "Ar-NO2 / R--N(--R)--O / RO-NO";
		fragdesc[77] = "Al-NO2";
		fragdesc[78] = "Ar-N=X / X-N=X";
		fragdesc[79] = "N+ (positively charged)";
		fragdesc[80] = "undefined";
		fragdesc[81] = "F attached to C1(sp3)";
		fragdesc[82] = "F attached to C2(sp3)";
		fragdesc[83] = "F attached to C3(sp3)";
		fragdesc[84] = "F attached to C1(sp2)";
		fragdesc[85] = "F attached to C2(sp2)-C4(sp2) / C1(sp) / C4(sp3) / X";
		fragdesc[86] = "Cl attached to C1(sp3)";
		fragdesc[87] = "Cl attached to C2(sp3)";
		fragdesc[88] = "Cl attached to C3(sp3)";
		fragdesc[89] = "Cl attached to C1(sp2)";
		fragdesc[90] = "Cl attached to C2(sp2)-C4(sp2) / C1(sp) / C4(sp3) / X";
		fragdesc[91] = "Br attached to C1(sp3)";
		fragdesc[92] = "Br attached to C2(sp3)";
		fragdesc[93] = "Br attached to C3(sp3)";
		fragdesc[94] = "Br attached to C1(sp2)";
		fragdesc[95] = "Br attached to C2(sp2)-C4(sp2) / C1(sp) / C4(sp3) / X";
		fragdesc[96] = "I attached to C1(sp3)";
		fragdesc[97] = "I attached to C2(sp3)";
		fragdesc[98] = "I attached to C3(sp3)";
		fragdesc[99] = "I attached to C1(sp2)";
		fragdesc[100] = "I attached to C2(sp2)-C4(sp2) / C1(sp) / C4(sp3) / X";
		fragdesc[101] = "fluoride ion";
		fragdesc[102] = "chloride ion";
		fragdesc[103] = "bromide ion";
		fragdesc[104] = "iodide ion";
		fragdesc[105] = "undefined";
		fragdesc[106] = "R-SH";
		fragdesc[107] = "R2S / RS-SR";
		fragdesc[108] = "R=S";
		fragdesc[109] = "R-SO-R";
		fragdesc[110] = "R-SO2-R";
		fragdesc[111] = ">Si<";
		fragdesc[112] = ">B- as in boranes";
		fragdesc[113] = "undefined";
		fragdesc[114] = "undefined";
		fragdesc[115] = "P ylids";
		fragdesc[116] = "R3-P=X";
		fragdesc[117] = "X3-P=X (phosphate)";
		fragdesc[118] = "PX3 (phosphite)";
		fragdesc[119] = "PR3 (phosphine)";
		fragdesc[120] = "C-P(X)2=X (phosphonate)";

		// fragments for AMR from Viswanadhan et al., 1989
		refracval[1] = 2.968;
		refracval[2] = 2.9116;
		refracval[3] = 2.8028;
		refracval[4] = 2.6205;
		refracval[5] = 3.015;
		refracval[6] = 2.9244;
		refracval[7] = 2.6329;
		refracval[8] = 2.504;
		refracval[9] = 2.377;
		refracval[10] = 2.5559;
		refracval[11] = 2.303;
		refracval[12] = 2.3006;
		refracval[13] = 2.9627;
		refracval[14] = 2.3038;
		refracval[15] = 3.2001;
		refracval[16] = 4.2654;
		refracval[17] = 3.9392;
		refracval[18] = 3.6005;
		refracval[19] = 4.487;
		refracval[20] = 3.2001;
		refracval[21] = 3.4825;
		refracval[22] = 4.2817;
		refracval[23] = 3.9556;
		refracval[24] = 3.4491;
		refracval[25] = 3.8821;
		refracval[26] = 3.7593;
		refracval[27] = 2.5009;
		refracval[28] = 2.5;
		refracval[29] = 3.0627;
		refracval[30] = 2.5009;
		refracval[31] = 0;
		refracval[32] = 2.6632;
		refracval[33] = 3.4671;
		refracval[34] = 3.6842;
		refracval[35] = 2.9372;
		refracval[36] = 4.019;
		refracval[37] = 4.777;
		refracval[38] = 3.9031;
		refracval[39] = 3.9964;
		refracval[40] = 3.4986;
		refracval[41] = 3.4997;
		refracval[42] = 2.7784;
		refracval[43] = 2.6267;
		refracval[44] = 2.5;
		refracval[45] = 0;
		refracval[46] = 0.8447;
		refracval[47] = 0.8939;
		refracval[48] = 0.8005;
		refracval[49] = 0.832;
		refracval[50] = 0.8;
		refracval[51] = 0.8188;
		refracval[52] = 0.9215;
		refracval[53] = 0.9769;
		refracval[54] = 0.7701;
		refracval[55] = 0;
		refracval[56] = 1.7646;
		refracval[57] = 1.4778;
		refracval[58] = 1.4429;
		refracval[59] = 1.6191;
		refracval[60] = 1.3502;
		refracval[61] = 1.945;
		refracval[62] = 0;
		refracval[63] = 0;
		refracval[64] = 11.1366;
		refracval[65] = 13.1149;
		refracval[66] = 2.6221;
		refracval[67] = 2.5;
		refracval[68] = 2.898;
		refracval[69] = 3.6841;
		refracval[70] = 4.2808;
		refracval[71] = 3.6189;
		refracval[72] = 2.5;
		refracval[73] = 2.7956;
		refracval[74] = 2.7;
		refracval[75] = 4.2063;
		refracval[76] = 4.0184;
		refracval[77] = 3.0009;
		refracval[78] = 4.7142;
		refracval[79] = 0;
		refracval[80] = 0;
		refracval[81] = 0.8725;
		refracval[82] = 1.1837;
		refracval[83] = 1.1573;
		refracval[84] = 0.8001;
		refracval[85] = 1.5013;
		refracval[86] = 5.6156;
		refracval[87] = 6.1022;
		refracval[88] = 5.9921;
		refracval[89] = 5.3885;
		refracval[90] = 6.1363;
		refracval[91] = 8.5991;
		refracval[92] = 8.9188;
		refracval[93] = 8.8006;
		refracval[94] = 8.2065;
		refracval[95] = 8.7352;
		refracval[96] = 13.9462;
		refracval[97] = 14.0792;
		refracval[98] = 14.073;
		refracval[99] = 12.9918;
		refracval[100] = 13.3408;
		refracval[101] = 0;
		refracval[102] = 0;
		refracval[103] = 0;
		refracval[104] = 0;
		refracval[105] = 0;
		refracval[106] = 7.8916;
		refracval[107] = 7.7935;
		refracval[108] = 9.4338;
		refracval[109] = 7.7223;
		refracval[110] = 5.7558;
		refracval[111] = 0;
		refracval[112] = 0;
		refracval[113] = 0;
		refracval[114] = 0;
		refracval[115] = 0;
		refracval[116] = 5.5306;
		refracval[117] = 5.5152;
		refracval[118] = 6.836;
		refracval[119] = 10.0101;
		refracval[120] = 5.2806;

		try {
			ap = AtomicProperties.getInstance();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private void FindUnassignedAtoms() {
		UnassignedAtoms = "";

		for (int i = 0; i <= m.getAtomCount() - 1; i++) {
			if (alogpfrag[i] == 0)
				UnassignedAtoms += Fragment[i] +", ";
//			UnassignedAtoms += Fragment[i] +" "+(i + 1) + ",";
		}
		// System.out.println(UnassignedAtoms);

	}

	
//	private Vector<String> getMethodOrderArray() {
//		
//		Vector <String> order=new Vector<String>();
//		//TODO: optimize order to match dragon!
//
//		//ring atoms:
//		order.add("CalcGroup025_026_028_029_031_032_034_035_043_044");
//		order.add("CalcGroup024_027_030_033_042");
//
//		//Carbon:
//		order.add("CalcGroup021_to_023_040");//#C
//		order.add("CalcGroup015");//=CH2
//		order.add("CalcGroup016_018_036_037_040");//=CH
//		order.add("CalcGroup017_019_020_038_to_041");//=C
//
//		order.add("CalcGroup001_005"); //CH4, CH3
//		order.add("CalcGroup002_006_007");//CH2
//		order.add("CalcGroup003_008_009_010");//CH
//		order.add("CalcGroup004_011_to_014");//C
//
//		//Halogens:
//		order.add("CalcGroup101_to_104");
//		order.add("CalcGroup081_to_085");
//		order.add("CalcGroup086_to_090");
//		order.add("CalcGroup091_to_095");
//		order.add("CalcGroup096_to_100");
//
//		//Hydrogens:
//		order.add("CalcGroup046_to_055");
//
//		//Sulfur:
//		order.add("CalcGroup106");
//		order.add("CalcGroup107");
//		order.add("CalcGroup108");
//		order.add("CalcGroup109");
//		order.add("CalcGroup110");
//
//		//Phosphorus:
//		order.add("CalcGroup116_117_120");
//		order.add("CalcGroup118_119");
//
//		//Nitrogen:
//		order.add("CalcGroup066_to_079");
//
//		//Oxygen:
//		order.add("CalcGroup056_57");
//		order.add("CalcGroup058_61");
//		order.add("CalcGroup059_060_063");
//
//		//Misc elements:
//		order.add("CalcGroup064");//Se
//		order.add("CalcGroup111");//Si
//		order.add("CalcGroup112");//B
//
//		return order;
//		
//	}

	//12/11/15 now manually assign the fragments in a specified order (previously the order was somewhat random)
	private void assignFrags(int i) {

		//Hydrogens:
		CalcGroup046_to_055(i);

		if (alogpfrag[i] != 0) return;

		//****************************************************************************************
		//ring atoms:
		CalcGroup025_026_028_029_031_032_034_035_043_044(i); 
		if (alogpfrag[i] != 0) return;
		
		CalcGroup024_027_030_033_042(i);
		if (alogpfrag[i] != 0) return;
		
		//****************************************************************************************
		//Carbon:
		CalcGroup021_to_023_040(i);//#C
		if (alogpfrag[i] != 0) return;
		
		CalcGroup015(i);//=CH2
		if (alogpfrag[i] != 0) return;
		
		CalcGroup016_018_036_037_040(i);//=CH
		if (alogpfrag[i] != 0) return;
		
		CalcGroup017_019_020_038_to_041(i);//=C
		if (alogpfrag[i] != 0) return;

		CalcGroup001_005(i); //CH4, CH3
		if (alogpfrag[i] != 0) return;

		CalcGroup002_006_007(i);//CH2
		if (alogpfrag[i] != 0) return;

		CalcGroup003_008_009_010(i);//CH
		if (alogpfrag[i] != 0) return;

		CalcGroup004_011_to_014(i);//C
		if (alogpfrag[i] != 0) return;

		//****************************************************************************************
		//Halogens:
		CalcGroup101_to_104(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup081_to_085(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup086_to_090(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup091_to_095(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup096_to_100(i);
		if (alogpfrag[i] != 0) return;

		//****************************************************************************************
		//Sulfur:
		CalcGroup106(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup107(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup108(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup109(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup110(i);
		if (alogpfrag[i] != 0) return;
		
		//****************************************************************************************
		//Phosphorus:
		CalcGroup116_117_120(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup118_119(i);
		if (alogpfrag[i] != 0) return;

		//****************************************************************************************
		//Nitrogen:
		CalcGroup066_to_079(i);
		if (alogpfrag[i] != 0) return;

		//****************************************************************************************
		//Oxygen:
		CalcGroup056_57(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup058_61(i);
		if (alogpfrag[i] != 0) return;

		CalcGroup059_060_063(i);
		if (alogpfrag[i] != 0) return;

		//****************************************************************************************
		//Misc elements:
		CalcGroup064(i);//Se
		if (alogpfrag[i] != 0) return;
		
		CalcGroup111(i);//Si
		if (alogpfrag[i] != 0) return;
		
		CalcGroup112(i);//B
		if (alogpfrag[i] != 0) return;
		

	}
	
	//this method needs molecules w/implicit hydrogens:
	
	public void Calculate(IAtomContainer m, DescriptorData dd, String[] Fragment,
			IRingSet rs) {
		this.m = m;

//		Vector<String>order=this.getMethodOrderArray();		
		
		this.dd = dd;
		this.Fragment = Fragment;
		this.rs = rs;

		alogpfrag = new int[m.getAtomCount()];

		
//		Method[] methods = this.getClass().getDeclaredMethods();
//		
//		for (int i=0;i<meth.length;i++) {
//			if (meth[i].getName().indexOf("CalcGroup")>-1) System.out.println(meth[i].getName());
//		}
		
		
		for (int i = 1; i <= 120; i++) {
			frags[i] = 0;
		}

		for (int i = 0; i <= m.getAtomCount() - 1; i++) {

			alogpfrag[i] = 0;
			try {

				if (Fragment[i] instanceof String) {
					this.assignFrags(i);
					
//					for (int j = 0; j < order.size(); j++) {	
//						
//						Method method = getMethod(order.get(j));
//
//						// skip method of atom is already assigned to a
//						// fragment:
//						if (alogpfrag[i] != 0
//								&& method.getName().indexOf(
//										"CalcGroup046_to_055") == -1) {							
//							continue;
//						}
//						
//						System.out.println(i+"\t"+method.getName());
//						Object[] objs = { new Integer(i) };
//						method.invoke(this, objs);
//						
//					}
					
//					System.out.println(i+"\t"+alogpfrag[i]);
				}

			} catch (Exception e) {
				e.printStackTrace();
			}
		} // end i atom loop

		dd.ALOGP = 0;
		dd.AMR = 0;

		for (int i = 1; i <= 120; i++) {
			dd.ALOGP += fragval[i] * frags[i];
			dd.AMR += refracval[i] * frags[i];

			if (debug) {
				if (frags[i] > 0)
					System.out.println(i + "\t" + frags[i] + "\t" + fragval[i]);
			}
		}
//		System.out.println("ALOGP\t"+dd.ALOGP);

		// System.out.println("");


		if (debug) {
			for (int i = 0; i <= m.getAtomCount() - 1; i++) {
				System.out.println((i) + "\t" + Fragment[i] + "\t"
						+ alogpfrag[i]+"\t"+fragdesc[alogpfrag[i]]);
			}
		}

		 

		this.FindUnassignedAtoms();

		// for now ignore the fact we may have missing frags:
//		if (!this.UnassignedAtoms.equals("")) {
//			this.ALogPwomissing=dd.ALOGP;			
//			dd.ALOGP = Double.NaN; // if have unassigned atoms dont report a val
//			dd.ALOGP2 = Double.NaN;
//			dd.AMR=Double.NaN;
//		} else {
			dd.ALOGP2 = dd.ALOGP * dd.ALOGP;
//		}

		
		
		if (debug)
			System.out.println("ALOGP=" + dd.ALOGP);
		if (debug)
			System.out.println("Hydrogen count=" + dd.nH);

		int HCountCalc = 0;

		for (int i = 46; i <= 55; i++) {
			HCountCalc += frags[i];
		}

		if (debug)
			System.out.println("Hydrogen count from frags=" + HCountCalc);

	}

//	private Method getMethod (String methodName) {
//		
//		Method[] methods = this.getClass().getDeclaredMethods();
//
//		for (int k=0;k<methods.length;k++) {
//			if (methods[k].getName().equals(methodName)) {
//				return methods[k];
//			}
//		}
//
//		return null;
//	}
	
	private void CalcGroup001_005(int i) {
		// C in CH3R

		//12/5/08 TMM: added to account for CH4:
		if (Fragment[i].equals("SCH4") || Fragment[i].equals("SCH4*")) {
			frags[1]++;
			alogpfrag[i] = 1;
			frags[46]+=4;
		}
		
		if (Fragment[i].equals("SsCH3")) {
			List ca = m.getConnectedAtomsList(m.getAtom(i));
			
			if (((IAtom)ca.get(0)).getSymbol().equals("C")) {
				frags[1]++;
				alogpfrag[i] = 1;
			} else {
				frags[5]++;
				alogpfrag[i] = 5;
			}

		}

	}

	private void CalcGroup002_006_007(int i) {
		// C in CH2RX

		if (Fragment[i].equals("SssCH2")) {

			List ca = m.getConnectedAtomsList(m.getAtom(i));

			int CarbonCount = 0;
			int HeteroCount = 0;

			// System.out.println("here");
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else
					HeteroCount++;
			}

			if (CarbonCount == 2 && HeteroCount == 0) {
				frags[2]++;
				alogpfrag[i] = 2;
			} else if (CarbonCount == 1 && HeteroCount == 1) {
				frags[6]++;
				alogpfrag[i] = 6;
			} else if (CarbonCount == 0 && HeteroCount == 2) {
				frags[7]++;
				alogpfrag[i] = 7;
			}

		}

	}

	private void CalcGroup003_008_009_010(int i) {

		if (Fragment[i].equals("SsssCH")) {

			List ca = m.getConnectedAtomsList(m.getAtom(i));

			int CarbonCount = 0;
			int HeteroCount = 0;

			// System.out.println("here");
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else
					HeteroCount++;
			}

			if (CarbonCount == 3 && HeteroCount == 0) {
				frags[3]++;
				alogpfrag[i] = 3;
			} else if (CarbonCount == 2 && HeteroCount == 1) {
				frags[8]++;
				alogpfrag[i] = 8;
			} else if (CarbonCount == 1 && HeteroCount == 2) {
				frags[9]++;
				alogpfrag[i] = 9;
			} else if (CarbonCount == 0 && HeteroCount == 3) {
				frags[10]++;
				alogpfrag[i] = 10;
			}

		}

	}

	private void CalcGroup004_011_to_014(int i) {
		// C in CH2RX

		if (Fragment[i].equals("SssssC")) {

			List ca = m.getConnectedAtomsList(m.getAtom(i));

			int CarbonCount = 0;
			int HeteroCount = 0;

			// System.out.println("here");
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else
					HeteroCount++;
			}

			if (CarbonCount == 4 && HeteroCount == 0) {
				frags[4]++;
				alogpfrag[i] = 4;
			} else if (CarbonCount == 3 && HeteroCount == 1) {
				frags[11]++;
				alogpfrag[i] = 11;
			} else if (CarbonCount == 2 && HeteroCount == 2) {
				frags[12]++;
				alogpfrag[i] = 12;
			} else if (CarbonCount == 1 && HeteroCount == 3) {
				frags[13]++;
				alogpfrag[i] = 13;
			} else if (CarbonCount == 0 && HeteroCount == 4) {
				frags[14]++;
				alogpfrag[i] = 14;
			}

		}

	}

	private void CalcGroup015(int i) {
		if (!Fragment[i].equals("SdCH2"))
			return;

		
		List ca = m.getConnectedAtomsList(m.getAtom(i));

		if (!((IAtom)ca.get(0)).getSymbol().equals("C"))
			return; // dont assign 15 to formaldehyde

		frags[15]++;
		alogpfrag[i] = 15;

		

	}

	private void CalcGroup016_018_036_037_040(int i) {

		IAtom ai = m.getAtom(i);

		if (!Fragment[i].equals("SdsCH") && !Fragment[i].equals("SdCH2"))
			return;

		List ca = m.getConnectedAtomsList(ai);

		boolean HaveCdX = false;
		boolean HaveCsX = false;
		boolean HaveCsAr = false;

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)) {
				if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
					HaveCsX = true;
				}

				if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)) {
					HaveCsAr = true;
				}

			} else if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.DOUBLE)) {
				if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
					HaveCdX = true;
				}

			}
		}
		

		if (HaveCdX) {
			
			//TMM 12/5/08- redid logic so that R-C(=X)X takes precedence over Ar-CH=X:
			if (HaveCsX) {
				frags[40]++;
				alogpfrag[i] = 40;
			} else {
				if (HaveCsAr) {
					frags[37]++;
					alogpfrag[i] = 37;
				} else {
					frags[36]++;
					alogpfrag[i] = 36;
				}
			}

//			if (HaveCsAr) {
//				frags[37]++;
//				alogpfrag[i] = 37;
//			} else {
//				if (HaveCsX) {//TMM 11/13/08: added possibility for group 40 here
//					frags[40]++;
//					alogpfrag[i] = 40;
//				} else {
//					frags[36]++;
//					alogpfrag[i] = 36;
//				}
//			}
		} else {
			if (HaveCsX) {
				frags[18]++;
				alogpfrag[i] = 18;
			} else {
				frags[16]++;
				alogpfrag[i] = 16;
			}
		}
	}

	private void CalcGroup017_019_020_038_to_041(int i) {

		IAtom ai = m.getAtom(i);

		if (!Fragment[i].equals("SdssC"))
			return;

		List ca = m.getConnectedAtomsList(ai);

		int RCount = 0;
		int XCount = 0;
		int NonAromaticXCount = 0;

		boolean HaveCdX = false;

		int AliphaticCount = 0;
		int AromaticCount = 0;

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
					RCount++;
				} else {
					XCount++;

					if (!((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)) {
						NonAromaticXCount++;
					}
				}

				if (!((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)) {
					AliphaticCount++;
				} else {
					AromaticCount++;
				}

			} else if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.DOUBLE)) {
				if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
					HaveCdX = true;
				}
			}
		}

		if (HaveCdX) {
			if (AromaticCount >= 1 && NonAromaticXCount == 0) { // Ar-C(=X)-R
				// TODO: add code to check if have R or X for nonaromatic
				// attachment to C?
				// if we do this check we would have missing fragment for
				// Ar-C(=X)-X
				// TODO: which fragment to use if we have Ar-C(=X)-Ar? Currently
				// this frag is used

				frags[39]++;
				alogpfrag[i] = 39;
			} else {

				if (RCount == 1 && XCount == 1) {
					frags[40]++;
					alogpfrag[i] = 40;
				} else if (RCount == 0 && XCount == 2) {
					frags[41]++;
					alogpfrag[i] = 41;
				} else {
					if (AromaticCount == 0) {
						frags[38]++;
						alogpfrag[i] = 38;
						return;
					}
				}

			}

		} else {
			if (RCount == 2 && XCount == 0) {
				frags[17]++;
				alogpfrag[i] = 17;
			} else if (RCount == 1 && XCount == 1) {
				frags[19]++;
				alogpfrag[i] = 19;
			} else if (RCount == 0 && XCount == 2) {
				frags[20]++;
				alogpfrag[i] = 20;
			}
		}

	}

//	private void CalcGroup021_to_023_040(int i) {
//
//		List ca = m.getConnectedAtomsList(m.getAtom(i));
//		IAtom ai = m.getAtom(i);
//
//		if (Fragment[i].equals("StCH")) {
//			frags[21]++;
//			alogpfrag[i] = 21;
//		} else if (Fragment[i].equals("SddC")) {
//			if (((IAtom)ca.get(0)).getSymbol().equals("C") && ((IAtom)ca.get(1)).getSymbol().equals("C")) {// R==C==R
//				frags[22]++;
//				alogpfrag[i] = 22;
//			} else if (!((IAtom)ca.get(0)).getSymbol().equals("C")
//					&& !((IAtom)ca.get(1)).getSymbol().equals("C")) {// X==C==X
//				frags[40]++;
//				alogpfrag[i] = 40;
//			}
//		} else if (Fragment[i].equals("StsC")) {
//
//			boolean HaveCtX = false;
//			boolean HaveCsX = false;
//
//			for (int j = 0; j <= ca.size() - 1; j++) {
//				if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.Single)) {
//					if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
//						HaveCsX = true;
//					}
//				} else if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.Triple)) {
//					if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
//						HaveCtX = true;
//					}
//				}
//
//			}
//
//			//TMM 11/13/08: removed !HaveCsX  
////			if (HaveCtX && !HaveCsX) {
//			if (HaveCtX) {
//				frags[40]++;
//				alogpfrag[i] = 40;
//			} else {
//				if (HaveCsX) {// #C-X
//					frags[23]++;
//					alogpfrag[i] = 23;
//				} else {// #C-R
//					frags[22]++;
//					alogpfrag[i] = 22;
//				}
//			}
//			
//
//		}
//
//	}

	//TMM 11/20/08: method redone
	private void CalcGroup021_to_023_040(int i) {

		List ca = m.getConnectedAtomsList(m.getAtom(i));
		IAtom ai = m.getAtom(i);

		if (Fragment[i].equals("SddC")) {
			if (((IAtom)ca.get(0)).getSymbol().equals("C") && ((IAtom)ca.get(1)).getSymbol().equals("C")) {// R==C==R
				frags[22]++;
				alogpfrag[i] = 22;
			} else if (!((IAtom)ca.get(0)).getSymbol().equals("C")
					&& !((IAtom)ca.get(1)).getSymbol().equals("C")) {// X==C==X
				frags[40]++;
				alogpfrag[i] = 40;
			}
		} else if (Fragment[i].equals("StsC") || Fragment[i].equals("StCH") ) {

			boolean HaveCtX = false;
			boolean HaveCsX = false;

			for (int j = 0; j <= ca.size() - 1; j++) {
				if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)) {
					if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
						HaveCsX = true;
					}
				} else if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.TRIPLE)) {
					if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
						HaveCtX = true;
					}
				}

			}

			if (HaveCtX) {
				frags[40]++;
				alogpfrag[i] = 40;
			} else {
				if (HaveCsX) {// #C-X
					frags[23]++;
					alogpfrag[i] = 23;
				} else {// #C-R
					if (Fragment[i].equals("StCH")) {
						frags[21]++;
						alogpfrag[i] = 21;
					} else if (Fragment[i].equals("StsC")) {
						frags[22]++;
						alogpfrag[i] = 22;
					}
				}
			}
		}
	}
	
	private void CalcGroup024_027_030_033_042(int i) {
		// 24: C in R--CH--R
		// 27: C in R--CH--X
		// 30: C in X--CH--X
		// 33: C in R--CH...X
		// 42: C in X--CH...X

		if (!Fragment[i].equals("SaaCH"))
			return;
		// System.out.println("here");

		IAtom ai = m.getAtom(i);
		List ca = m.getConnectedAtomsList(m.getAtom(i));

		if (((IAtom)ca.get(0)).getSymbol().equals("C") && ((IAtom)ca.get(1)).getSymbol().equals("C")) {
			frags[24]++;
			alogpfrag[i] = 24;
			return;
		}

		// check if both hetero atoms have at least one double bond
		List bonds = m.getConnectedBondsList(((IAtom)ca.get(0)));

		boolean HaveDouble1 = false;

		for (int k = 0; k <= bonds.size() - 1; k++) {
			if (((IBond)bonds.get(k)).getOrder().equals(IBond.Order.DOUBLE)) {
				HaveDouble1 = true;
				break;
			}

		}

		bonds = m.getConnectedBondsList(((IAtom)ca.get(1)));

		boolean HaveDouble2 = false;

		for (int k = 0; k <= bonds.size() - 1; k++) {
			if (((IBond)bonds.get(k)).getOrder().equals(IBond.Order.DOUBLE)) {
				HaveDouble2 = true;
				break;
			}

		}

		if (!((IAtom)ca.get(0)).getSymbol().equals("C") && !((IAtom)ca.get(1)).getSymbol().equals("C")) {
			if (HaveDouble1 && HaveDouble2) { // X--CH--X
				frags[30]++;
				alogpfrag[i] = 30;
			} else { // X--CH...X
				frags[42]++;
				alogpfrag[i] = 42;
			}

		} else if ((((IAtom)ca.get(0)).getSymbol().equals("C") && !((IAtom)ca.get(1)).getSymbol().equals(
				"C"))
				|| (!((IAtom)ca.get(0)).getSymbol().equals("C") && ((IAtom)ca.get(1)).getSymbol().equals(
						"C"))) {

			if (HaveDouble1 && HaveDouble2) { // R--CH--X
				frags[27]++;
				alogpfrag[i] = 27;
			} else {// R--CH...X
				frags[33]++;
				alogpfrag[i] = 33;

			}

		}

	}

	private void CalcGroup025_026_028_029_031_032_034_035_043_044(int i) { 
		//this sub was modified on 8/15/07 to resolve problems for 
		// groups 26&34 for compounds like 132-64-9
		
		// 25: R--CR--R
		// 26: R--CX--R
		// 28: R--CR--X
		// 29: R--CX--X
		// 31: X--CR--X
		// 32: X--CX--X
		// 34: R--CR...X
		// 35: R--CX...X
		// 43: X--CR...X
		// 44: X--CX...X

		if (!Fragment[i].equals("SsaaC") && !Fragment[i].equals("SaaaC"))
			return;

		IAtom ai = m.getAtom(i);
		List ca = m.getConnectedAtomsList(m.getAtom(i));

		IAtom[] sameringatoms = new IAtom[2];
		IAtom nonringatom = null;

		int sameringatomscount = 0;
		for (int j = 0; j <= ca.size() - 1; j++) {

			IAtom caj=(IAtom) ca.get(j);

			if (EStateFragmentDescriptor.InSameAromaticRing(m, ai,caj, rs)) {
				sameringatomscount++;				
			}

		}

		if (sameringatomscount == 3) {

			int XCount = 0;

			for (int j = 0; j < ca.size(); j++) {
				IAtom caj = (IAtom) ca.get(j);
				if (!caj.getSymbol().equals("C")) {
					XCount++;
				}
			}

			if (XCount == 0) {
				frags[25]++;
				alogpfrag[i] = 25;
			} else if (XCount == 1) {

				IAtom xatom = null;
				for (int j = 0; j < ca.size(); j++) {
					IAtom caj = (IAtom) ca.get(j);
					if (!caj.getSymbol().equals("C")) {
						xatom = caj;
						break;
					}
				}

				List bonds = m.getConnectedBondsList(xatom);

				boolean HaveDouble = false;

				for (int k = 0; k <= bonds.size() - 1; k++) {
					if (((IBond) bonds.get(k)).getOrder().equals(IBond.Order.DOUBLE)) {
						HaveDouble = true;
						break;
					}
				}

				if (HaveDouble) {// R--CR--X
					frags[28]++;
					alogpfrag[i] = 28;
				} else {
					// R--CR...X
					frags[34]++;
					alogpfrag[i] = 34;
				}

			} else if (XCount==2){
				
				int HaveDoubleCount=0;
				
				for (int j = 0; j < ca.size(); j++) {
					IAtom caj = (IAtom) ca.get(j);
					if (!caj.getSymbol().equals("C")) {
						List bonds = m.getConnectedBondsList(caj);
						
						for (int k = 0; k <= bonds.size() - 1; k++) {
							if (((IBond) bonds.get(k)).getOrder().equals(IBond.Order.DOUBLE)) {
								HaveDoubleCount++;
								break;
							}
						}
						
					}
				}

				if (HaveDoubleCount==2) {
					//X--CR--X
					frags[31]++;
					alogpfrag[i] = 31;
				} else if (HaveDoubleCount==1) {
					//X--CR...X
					frags[43]++;
					alogpfrag[i] = 43;
					
				} else {
					System.out.println("XCount==2, HaveDoubleCount==0");
				}

				
			} else { //XCount==3
//				System.out.println("Have XCount>2 for sameringatomscount==3, need to implement");

				int HaveDoubleCount=0;
				
				for (int j = 0; j < ca.size(); j++) {
					IAtom caj = (IAtom) ca.get(j);
					if (!caj.getSymbol().equals("C")) {
						List bonds = m.getConnectedBondsList(caj);
						
						for (int k = 0; k <= bonds.size() - 1; k++) {
							if (((IBond) bonds.get(k)).getOrder().equals(IBond.Order.DOUBLE)) {
								HaveDoubleCount++;
								break;
							}
						}
						
					}
				}

				if (HaveDoubleCount>=2) {
					// X--CX--X
					frags[32]++;
					alogpfrag[i] = 32;

				} else if (HaveDoubleCount==1) {
					// X--CX...X
					frags[44]++;
					alogpfrag[i] = 44;

				} else {
					System.out.println("XCount==3, HaveDoubleCount==0");
				}


			}

		} else if (sameringatomscount == 2) {
			int count = 0;
			for (int j = 0; j <= ca.size() - 1; j++) {

				IAtom caj = (IAtom) ca.get(j);

				if (EStateFragmentDescriptor.InSameAromaticRing(m, ai,caj, rs)) {
					sameringatoms[count] = caj;
					count++;
				} else {
					nonringatom = caj;
				}
				
			}

			// check if both hetero atoms have at least one double bond
			List bonds = m.getConnectedBondsList(sameringatoms[0]);

			boolean HaveDouble1 = false;

			for (int k = 0; k <= bonds.size() - 1; k++) {
				if (((IBond) bonds.get(k)).getOrder().equals(IBond.Order.DOUBLE)) {
					HaveDouble1 = true;
					break;
				}

			}

			bonds = m.getConnectedBondsList(sameringatoms[1]);

			boolean HaveDouble2 = false;

			for (int k = 0; k <= bonds.size() - 1; k++) {
				if (((IBond) bonds.get(k)).getOrder().equals(IBond.Order.DOUBLE)) {
					HaveDouble2 = true;
					break;
				}

			}

			if (!sameringatoms[0].getSymbol().equals("C")
					&& !sameringatoms[1].getSymbol().equals("C")) {
				if (HaveDouble1 && HaveDouble2) { // X--CR--X
					if (nonringatom.getSymbol().equals("C")) {
						frags[31]++;
						alogpfrag[i] = 31;
					} else { // X--CX--X
						frags[32]++;
						alogpfrag[i] = 32;
					}

				} else {

					if (nonringatom.getSymbol().equals("C")) { // X--CR..X
						frags[43]++;
						alogpfrag[i] = 43;

					} else { // X--CX...X
						frags[44]++;
						alogpfrag[i] = 44;
					}

				}
			} else if (sameringatoms[0].getSymbol().equals("C")
					&& sameringatoms[1].getSymbol().equals("C")) {

				if (nonringatom.getSymbol().equals("C")) {// R--CR--R
					frags[25]++;
					alogpfrag[i] = 25;
				} else {
						// R--CX--R
						frags[26]++;
						alogpfrag[i] = 26;
				}

			} else if ((sameringatoms[0].getSymbol().equals("C") && !sameringatoms[1]
					.getSymbol().equals("C"))
					|| (!sameringatoms[0].getSymbol().equals("C") && sameringatoms[1]
							.getSymbol().equals("C"))) {

				if (HaveDouble1 && HaveDouble2) { // R--CR--X
					if (nonringatom.getSymbol().equals("C")) {
						frags[28]++;
						alogpfrag[i] = 28;
					} else { // R--CX--X
						frags[29]++;
						alogpfrag[i] = 29;
					}

				} else {

					if (nonringatom.getSymbol().equals("C")) { // R--CR...X
						frags[34]++;
						alogpfrag[i] = 34;

					} else { // R--CX...X
						frags[35]++;
						alogpfrag[i] = 35;
					}

				}

			}
		}

	}

	private void CalcGroup046_to_055(int i) {

		
		IAtom ai = m.getAtom(i);
		
		List ca = m.getConnectedAtomsList(m.getAtom(i));

		if (ai.getImplicitHydrogenCount() == 0)
			return;

		// !ai.getFlag(CDKConstants.ISAROMATIC)) {

		if (!ai.getSymbol().equals("C")) {
			frags[50] += ai.getImplicitHydrogenCount();
			return;
		}

		
		List bonds = m.getConnectedBondsList(ai);

		int doublebondcount = 0;
		int triplebondcount = 0;

		String hybrid = "";

		for (int j = 0; j <= bonds.size() - 1; j++) {
			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.DOUBLE)) {
				doublebondcount++;
			}

			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.TRIPLE)) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1 && triplebondcount == 0) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		
		
		int OxNum = 0;
		int XCount = 0;

		for (int j = 0; j <= ca.size() - 1; j++) {

			String s = ((IAtom)ca.get(j)).getSymbol();

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			// System.out.println(ap.GetNormalizedElectronegativity(((IAtom)ca.get(j)).getSymbol()));

			if (this.IsElectronegative(((IAtom)ca.get(j)))) {
				List bonds2 = m.getConnectedBondsList(((IAtom)ca.get(j)));

				boolean HaveDouble = false;

				for (int k = 0; k <= bonds2.size() - 1; k++) {
					if (((IBond)bonds2.get(k)).getOrder().equals(IBond.Order.DOUBLE)) {
						HaveDouble = true;
						break;
					}

				}
				if (HaveDouble && ((IAtom)ca.get(j)).getSymbol().equals("N")
						&& ((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC))
					OxNum += 2; // C-N bond order for pyridine type N's is
				// considered to be 2
				else
					OxNum += m.getBond(ai, ((IAtom)ca.get(j))).getOrder().numeric();
			} else {
				// if (!((IAtom)ca.get(j)).getSymbol().equals("C")) System.out.println("Found
				// non electroneg:"+"\t"+((IAtom)ca.get(j)).getSymbol());

			}

			List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(j)));

			for (int k = 0; k <= ca2.size() - 1; k++) {
				String s2 = ((IAtom)ca2.get(k)).getSymbol();

				if (!s2.equals("C") && ((IAtom)ca.get(j)).getSymbol().equals("C")) {
					XCount++;
				}
			}

		}// end j loop

//		 System.out.println((i+1)+"\t"+Fragment[i]+"\t"+hybrid+"\t"+XCount+"\t"+OxNum);

		if (OxNum == 0) {
			if (hybrid.equals("sp3")) {
				if (XCount == 0) {
					frags[46] += ai.getImplicitHydrogenCount();
					return;
				} else {

					// first check for alpha carbon:
					if (ai.getSymbol().equals("C")
							&& !ai.getFlag(CDKConstants.ISAROMATIC)) {
						for (int j = 0; j <= ca.size() - 1; j++) {

							if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)
									&& ((IAtom)ca.get(j)).getSymbol().equals("C")) { // single
																		// bonded
								List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(j)));

								for (int k = 0; k <= ca2.size() - 1; k++) {

									if (!((IAtom)ca2.get(k)).getSymbol().equals("C")) {
										if (m.getBond(((IAtom)ca.get(j)), ((IAtom)ca2.get(k))).getOrder().numeric() >= 2) {
											frags[51] += ai.getImplicitHydrogenCount();
											return;
										}

										if (((IAtom)ca.get(j))
												.getFlag(CDKConstants.ISAROMATIC)
												&& ((IAtom)ca2.get(k))
														.getFlag(CDKConstants.ISAROMATIC)) {
											if (EStateFragmentDescriptor
													.InSameAromaticRing(m,
															((IAtom)ca.get(j)), ((IAtom)ca2.get(k)), rs)) {
												frags[51] += ai
														.getImplicitHydrogenCount();
												return;
											}
										}
									} // end !((IAtom)ca2.get(k)).getSymbol().equals("C"))

								} // end k loop

							} // end if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE))
								// {
						}// end j loop
					} // end if(ai.getSymbol().equals("C") &&

					if (XCount == 1) {
						frags[52] += ai.getImplicitHydrogenCount();
						return;
					} else if (XCount == 2) {
						frags[53] += ai.getImplicitHydrogenCount();
						return;
					} else if (XCount == 3) {
						frags[54] += ai.getImplicitHydrogenCount();
						return;
					} else if (XCount >= 4) {
						frags[55] += ai.getImplicitHydrogenCount();
						return;
					}
				}

			} else if (hybrid.equals("sp2")) {
				frags[47] += ai.getImplicitHydrogenCount();
				return;
			} else if (hybrid.equals("sp")) {
				frags[48] += ai.getImplicitHydrogenCount();
				return;
			}

		//TMM, 12/15/08: made sure 49 gets assigned for sp and OxNum 1-3
		} else if (hybrid.equals("sp") && OxNum>=1 && OxNum<=3) { 
			frags[49] += ai.getImplicitHydrogenCount();
			return;
		}else if (OxNum == 1) {
			if (hybrid.equals("sp3")) {
				frags[47] += ai.getImplicitHydrogenCount();
				return;
			} else if (hybrid.equals("sp2")) {
				frags[48] += ai.getImplicitHydrogenCount();
				return;
			} 
		} else if ((OxNum == 2 && hybrid.equals("sp3"))) {
			frags[48] += ai.getImplicitHydrogenCount();
			return;
		} else if ((OxNum == 3 && hybrid.equals("sp3"))
				|| (OxNum >= 2 && hybrid.equals("sp2"))) {
			frags[49] += ai.getImplicitHydrogenCount();
			return;
		} 

		// check for alpha carbon:
		if (ai.getSymbol().equals("C") && !ai.getFlag(CDKConstants.ISAROMATIC)) {
			for (int j = 0; j <= ca.size() - 1; j++) {

				if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)
						&& ((IAtom)ca.get(j)).getSymbol().equals("C")) { // single bonded
					List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(j)));

					for (int k = 0; k <= ca2.size() - 1; k++) {

						if (!((IAtom)ca2.get(k)).getSymbol().equals("C")) {
							if (m.getBond(((IAtom)ca.get(j)), ((IAtom)ca2.get(k))).getOrder().numeric() >= 2) {
								frags[51] += ai.getImplicitHydrogenCount();
								return;
							}

							if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)
									&& ((IAtom)ca2.get(k)).getFlag(CDKConstants.ISAROMATIC)) {
								if (EStateFragmentDescriptor
										.InSameAromaticRing(m, ((IAtom)ca.get(j)), ((IAtom)ca2.get(k)),
												rs)) {
									frags[51] += ai.getImplicitHydrogenCount();
									return;
								}
							}
						} // end !((IAtom)ca2.get(k)).getSymbol().equals("C"))

					} // end k loop

				} // end if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)) {
			}// end j loop
		} // end if(ai.getSymbol().equals("C") &&

	}

	private void CalcGroup056_57(int i) {
		// 56: O in =O
		// 57: O in phenol, enol, and carboxyl
		// enol : compound containing a hydroxyl group bonded to a carbon atom
		// that in turn forms a double bond with another carbon atom.
		// enol = HO-C=C- 
		// carboxyl= HO-C(=O)-

		if (!Fragment[i].equals("SsOH"))
			return;

		List ca = m.getConnectedAtomsList(m.getAtom(i));

		if (((IAtom)ca.get(0)).getFlag(CDKConstants.ISAROMATIC)) { // phenol
			frags[57]++;
			alogpfrag[i] = 57;
			return;
		}

		List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(0)));

		for (int k = 0; k <= ca2.size() - 1; k++) {

			if (m.getBond(((IAtom)ca2.get(k)), ((IAtom)ca.get(0))).getOrder().equals(IBond.Order.DOUBLE)) {

//				if ((((IAtom)ca.get(0)).getSymbol().equals("C") && ((IAtom)ca.get(0)).getHydrogenCount() == 0)) {
				if (((IAtom)ca.get(0)).getSymbol().equals("C")) {//TMM 11/13/08: removed hydrogencount=0 constraint and added S as pseudo carboxyl 
					if (((IAtom)ca2.get(k)).getSymbol().equals("O") || ((IAtom)ca2.get(k)).getSymbol().equals("S")	|| ((IAtom)ca2.get(k)).getSymbol().equals("C")) {
						frags[57]++;
						alogpfrag[i] = 57;
						return;
					}
				}

			}

		}

		frags[56]++;
		alogpfrag[i] = 56;

	}

	private void CalcGroup058_61(int i) {
		List ca = m.getConnectedAtomsList(m.getAtom(i));

		// 58: O in =O
		// 61: --O in nitro, N-oxides
		// 62: O in O-

		if (Fragment[i].equals("SsOm")) {

			if (((IAtom)ca.get(0)).getSymbol().equals("N") && ((IAtom)ca.get(0)).getFormalCharge() == 1) {
				frags[61]++;
				alogpfrag[i] = 61;
			} else {
				frags[62]++;
				alogpfrag[i] = 62;
			}

		} else if (Fragment[i].equals("SdO")) {

			if (((IAtom)ca.get(0)).getSymbol().equals("N") && ((IAtom)ca.get(0)).getFormalCharge() == 1) {
				frags[61]++;
				alogpfrag[i] = 61;
			} else {
				frags[58]++;
				alogpfrag[i] = 58;
			}
		}

	}

	private void CalcGroup059_060_063(int i) {
		// O in Al-O-Ar, Ar2O, R...O...R, ROC=X
		// ... = aromatic single bonds
		
		if (!Fragment[i].equals("SssO") && !Fragment[i].equals("SaaO")
				&& !Fragment[i].equals("SsOH"))
			return;

		// Al-O-Ar, Ar2O
		List ca = m.getConnectedAtomsList(m.getAtom(i));

		IAtom ai = m.getAtom(i);

				
		
		if (Fragment[i].equals("SssO")) {
			if (((IAtom)ca.get(0)).getFlag(CDKConstants.ISAROMATIC)
					|| ((IAtom)ca.get(1)).getFlag(CDKConstants.ISAROMATIC)) {

				
				frags[60]++;
				alogpfrag[i] = 60;
				return;
				

			} else {
				// check for ROC=X
				
				
				for (int j = 0; j <= ca.size() - 1; j++) {
					if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
						List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(j)));

						for (int k = 0; k <= ca2.size() - 1; k++) {
							if (!((IAtom)ca2.get(k)).getSymbol().equals("C")) {
								if (m.getBond(((IAtom)ca.get(j)), ((IAtom)ca2.get(k))).getOrder().equals(IBond.Order.DOUBLE)) {
									frags[60]++;
									alogpfrag[i] = 60;
									return;
								}
							}
						}
					}

				}

				
				if ((((IAtom)ca.get(0)).getSymbol().equals("O") && m.getBond(ai, ((IAtom)ca.get(0)))
						.getOrder().equals(IBond.Order.SINGLE))
						|| (((IAtom)ca.get(1)).getSymbol().equals("O") && m.getBond(ai,
								((IAtom)ca.get(1))).getOrder().equals(IBond.Order.SINGLE))) {
					frags[63]++;
					alogpfrag[i] = 63;
				} else {
					frags[59]++;
					alogpfrag[i] = 59;
				}
				

			}
		} else if (Fragment[i].equals("SaaO")) { // R...O...R
			frags[60]++;
			alogpfrag[i] = 60;
		}

	}

	private void CalcGroup064(int i) {
		if (m.getAtom(i).getSymbol().equals("Se")) {
			frags[64]++;
			alogpfrag[i]=64;
		}
		
		
	}
	
	private void CalcGroup066_to_079(int i) {

		int NAr = 0;
		int NAl = 0;
		IAtom ai = m.getAtom(i);
		if (!ai.getSymbol().equals("N"))
			return;

		List ca = m.getConnectedAtomsList(m.getAtom(i));

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC))
				NAr++;
			else
				NAl++;
		}

		// first check if have RC(=O)N or NX=X

		for (int j = 0; j <= ca.size() - 1; j++) {

			// if (Fragment[i].equals("SsNH2")) break;
			if (Fragment[i].indexOf("d") > -1)
				break; // dont assign 72 if N has double bond

			// boolean HaveCsX=false;

			if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().numeric() > 1)
				continue;

			List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(j)));

			if (((IAtom)ca.get(j)).getSymbol().equals("P"))
				continue;

			int CdXCount = 0;

			for (int k = 0; k <= ca2.size() - 1; k++) {

				//TMM 11/13/08 - commented out following if:
//				if (((IAtom) ca.get(j)).getSymbol().equals("N")
//						&& ((IAtom) ca2.get(k)).getSymbol().equals("N"))
//					continue;

				if (((IAtom) ca.get(j)).getFormalCharge() != 0)
					continue;

				if (m.getAtomNumber(((IAtom) ca2.get(k))) == i)
					continue;// 11/13/08

				// if (((((IAtom)ca.get(j)).getSymbol().equals("C") &&
				// ((IAtom)ca.get(j))
				// .getHydrogenCount() == 0) && (((IAtom)ca2.get(k)).getSymbol()
				// .equals("O") || ((IAtom)ca2.get(k)).getSymbol().equals("S")))
				// || (!((IAtom)ca.get(j)).getSymbol().equals("C") &&
				// !((IAtom)ca2.get(k))
				// .getSymbol().equals("C"))) {

				// TMM 11/13/08: removed hydrogenCount constraint
				if ((((IAtom) ca.get(j)).getSymbol().equals("C") && (((IAtom) ca2
						.get(k)).getSymbol().equals("O") || ((IAtom) ca2.get(k))
						.getSymbol().equals("S")))
						|| (!((IAtom) ca.get(j)).getSymbol().equals("C") && !((IAtom) ca2
								.get(k)).getSymbol().equals("C"))) {

					if (!((IAtom) ca2.get(k)).getFlag(CDKConstants.ISAROMATIC)
							&& !((IAtom) ca.get(j))
									.getFlag(CDKConstants.ISAROMATIC)
							&& !ai.getFlag(CDKConstants.ISAROMATIC)) {

						if (m.getBond(((IAtom) ca.get(j)),
										((IAtom) ca2.get(k))).getOrder().equals(IBond.Order.DOUBLE)) {

							CdXCount++;
						} else if (m.getBond(((IAtom) ca.get(j)),
								((IAtom) ca2.get(k))).getOrder().equals(IBond.Order.SINGLE)) {
							// HaveCsX = true;
						}
					}

				}

			}

			if (CdXCount == 1) {
				frags[72]++;
				alogpfrag[i] = 72;
				return;
			}

			// if (HaveCdX && !HaveCsX) { // later 68 might be assigned if false

		}

		if (Fragment[i].equals("SsNH2")) {

			if (((IAtom)ca.get(0)).getFlag(CDKConstants.ISAROMATIC)
					|| !((IAtom)ca.get(0)).getSymbol().equals("C")) {
				frags[69]++;
				alogpfrag[i] = 69;
			} else {
				frags[66]++;
				alogpfrag[i] = 66;
			}
		} else if (Fragment[i].equals("SaaNH") || Fragment[i].equals("SsaaN")
				|| Fragment[i].equals("SaaaN")) { // R...NH...R
			frags[73]++;
			alogpfrag[i] = 73;
		} else if (Fragment[i].equals("SssNH")) {

			if (NAr == 2 && NAl == 0) { // Ar2NH
				frags[73]++;
				alogpfrag[i] = 73;
			} else if (NAr == 1 && NAl == 1) { // Ar-NH-Al
				frags[70]++;
				alogpfrag[i] = 70;

			} else if (NAr == 0 && NAl == 2) { // Al2NH
				frags[67]++;
				alogpfrag[i] = 67;

			}

		} else if (Fragment[i].equals("SsssN")) {
			if ((NAr == 3 && NAl == 0) || (NAr == 2 && NAl == 1)) { // Ar3N &
				// Ar2NAl
				frags[73]++;
				alogpfrag[i] = 73;
			} else if (NAr == 1 && NAl == 2) {
				frags[71]++;
				alogpfrag[i] = 71;
			} else if (NAr == 0 && NAl == 3) {
				frags[68]++;
				alogpfrag[i] = 68;
			}
		} else if (Fragment[i].equals("SaaN")) {
			frags[75]++;
			alogpfrag[i] = 75;

		} else if (Fragment[i].equals("SdssNp")) {

			boolean HaveSsOm = false;
			boolean HaveSdO = false;
			boolean Ar = false;

			for (int j = 0; j <= ca.size() - 1; j++) {
				if (Fragment[m.getAtomNumber(((IAtom)ca.get(j)))].equals("SsOm")) {
					HaveSsOm = true;
				} else if (Fragment[m.getAtomNumber(((IAtom)ca.get(j)))].equals("SdO")) {
					HaveSdO = true;
				} else {
					if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)) {
						Ar = true;
					}
				}
			}

			if (HaveSsOm && HaveSdO && Ar) {
				frags[76]++;
				alogpfrag[i] = 76;
			} else if (HaveSsOm && HaveSdO && !Ar) {
				frags[77]++;
				alogpfrag[i] = 77;
			} else {
				frags[79]++;
				alogpfrag[i] = 79;
			}

		} else if (Fragment[i].equals("StN")) {
			//12/5/08, TMM: remove requirement that attached atom = C:
//			if (((IAtom)ca.get(0)).getSymbol().equals("C")) { // R#N
				frags[74]++;
				alogpfrag[i] = 74;
//			}
		} else if (Fragment[i].equals("SdNH") || Fragment[i].equals("SdsN")) {

			// test for RO-NO
			if (Fragment[i].equals("SdsN")) {
				if (((IAtom)ca.get(0)).getSymbol().equals("O")
						&& ((IAtom)ca.get(1)).getSymbol().equals("O")) {
					frags[76]++;
					alogpfrag[i] = 76;
					return;
				}
			}

			boolean flag1 = false;
			boolean flag2 = false;

			for (int j = 0; j <= ca.size() - 1; j++) {
				if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.DOUBLE)) {
					if (((IAtom)ca.get(j)).getSymbol().equals("C")
							|| (((IAtom)ca.get(j)).getSymbol().equals("N") && Fragment[i]
									.equals("SdNH"))) {
						frags[74]++;
						alogpfrag[i] = 74;
						return;
					} else {
						flag1 = true;
					}
				} else {
					if (!((IAtom)ca.get(j)).getSymbol().equals("C")
							|| ((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)) {
						flag2 = true;
					}
				}

				if (flag1 && flag2) { // X-N=X or Ar-N=X
					frags[78]++;
					alogpfrag[i] = 78;
					return;
				} else if (flag1) {
					// System.out.println("missing group: R-N=X");
					// assign HN=X and RN=X to 78:
					frags[78]++;
					alogpfrag[i] = 78;
					return;

				}

			}

		} else if (Fragment[i].indexOf("p") > -1) {
			//TMM: 11/14/08
			if (Fragment[i].equals("SddNp")) { //part of C=[N+]=[N-]
				frags[74]++;
				alogpfrag[i] = 74;
			} else {
				if (!Fragment[i].equals("SsaaNp")) {
					frags[79]++;
					alogpfrag[i] = 79;
				}
			}
//			System.out.println(Fragment[i]);
		}
		
//		TMM: 11/14/08
		if (Fragment[i].equals("SdNm")) {//part of C=[N+]=[N-]
			frags[74]++;
			alogpfrag[i] = 74;
		}
		
		
		// R--N(--R)--O
		if (Fragment[i].equals("SdaaN")) {
			
//			System.out.println("here");
			for (int j=0;j<=ca.size()-1;j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("O") && m.getBond(ai,((IAtom)ca.get(j))).getOrder().equals(IBond.Order.DOUBLE)) {
					boolean InSameRing=EStateFragmentDescriptor.InSameAromaticRing(m, ai,(IAtom)ca.get(j), rs);
					if (!InSameRing) {
						frags[76]++;
						alogpfrag[i] = 76;
						break;
					}
				}
			}
		}
//		TMM: 11/14/08: if have charged version of pyridine-N-oxide
		if (Fragment[i].equals("SsaaNp")) {
			for (int j=0;j<=ca.size()-1;j++) {
				IAtom caj=(IAtom)ca.get(j);
				int num=m.getAtomNumber(caj);
				if (Fragment[num].equals("SsOm")) {
						frags[76]++;
						alogpfrag[i] = 76;
						break;
				}
			}
		}


	}

	private void CalcGroup081_to_085(int i) {
		//TMM-11/20/08: following requirement commented out:
//		if (!Fragment[i].equals("SsF"))
//			return;

		IAtom ai = m.getAtom(i);
		if (!ai.getSymbol().equals("F")) return;//added TMM-11/20/08
		
		List ca = m.getConnectedAtomsList(m.getAtom(i));

		List bonds = m.getConnectedBondsList(((IAtom)ca.get(0)));

		int doublebondcount = 0;
		int triplebondcount = 0;

		String hybrid = "";

		for (int j = 0; j <= bonds.size() - 1; j++) {
			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.DOUBLE)) {
				doublebondcount++;
			}

			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.TRIPLE)) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(0)));

		int OxNum = 0;

		for (int k = 0; k <= ca2.size() - 1; k++) {

			String s = ((IAtom)ca2.get(k)).getSymbol();

			// // F,O,Cl,Br,N

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			if (this.IsElectronegative(((IAtom)ca2.get(k)))) {
				OxNum += m.getBond(((IAtom)ca.get(0)), ((IAtom)ca2.get(k))).getOrder().numeric();
			}

		}

		if ((hybrid.equals("sp2") && OxNum > 1)
				|| (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4)
				|| !((IAtom)ca.get(0)).getSymbol().equals("C")) {
			frags[85]++;
			alogpfrag[i] = 85;
		} else if (hybrid.equals("sp3") && OxNum == 1) {
			frags[81]++;
			alogpfrag[i] = 81;
		} else if (hybrid.equals("sp3") && OxNum == 2) {
			frags[82]++;
			alogpfrag[i] = 82;
		} else if (hybrid.equals("sp3") && OxNum == 3) {
			frags[83]++;
			alogpfrag[i] = 83;
		} else if (hybrid.equals("sp2") && OxNum == 1) {
			frags[84]++;
			alogpfrag[i] = 84;
		}

	}

	private void CalcGroup086_to_090(int i) {
		//TMM-11/20/08: following requirement commented out:
//		if (!Fragment[i].equals("SsCl"))
//			return;
		IAtom ai = m.getAtom(i);
		if (!ai.getSymbol().equals("Cl")) return;//added TMM-11/20/08

		List ca = m.getConnectedAtomsList(m.getAtom(i));

		List bonds = m.getConnectedBondsList(((IAtom)ca.get(0)));

		int doublebondcount = 0;
		int triplebondcount = 0;

		String hybrid = "";

		for (int j = 0; j <= bonds.size() - 1; j++) {
			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.DOUBLE)) {
				doublebondcount++;
			}

			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.TRIPLE)) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(0)));

		int OxNum = 0;

		for (int k = 0; k <= ca2.size() - 1; k++) {

			String s = ((IAtom)ca2.get(k)).getSymbol();

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			if (this.IsElectronegative(((IAtom)ca2.get(k)))) {
				// // F,O,Cl,Br,N
				OxNum += m.getBond(((IAtom)ca.get(0)), ((IAtom)ca2.get(k))).getOrder().numeric();
			}
		}

		if ((hybrid.equals("sp2") && OxNum > 1)
				|| (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4)
				|| !((IAtom)ca.get(0)).getSymbol().equals("C")) {
			frags[90]++;
			alogpfrag[i] = 90;
		} else if (hybrid.equals("sp3") && OxNum == 1) {
			frags[86]++;
			alogpfrag[i] = 86;
		} else if (hybrid.equals("sp3") && OxNum == 2) {
			frags[87]++;
			alogpfrag[i] = 87;
		} else if (hybrid.equals("sp3") && OxNum == 3) {
			frags[88]++;
			alogpfrag[i] = 88;
		} else if (hybrid.equals("sp2") && OxNum == 1) {
			frags[89]++;
			alogpfrag[i] = 89;
		}

	}

	private void CalcGroup091_to_095(int i) {
		//TMM-11/20/08: following requirement commented out:
//		if (!Fragment[i].equals("SsBr"))
//			return;
		IAtom ai = m.getAtom(i);
		if (!ai.getSymbol().equals("Br")) return;//added TMM-11/20/08


		List ca = m.getConnectedAtomsList(m.getAtom(i));

		List bonds = m.getConnectedBondsList(((IAtom)ca.get(0)));

		int doublebondcount = 0;
		int triplebondcount = 0;

		String hybrid = "";

		for (int j = 0; j <= bonds.size() - 1; j++) {
			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.DOUBLE)) {
				doublebondcount++;
			}

			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.TRIPLE)) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(0)));

		int OxNum = 0;

		for (int k = 0; k <= ca2.size() - 1; k++) {

			String s = ((IAtom)ca2.get(k)).getSymbol();

			// // F,O,Cl,Br,N

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			if (this.IsElectronegative(((IAtom)ca2.get(k)))) {
				OxNum += m.getBond(((IAtom)ca.get(0)), ((IAtom)ca2.get(k))).getOrder().numeric();
			}

		}

		if ((hybrid.equals("sp2") && OxNum > 1)
				|| (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4)
				|| !((IAtom)ca.get(0)).getSymbol().equals("C")) {
			frags[95]++;
			alogpfrag[i] = 95;
		} else if (hybrid.equals("sp3") && OxNum == 1) {
			frags[91]++;
			alogpfrag[i] = 91;
		} else if (hybrid.equals("sp3") && OxNum == 2) {
			frags[92]++;
			alogpfrag[i] = 92;
		} else if (hybrid.equals("sp3") && OxNum == 3) {
			frags[93]++;
			alogpfrag[i] = 93;
		} else if (hybrid.equals("sp2") && OxNum == 1) {
			frags[94]++;
			alogpfrag[i] = 94;
		}

	}

	private void CalcGroup096_to_100(int i) {

		//TMM-11/20/08: following requirement commented out:
//		if (!Fragment[i].equals("SsI"))
//			return;
		IAtom ai = m.getAtom(i);
		if (!ai.getSymbol().equals("I")) return;//added TMM-11/20/08

		List ca = m.getConnectedAtomsList(m.getAtom(i));

		List bonds = m.getConnectedBondsList(((IAtom)ca.get(0)));

		int doublebondcount = 0;
		int triplebondcount = 0;

		String hybrid = "";

		for (int j = 0; j <= bonds.size() - 1; j++) {
			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.DOUBLE)) {
				doublebondcount++;
			}

			if (((IBond)bonds.get(j)).getOrder().equals(IBond.Order.TRIPLE)) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		List ca2 = m.getConnectedAtomsList(((IAtom)ca.get(0)));

		int OxNum = 0;

		for (int k = 0; k <= ca2.size() - 1; k++) {

			String s = ((IAtom)ca2.get(k)).getSymbol();

			// // F,O,Cl,Br,N

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			if (this.IsElectronegative(((IAtom)ca2.get(k)))) {
				OxNum += m.getBond(((IAtom)ca.get(0)), ((IAtom)ca2.get(k))).getOrder().numeric();
			}

		}

		if ((hybrid.equals("sp2") && OxNum > 1)
				|| (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4)
				|| !((IAtom)ca.get(0)).getSymbol().equals("C")) {
			frags[100]++;
			alogpfrag[i] = 100;
		} else if (hybrid.equals("sp3") && OxNum == 1) {
			frags[96]++;
			alogpfrag[i] = 96;
		} else if (hybrid.equals("sp3") && OxNum == 2) {
			frags[97]++;
			alogpfrag[i] = 97;
		} else if (hybrid.equals("sp3") && OxNum == 3) {
			frags[98]++;
			alogpfrag[i] = 98;
		} else if (hybrid.equals("sp2") && OxNum == 1) {
			frags[99]++;
			alogpfrag[i] = 99;
		}

	}

	private void CalcGroup101_to_104(int i) {
		IAtom ai = m.getAtom(i);
		String s = ai.getSymbol();

		if (ai.getFormalCharge() == -1) {
			if (s.equals("F")) {
				frags[101]++;
				alogpfrag[i] = 101;
			} else if (s.equals("Cl")) {
				frags[102]++;
				alogpfrag[i] = 102;
			} else if (s.equals("Br")) {
				frags[103]++;
				alogpfrag[i] = 103;
			} else if (s.equals("I")) {
				frags[104]++;
				alogpfrag[i] = 104;
			}

		}

	}

	private void CalcGroup106(int i) {

		// S in SH
		if (Fragment[i].equals("SsSH")) {
			frags[106]++;
			alogpfrag[i] = 106;
		}
	}

	private void CalcGroup107(int i) {

		// S in R2S, RS-SR
		// R = any group linked through C
		// if (!Fragment[i].equals("SssS")) return;

		// In ALOGP, for malathion PSC is consider to have group 107 (even
		// though has P instead of R)

		// for lack of fragment, use this fragment for SaaS
		// for 000000-28-1, dragon assigns 107 for SssssssS
		
		//TMM 11/14/08 added SssssssS to following if statement:
		if (Fragment[i].equals("SssS") || Fragment[i].equals("SaaS") || Fragment[i].equals("SssssssS")) {
			frags[107]++;
			alogpfrag[i] = 107;
		}

		// IAtom [] ca=m.getConnectedAtomsList(m.getAtom(i));
		//		
		// if ((((IAtom)ca.get(0)).getSymbol().equals("C") && ((IAtom)ca.get(1)).getSymbol().equals("C"))
		// ||
		// (((IAtom)ca.get(0)).getSymbol().equals("C") && ((IAtom)ca.get(1)).getSymbol().equals("S")) ||
		// (((IAtom)ca.get(0)).getSymbol().equals("S") && ((IAtom)ca.get(1)).getSymbol().equals("C"))) {
		// frags[107]++;
		// alogpfrag[i]=107;
		// }

	}

	private void CalcGroup108(int i) {

		// S in R=S

		// In ALOGP, for malathion P=S is consider to have group 108 (even
		// though has P instead of R)
		if (Fragment[i].equals("SdS")) {
			frags[108]++;
			alogpfrag[i] = 108;
		}

	}

	private void CalcGroup109(int i) {

		// for now S in O-S(=O)-O is assigned to this group
		// (it doesn't check which atoms are singly bonded to S
		if (!Fragment[i].equals("SdssS"))
			return;

		List ca = m.getConnectedAtomsList(m.getAtom(i));
		IAtom ai = m.getAtom(i);
		int SdOCount = 0;
		int SsCCount = 0;

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
					SsCCount++;
				}
			} else if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.DOUBLE)) {
				if (((IAtom)ca.get(j)).getSymbol().equals("O")) {
					SdOCount++;
				}
			}
		}
		if (SdOCount == 1) { // for now dont check if SsCCount==2
			frags[109]++;
			alogpfrag[i] = 109;
		}
	}

	private void CalcGroup110(int i) {
		if (!Fragment[i].equals("SddssS"))
			return;

		List ca = m.getConnectedAtomsList(m.getAtom(i));
		IAtom ai = m.getAtom(i);
		int SdOCount = 0;
		int SsCCount = 0;

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
					SsCCount++;
				}
			} else if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.DOUBLE)) {
				if (((IAtom)ca.get(j)).getSymbol().equals("O")) {
					SdOCount++;
				}
			}
		}
		if (SdOCount == 2) { // for now dont check if SsCCount==2
			frags[110]++;
			alogpfrag[i] = 110;
		}

	}

	private void CalcGroup111(int i) {

		if (m.getAtom(i).getSymbol().equals("Si")) {
			// if (Fragment[i].equals("SssssSi")) {
			frags[111]++;
			alogpfrag[i] = 111;
		}
	}

	//TMM: added 11/20/08
	private void CalcGroup112(int i) {

		if (m.getAtom(i).getSymbol().equals("B")) {//for now assign all B's to this frag:
			frags[112]++;
			alogpfrag[i] = 112;
		}
	}
	private void CalcGroup116_117_120(int i) {

		// S in R=S

		List ca = m.getConnectedAtomsList(m.getAtom(i));
		IAtom ai = m.getAtom(i);

		int XCount = 0;
		int RCount = 0;
		boolean PdX = false;

		//TMM 11/13/08 : added !Fragment[i].equals("SdssPH")
		if (!Fragment[i].equals("SdsssP") && !Fragment[i].equals("SdssPH"))
			return;

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
					RCount++;
				} else {
					XCount++;
				}
			} else if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.DOUBLE)) {
				if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
					PdX = true;
				}
			}
		}

		if (PdX) {
			//TMM-11/20/08 redone following
			if (XCount == 3) {
				frags[117]++;
				alogpfrag[i] = 117;
			} else if (XCount == 2) {
				frags[120]++;
				alogpfrag[i] = 120;
			} else if (XCount < 2) { // see 2074-67-1 and 16391-07-4 for cases with PH
				frags[116]++;
				alogpfrag[i] = 116;
			}

		}

	}

	private void CalcGroup118_119(int i) {
		if (!Fragment[i].equals("SsssP"))
			return;

		List ca = m.getConnectedAtomsList(m.getAtom(i));
		IAtom ai = m.getAtom(i);
		int XCount = 0;
		int RCount = 0;

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (m.getBond(ai, ((IAtom)ca.get(j))).getOrder().equals(IBond.Order.SINGLE)) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
					RCount++;
				} else {
					XCount++;
				}
			}
		}

		if (XCount == 3) {
			frags[118]++;
			alogpfrag[i] = 118;
		} else if (RCount == 3) {
			frags[119]++;
			alogpfrag[i] = 119;
		}

	}

	public void WriteFragTable(FileWriter fw) {
		try {
			fw
					.write("<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">\n");

			fw.write("<caption>ALOGP fragments</caption>\n");

			fw.write("<tr bgcolor=\"#D3D3D3\">\n");
			fw.write("\t<th>Fragment</th>\n");
			fw.write("\t<th>Description</th>\n");
			fw.write("\t<th>Count</th>\n");
			fw.write("\t<th>Fragval</th>\n");
			fw.write("</tr>\n");

			for (int i = 1; i <= 120; i++) {
				if (frags[i] > 0) {

					fw.write("<tr>\n");
					fw.write("\t<td>" + i + "</td>\n");
					fw.write("\t<td>" + fragdesc[i] + "</td>\n");
					fw.write("\t<td>" + frags[i] + "</td>\n");
					fw.write("\t<td>" + fragval[i] + "</td>\n");
					fw.write("</tr>\n");

				}
			}
			fw.write("</table>\n");

		} catch (Exception e) {

		}

	}

	private boolean IsElectronegative(IAtom a) {
		if (!a.getSymbol().equals("C"))
			return true;
		// if (ap.GetNormalizedElectronegativity(a.getSymbol()) > 1 ||
		// a.getSymbol().equals("Si")) return true;
		else
			return false;
	}

	// public void figureoutwhichfragments(double matchval) {
	//		
	// double [] Cfragval=new double [44];
	// double [] Hfragval=new double [10];
	// double [] Ofragval=new double [8];
	//		
	// double alogp=0;
	//		
	// for (int i=0;i<=43;i++) {
	// Cfragval[i]=fragval[i+1];
	// }
	//		
	// for (int i=0;i<=9;i++) {
	// Hfragval[i]=fragval[i+46];
	// }
	//		
	// for (int i=0;i<=7;i++) {
	// Ofragval[i]=fragval[i+56];
	// }
	//		
	// int [][] choices=new int [m.getAtomCount()][44];
	//		
	// for (int i = 0; i <= 43; i++) {
	// for (int j = 0; j <= 9; j++) {
	// for (int k = 0; k <= 9; k++) {
	// for (int l= 0; l <= 7; l++) {
	//						for (int m=0;m<=7;m++) {
	//							System.out.println(i+"\t"+j+"\t"+k+"\t"+l+"\t"+m+"\t"+alogp);
	//						}						
	//					}					
	//				}
	//			}
	//		}
	//		
	//		
	//	}

}// end class

