package gov.epa.TEST.Descriptors.DescriptorCalculations;


public class KowDLL {

	/**
	 * @param args
	 */
	
	double KowExp;
	
	public double Calculate (String Smiles) {
	
		KowExp=-9999;
		
		
		String ChemicalName="bob";
		float EstValue=-9999;		
		float ExpKowVal=-9999;
		String DetailResults="";
		int numLines=-9999;
		String ErrorMess="";
		
		String result="";
		
		try {
			result=this.GetSrcKow32(Smiles,ChemicalName, EstValue, ExpKowVal,DetailResults,numLines, ErrorMess);
//			System.out.println(result);
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
//		ArrayList al=new ArrayList();
		if (result.equals("")) return -9999;
		
		String s=result;
		
		String strKowCalc="";
		
		while (s.indexOf("\n")>-1) {
			String line=s.substring(0,s.indexOf("\n"));
//			System.out.println("*"+line);
//			al.add(line);
			
			if (line.indexOf("Log Kow")>-1) {
				strKowCalc=line.substring(line.indexOf("=")+1,s.indexOf("\n"));
				strKowCalc=strKowCalc.trim();				
			}
			
			if (line.indexOf("Exp Log P:")>-1) {
				String strKowExp=line.substring(line.indexOf(":")+1,line.length());
				strKowExp=strKowExp.trim();
//				System.out.println("Expval="+strKowExp);
			}
			
			
			s=s.substring(s.indexOf("\n")+1,s.length());
		}
		
		double KowCalc=Double.parseDouble(strKowCalc);
//		System.out.println(KowCalc);
		
		return KowCalc;
		
	}
	
	
	
	public native String GetSrcKow32(String cSmilePass,
		       String cChemical, float EstValue, float ExpKowVal,
		       String DetailResults, int numLines, String ErrorMess);
		
	static
	{
	 	System.loadLibrary("KOWWrapper");  
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		KowDLL k=new KowDLL();
		String Smiles = "CCCCO";
		k.Calculate(Smiles);
		
				
	}

}
