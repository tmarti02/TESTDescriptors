import static org.junit.Assert.assertEquals;

import java.io.InputStream;
import java.util.Scanner;
import java.util.Vector;

import org.junit.Test;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorsFromSmiles;
public class DescriptorsTest {

	@Test
	public void compareToTESTGUI() throws Exception {	
		
		try {
			
			InputStream ins=this.getClass().getClassLoader().getResourceAsStream("descriptors input.tsv");
			InputStream ins2=this.getClass().getClassLoader().getResourceAsStream("descriptors output TEST5.1.tsv");
			
			
			Scanner scanner=new Scanner(ins);
			Scanner scanner2=new Scanner(ins2);
			String headerGUI=scanner2.nextLine();
			
			while (scanner.hasNext()) {
				String Line=scanner.nextLine();
				String Line2=scanner2.nextLine();
				
				String []vals=Line.split("\t");
				
				DescriptorData dd=DescriptorsFromSmiles.goDescriptors(vals[0], vals[1], null);
				
				Vector<String>calcDesc=dd.getDescriptorValues();
				String [] guiDescVector=Line2.split("\t");
				
				for (int i=0;i<calcDesc.size();i++) {
					double val1=Double.parseDouble(calcDesc.get(i));
					double val2=Double.parseDouble(guiDescVector[i+1]);
//					System.out.println(val1+"\t"+val2);
					assertEquals(val1, val2,1e-5);
				}
				
				
			}
			
			
			scanner.close();
			scanner2.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		
		
		//TODO
		assertEquals(0, 0);
		
	}
	
	
}

