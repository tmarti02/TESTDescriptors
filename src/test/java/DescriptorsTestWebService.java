import static org.junit.Assert.assertEquals;

import java.io.InputStream;
import java.util.Hashtable;
import java.util.Scanner;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import org.json.JSONObject;
import org.junit.Test;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.async.Callback;
import com.mashape.unirest.http.exceptions.UnirestException;


public class DescriptorsTestWebService {

	public static Hashtable<String,String> htDescriptors=new Hashtable<>();
	
	static int count=0;
	
	@Test
	public void compareToTESTGUI() throws Exception {	

		try {

			InputStream ins=this.getClass().getClassLoader().getResourceAsStream("descriptors input.tsv");
			InputStream ins2=this.getClass().getClassLoader().getResourceAsStream("descriptors output TEST5.1.tsv");

			Scanner scanner=new Scanner(ins);
			Scanner scanner2=new Scanner(ins2);
			String headerGUI=scanner2.nextLine();

			Unirest.setTimeouts(10000, 10000);
						
			while (scanner.hasNext()) {
				String Line=scanner.nextLine();
				String Line2=scanner2.nextLine();

				String []vals=Line.split("\t");

				String smiles=vals[0];
				String CID=vals[1];
				
				if (smiles.contains("\\")) continue;
				

//				Thread.sleep(100);
				
				runRequest(smiles,CID);	
//				String result=waitUntilStringFoundInResponse(smiles, 100);
				
				count++;
				
//				System.out.println("****"+"\t"+count+"\t"+CID+"\t"+result);

//				//				System.out.println(response.getBody());
//				System.out.println("***descriptors="+strDescriptors);
//
//
//				//				String strDescriptors=DescriptorsFromSmiles.goDescriptors(vals[0], null,null);				
//
//				String [] newVals=strDescriptors.split("\t");
//				String [] guiVals=Line2.split("\t");
//
//				for (int i=0;i<newVals.length;i++) {
//					double val1=Double.parseDouble(newVals[i]);
//					double val2=Double.parseDouble(guiVals[i+1]);
//					//					System.out.println(val1+"\t"+val2);
//					assertEquals(val1, val2,1e-5);
//				}
				//				if (true) break;

			}


			scanner.close();
			scanner2.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}


		//TODO
		assertEquals(0, 0);

	}

	String runRequest (String smiles,String CID) {
		
		String result=null;

		Future<HttpResponse<JsonNode>> future = Unirest.get("http://localhost:8080/TEST-descriptors/calculation?smiles="+smiles)		        
				.asJsonAsync(new Callback<JsonNode>() {

					public void completed(HttpResponse<JsonNode> response) {
						JSONObject myObj = response.getBody().getObject();
						String strDescriptors=myObj.getString("descriptorValues");
						htDescriptors.put(CID, strDescriptors);
						System.out.println(++count+"\t"+CID+"\t"+strDescriptors);
					}

					public void cancelled() {
						// Do something if the request is cancelled
					}

					@Override
					public void failed(UnirestException e) {
						System.out.println(CID+"\t***FAILED");

					}

				});
		return null;

	}
	
	public String waitUntilStringFoundInResponse(String smiles, int TIMEOUT) throws Exception {
	    String result = null;
	    int i = 0;
	        
	    	HttpResponse<JsonNode> response = Unirest.get("http://localhost:8080/TEST-descriptors/calculation?smiles="+smiles)
	    			  .asJson();

	    	while (i < TIMEOUT) {
	    	
			try {
				JSONObject myObj = response.getBody().getObject();
				result=myObj.getString("descriptorValues");
				return result;
			} catch (Exception ex) {
				
			}
			
            TimeUnit.SECONDS.sleep(1);
            ++i;
            
            System.out.println(i);
            
            if (i == TIMEOUT) {
            	throw new TimeoutException("Timed out after waiting for " + i + " seconds");
            }
	        
	    }
	    return result;
	}

}

