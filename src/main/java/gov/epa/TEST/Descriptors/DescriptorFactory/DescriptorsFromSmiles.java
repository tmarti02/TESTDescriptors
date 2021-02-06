package gov.epa.TEST.Descriptors.DescriptorFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import gov.epa.TEST.Descriptors.DatabaseUtilities.DescriptorDatabaseUtilities;
import gov.epa.TEST.Descriptors.DatabaseUtilities.RecordDescriptors;
import gov.epa.TEST.Descriptors.DatabaseUtilities.RecordDescriptorsMetadata;
import gov.epa.TEST.Descriptors.DatabaseUtilities.SQLite_CreateTable;
import gov.epa.TEST.Descriptors.DatabaseUtilities.SQLite_Utilities;
import gov.epa.TEST.Descriptors.DescriptorUtilities.HueckelAromaticityDetector;

import java.io.File;
import java.io.FileWriter;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashMap;


public class DescriptorsFromSmiles {
	public static boolean dashboardStructuresAvailable=true;



	private static final Logger logger = LogManager.getLogger(DescriptorsFromSmiles.class);


	private static void writeJSON(String filePath,Object object) {
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		//		builder.disableHtmlEscaping();
		Gson gson = builder.create();
		try {
			FileWriter fw=new FileWriter(filePath);
			fw.write(gson.toJson(object));
			fw.close();

		} catch (Exception ex) {
			logger.error("Error writing "+filePath);
		}
	}

	private static void printJSON(Object object) {
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		//		builder.disableHtmlEscaping();
		Gson gson = builder.create();
		System.out.println(gson.toJson(object));
	}

	

	public static String goDescriptors(String smiles,String DSSTox_Structure_Id,Connection conn) {


		String descriptorSoftware=RecordDescriptorsMetadata.softwareTEST;

		AtomContainer ac=loadSMILES(smiles);
		
		String error = (String) ac.getProperty("Error");
		HueckelAromaticityDetector.debug=false;//suppress logging

		try {
			
			if (conn!=null) {
				String strDescriptors=RecordDescriptors.lookupDescriptorsInDatabase(conn,DSSTox_Structure_Id,descriptorSoftware);
				
				if (strDescriptors!=null) {
//					System.out.println("Found in db:"+strDescriptors);
					return strDescriptors;
				}
			}

			int descresult=-1;

			DescriptorData dd=new DescriptorData();


			long t1 = System.currentTimeMillis();

			DescriptorFactory df = new DescriptorFactory(false);
			df.Calculate3DDescriptors=false;

			if(error.isEmpty()) {
				descresult = df.CalculateDescriptors(ac, dd, true);
			}

			long t2 = System.currentTimeMillis();

			logger.debug("Calculated descriptors in {}s", (t2 - t1) / 1000.);

			if (descresult == -1) {
				dd.Error=df.errorMsg;//store error message in dd
				//TODO should we have stored error in dd to begin with?

				logger.error("Error processing record {}: {}", DSSTox_Structure_Id,dd.Error);
				ac.setProperty("Error", dd.Error);
				ac.setProperty("ErrorCode", CheckAtomContainer.ERROR_CODE_DESCRIPTOR_CALCULATION_ERROR);
			} else {
				dd.Error = "OK";
			}		

			if (error.isEmpty()) {
				if(df.errorMsg.isEmpty()) dd.Error="OK";
				else dd.Error=df.errorMsg;
			} else {
				dd.Error=error;
			}

			//Store in db:
							
			RecordDescriptors r=new RecordDescriptors();
			r.Descriptors=dd.toStringDescriptorValues();
			r.DescriptorSoftware=descriptorSoftware;
			r.DSSTox_Structure_Id=DSSTox_Structure_Id;
			
//			System.out.println(r.Descriptors);
						
			if (conn!=null) r.addRecord(conn);

			////////////////////////////////////////
			return r.Descriptors;

		} catch (Exception ex) {
			logger.error("Error processing record {}: {} {}", DSSTox_Structure_Id, ex.getClass(), ex.getMessage());
			logger.catching(ex);
			DescriptorData dd = new DescriptorData();
			dd.Error = "Error processing record " + DSSTox_Structure_Id + ", error=" + ex.getMessage();
			ac.setProperty("Error", dd.Error);
			ac.setProperty("ErrorCode", CheckAtomContainer.ERROR_CODE_DESCRIPTOR_CALCULATION_ERROR);

			return "Error:"+dd.Error;
		}

	}

	
	
	public static String goDescriptors(String smiles) {
		return goDescriptors(smiles, null, null);
	}


	

	
	private static AtomContainer prepareSmilesMolecule(String Smiles) {

		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

		AtomContainer m = null;
		try {
			// m=sp.parseSmiles(Smiles);
			m = (AtomContainer) sp.parseSmiles(Smiles);



			m.setProperty("Error", "");

			if (Smiles.indexOf(".") > -1) {
				m.setProperty("Error", "Molecules can only contain one fragment");
				m.setProperty("ErrorCode", CheckAtomContainer.ERROR_CODE_APPLICABILITY_DOMAIN_ERROR);
			} else {
				CheckAtomContainer.checkAtomContainer(m);	
			}


		} catch (org.openscience.cdk.exception.InvalidSmilesException e) {
			m = new AtomContainer();
			String error = e.getMessage() + ", SMILES=" + Smiles;
			m.setProperty("Error", error);
			m.setProperty("ErrorCode", CheckAtomContainer.ERROR_CODE_STRUCTURE_ERROR);
			logger.error(error);
		}
		return m;
	}

	public static AtomContainer loadSMILES(String smiles) {

		AtomContainer m = prepareSmilesMolecule(smiles);

		String ID;

		if (m.getProperty("ErrorCode")==CheckAtomContainer.ERROR_CODE_STRUCTURE_ERROR) {
			ID = "C_" + System.currentTimeMillis();
			m.setProperty("CAS", ID);
			return m;
		}

		return m;
	}

	
	

	String convertDescriptorDataToLine(DescriptorData dd) {
	
		LinkedHashMap<String, String>ddMap=dd.convertToLinkedHashMap();
	
		String Error=dd.Error;
	
		if (!Error.contentEquals("OK")) {
			return dd.ID+"\t"+Error;			
		}
	
		String Line="";
	
	
		for( String key : ddMap.keySet() ){
	
			if (key.contentEquals("Index")) continue;
			if (key.contentEquals("Query")) continue;
			if (key.contentEquals("SmilesRan")) continue;
			if (key.contentEquals("Error")) continue;
	
			String value = ddMap.get(key);
			Line+=value+"\t";		  
			//		System.out.println(key+"\t"+value);
		}
		Line=Line.trim();
	
	
		return Line;
	
	}
	
	

	
	public static void main(String[] args) {
//		String smiles="OC(=O)C=C";
//		String inchiKey1="NIXOWILDQLNWCW";
		
		String smiles="OC(=O)C=CXXX";
		String inchiKey1="XXXXXXXXXXXX";


		String strDescriptors=goDescriptors(smiles,"CID123",null);
					
//		System.out.println("Vals returned:"+strDescriptors);
			
	}


}
