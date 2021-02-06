package gov.epa.TEST.Descriptors.DatabaseUtilities;

import java.nio.charset.Charset;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


public class DescriptorDatabaseUtilities {
	
	private static final Logger logger = LogManager.getLogger(DescriptorDatabaseUtilities.class);
	
	static final String mySQL_DB_URL = "jdbc:mysql://localhost/";
	static final String USER = "myuser";
	static final String PASS = "Epa778899";
		
//	public static final String filepathDB="databases/descriptors.db";
	public static final String filepathDB="databases/datasets.db";
	
//	public static final String fieldNameKey="InChIKey1QSARReady";
			
	
	public static void createDescriptorsDB(Connection connDescriptors) throws SQLException {
		Statement stat=connDescriptors.createStatement();
		SQLite_CreateTable.create_table(stat,RecordDescriptors.tableNameDescriptors,RecordDescriptors.getFields());
		RecordDescriptors.createIndex(stat);	
		SQLite_CreateTable.create_table(stat,RecordDescriptorsMetadata.tableNameMetaData,RecordDescriptorsMetadata.getFields());
	}
		
	public static Connection getConnection (String descriptorSoftware,String descriptorNames,String filepathDB) {
		Connection connDescriptors=null;

		try {
			connDescriptors=SQLite_Utilities.getConnection(filepathDB);
			createDescriptorsDB(connDescriptors);		
			RecordDescriptorsMetadata.addHeaderRecord(connDescriptors,descriptorSoftware,descriptorNames);//TODO move to different table
			
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return connDescriptors;
		
	}

	

	
//	/**
//	 * For mysql db
//	 * 
//	 * @return
//	 */
//	private static int [] getFieldLengthsDescriptors() {
//		
//		ArrayList<Integer> fields=new ArrayList<Integer>();
//		fields.add(15);//CAS
////		fields.add("InChi");
//		fields.add(30);//InChiKey
//		fields.add(25000);//Descriptors- should we zip it?
//		fields.add(100);//error
//				
//		int[] intLength = new int[fields.size()];
//		
//		for (int i=0;i<fields.size();i++) {
//			intLength[i]=fields.get(i);
//		}
//		
//		return intLength;
//
//	}
		
	

	
	
}

