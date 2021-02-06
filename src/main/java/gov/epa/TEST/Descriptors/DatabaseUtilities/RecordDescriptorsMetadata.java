package gov.epa.TEST.Descriptors.DatabaseUtilities;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public class RecordDescriptorsMetadata {

	public String DescriptorSoftware;
	public String DescriptorNames;
	
	public static final String fieldNameDescriptorSoftware="DescriptorSoftware";
	public static final String fieldNameDescriptorNames="DescriptorNames";
	public static final String tableNameMetaData="MetaData";
	
	public static final String softwareTEST="T.E.S.T. 5.1";
	public static final String softwarePadel="Padelpy webservice single";
	
	public static String [] getFields() {
		String[] strFields = {fieldNameDescriptorSoftware,fieldNameDescriptorNames};
		return strFields;
	}
	
	public static String getDescriptorNames(Statement stat,String software) {
		
		ResultSet rsHeader=SQLite_GetRecords.getRecords(stat, RecordDescriptorsMetadata.tableNameMetaData, fieldNameDescriptorSoftware, software);
											
	
		try {
			return rsHeader.getString(2);
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}

	}
	
	
	public void addRecord(Connection conn) {				
		String [] values= {DescriptorSoftware,DescriptorNames};				
		SQLite_CreateTable.addDataToTable(tableNameMetaData, getFields(), values,conn);
	}

	
	private static void addDataRecord(Connection conn,String DescriptorSoftware,String DescriptorNames) {				
		String [] values= {DescriptorSoftware,DescriptorNames};				
		SQLite_CreateTable.addDataToTable(tableNameMetaData, getFields(), values,conn);
	}
	
	public static ResultSet getMetaDataRecords(Statement stat,String descriptorSoftware) {

		try {
			String query="select * from "+tableNameMetaData+" where "+fieldNameDescriptorSoftware+"=\""+descriptorSoftware+"\";";
			ResultSet rs = stat.executeQuery(query);
			return rs;
		} catch (Exception ex) {
			ex.printStackTrace();
			return null;
		}
	}
	
	
	//TODO move to different table in db
		public static void addHeaderRecord(Connection conn,String descriptorSoftware,String descriptorNames) {
						
			try {
				Statement stat=conn.createStatement();
				
				ResultSet rs=getMetaDataRecords(stat, descriptorSoftware);
							
				if (!rs.next()) {
					addDataRecord(conn,descriptorSoftware,descriptorNames);
				}
			
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
}

