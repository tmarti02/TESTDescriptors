package gov.epa.TEST.Descriptors.DatabaseUtilities;

import java.nio.charset.Charset;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RecordDescriptors {
	public static final String tableNameDescriptors="Descriptors";
	
	public String DescriptorSoftware;
	public String DSSTox_Structure_Id;
	public String Descriptors;
	
	public static final String fieldNameSoftware="DescriptorSoftware";
	public static final String fieldNameKey="DSSTox_Structure_Id";
	public static final String fieldNameDescriptors="Descriptors";
	
	public static final String fieldNameDescriptorNames="DescriptorNames";
	public static final boolean compressDescriptorsInDB=false;	
	
	private static final Logger logger = LogManager.getLogger(RecordDescriptors.class);
	
	public static String [] getFields() {
		
		ArrayList<String> fields=new ArrayList<String>();
		
		fields.add(fieldNameSoftware);		
		fields.add(fieldNameKey);
		fields.add(fieldNameDescriptors);
				
		String[] strFields = new String[fields.size()];
		strFields = fields.toArray(strFields);

		return strFields;

	}
	
	private static String getDescriptors(Statement stat,String CID,String software) {
//		long t1=System.currentTimeMillis();

		try {
			ResultSet rs = getRecords(stat,CID,software);
			if (!rs.next()) return null;
			return getDescriptors(rs);
			
		} catch (Exception ex) {
			ex.printStackTrace();
			return null;
		}
	}

	private static ResultSet getRecords(Statement stat,String CID,String software) {

		try {
			String query="select * from "+tableNameDescriptors+" where "+fieldNameKey+" = \""+CID+"\" and "+fieldNameSoftware+"=\""+software+"\";";
			ResultSet rs = stat.executeQuery(query);
			//			ResultSetMetaData rsmd = rs.getMetaData();
//			System.out.println(query);

			return rs;
		} catch (Exception ex) {
			ex.printStackTrace();
			return null;
		}
	}

	
	public static void createIndex (Statement stat) {
	
		try {
			String sqlAddIndex="CREATE INDEX if not exists idx ON "+tableNameDescriptors+" ("+fieldNameKey+")";
			stat.executeUpdate(sqlAddIndex);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}


	private static int getColNum(ResultSet rs, String fieldName) throws SQLException {
		int colNum=-1;
		for (int i=1;i<=rs.getMetaData().getColumnCount();i++) {
			if (rs.getMetaData().getColumnName(i).equals(fieldName)) {
				colNum=i;
				break;
			}
		}
		return colNum;
	}

	
	public static String lookupDescriptorsInDatabase(Connection conn,String DSSTox_Structure_Id,String descriptorSoftware) {
		//		if (usePreviousDescriptors && InChiKey!=null) {

		if (DSSTox_Structure_Id==null) return null;

		try {
			long t1=System.currentTimeMillis();
			Statement statDescriptors=conn.createStatement();
			ResultSet rs = getRecords(statDescriptors,DSSTox_Structure_Id,descriptorSoftware);
			
			if (rs.isClosed()) return null;
			
			String descriptors=getDescriptors(rs);
			long t2=System.currentTimeMillis();

			if (descriptors!=null) {
				logger.info("Time to load descriptors from db:"+(t2-t1));
				return descriptors;
			}

		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return null;
	}

	private static String getDescriptors(ResultSet rs) {

		try {
			
			String strFieldName=fieldNameDescriptors;
			
			int colNum = getColNum(rs, strFieldName);

			String strDescriptors=null;
			
			if (compressDescriptorsInDB) {			
				byte [] bytes=rs.getBytes(colNum);				
//				System.out.println(bytes.length);
				strDescriptors=StringCompression.decompress(bytes, Charset.forName("UTF-8"));
			} else {
				strDescriptors=rs.getString(colNum);
			}

//			System.out.println(strJSON);
			return strDescriptors;

		} catch (Exception ex) {
			ex.printStackTrace();
			return null;
		}
	}
	
	public static String create_sql_insert(String[] fields, String table) {
		String s = "insert into "+table+" values (";

		for (int i = 1; i <= fields.length; i++) {
			s += "?";
			if (i < fields.length)
				s += ",";
		}
		s += ");";
		return s;
	}

	public void addRecord(Connection conn) {
		
		try {

			String [] fields=RecordDescriptors.getFields();
						
			String s = create_sql_insert(fields, tableNameDescriptors);			
			PreparedStatement prep= conn.prepareStatement(s);

			int i=1;

			prep.setString(i++, this.DescriptorSoftware);
			prep.setString(i++, this.DSSTox_Structure_Id);
			
			if (compressDescriptorsInDB) {
				byte [] bytesJSON=StringCompression.compress(this.Descriptors,Charset.forName("UTF-8"));
				prep.setBytes(i++, bytesJSON);
			} else {				
				prep.setString(i++, this.Descriptors);
			}
										
			prep.executeUpdate();
//			

		} catch (Exception ex) {
			System.out.println(ex.getMessage());

			//			ex.printStackTrace();
		}
		
		
	}
	
}
