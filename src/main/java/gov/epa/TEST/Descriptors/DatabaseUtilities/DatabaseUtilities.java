package gov.epa.TEST.Descriptors.DatabaseUtilities;

import java.nio.charset.Charset;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

import org.openscience.cdk.AtomContainer;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;

public class DatabaseUtilities {
	
	private static final Logger logger = LogManager.getLogger(DatabaseUtilities.class);
	
	private static final boolean compressDescriptorsInDB=false;
	
	private static String [] getFieldsDescriptors() {
		
		ArrayList<String> fields=new ArrayList<String>();

		fields.add("CAS");
//		fields.add("InChi");
		fields.add("InChiKey");
		fields.add("Descriptors");
		fields.add("error");
				
		String[] strFields = new String[fields.size()];
		strFields = fields.toArray(strFields);

		return strFields;

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
	
	static void createDescriptorsDB(Connection connDescriptors,String primaryKey) throws SQLException {
//		String [] fields= {"CAS","Descriptors"};//Add inchii key???
		Statement stat=connDescriptors.createStatement();
		String [] fields=getFieldsDescriptors();
		create_table(stat,"Descriptors",fields);
	}
	
	
//	static void createDescriptorsTableMySQL(Connection connDescriptors,String primaryKey) throws SQLException {
////		String [] fields= {"CAS","Descriptors"};//Add inchii key???
//		Statement stat=connDescriptors.createStatement();
//		String [] fields=getFieldsDescriptors();
//		int [] lengths=getFieldLengthsDescriptors();
//		MySQL_DB.create_table_mysql(stat,"Descriptors",fields,lengths,primaryKey);
//	}
	
	
	public static void create_table (Statement stat,String table,String []fields) {
		
		try {
	
			String sql = "create table if not exists " + table + " (";
	
			int count = 0;// number of fields
	
	
			for (int i = 0; i < fields.length; i++) {
				sql += fields[i] + " TEXT,";
				count++;
			}
	
	
			// Trim off trailing comma:
			if (sql.substring(sql.length() - 1, sql.length()).equals(",")) {
				sql = sql.substring(0, sql.length() - 1);
			}
	
			sql += ");";
	
			//			System.out.println(sql);
	
			stat.executeUpdate(sql);
	
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	
	
	
	
	}
	
	private static int [] getFieldLengthsDescriptors() {
		
		ArrayList<Integer> fields=new ArrayList<Integer>();
		fields.add(15);//CAS
//		fields.add("InChi");
		fields.add(30);//InChiKey
		fields.add(25000);//Descriptors- should we zip it?
		fields.add(100);//error
				
		int[] intLength = new int[fields.size()];
		
		for (int i=0;i<fields.size();i++) {
			intLength[i]=fields.get(i);
		}
		
		return intLength;

	}
	
	public static DescriptorData getDescriptors(ResultSet rs,boolean compressDescriptorsInDB) {

		try {
			
			String strDescriptors="Descriptors";
			
			int colNum = getColNum(rs, strDescriptors);

			String strJSON=null;
			if (compressDescriptorsInDB) {
				byte [] bytes=rs.getBytes(colNum);
				strJSON=StringCompression.decompress(bytes, Charset.forName("UTF-8"));
			} else {
				strJSON=rs.getString(colNum);
			}

//			System.out.println(strJSON);
			
			
			if (strJSON==null || strJSON.isEmpty()) return null;
			
			DescriptorData dd=null;//TODO
			return dd;

		} catch (Exception ex) {
			ex.printStackTrace();
			return null;
		}
	}
	
	public static DescriptorData lookupDescriptorsInDatabase(AtomContainer ac,Connection conn,String fieldValue,String fieldName) {
		//		if (usePreviousDescriptors && InChiKey!=null) {

		if (fieldValue==null) return null;

		try {
			long t1=System.currentTimeMillis();
			Statement statDescriptors=conn.createStatement();
			DescriptorData dd=DatabaseUtilities.getDescriptors(statDescriptors, fieldName,fieldValue,compressDescriptorsInDB);
			long t2=System.currentTimeMillis();

			if (dd!=null) {
				logger.info("Time to load descriptors from db:"+(t2-t1));
				return dd;
			}

		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return null;
	}
	
	public static DescriptorData getDescriptors(Statement stat,String searchField,String searchValue,boolean compressDescriptorsInDB) {
//		long t1=System.currentTimeMillis();

		try {
			
			ResultSet rs = getRecords(stat,"Descriptors", searchField,searchValue);
			
			if (!rs.next()) return null;
			
			return getDescriptors(rs,compressDescriptorsInDB);
			
		} catch (Exception ex) {
			ex.printStackTrace();
			return null;
		}

	}
	
	public static ResultSet getRecords(Statement stat,String tableName,String keyField,String keyValue) {

		try {
			String query="select * from "+tableName+" where "+keyField+" = \""+keyValue+"\";";
			//			System.out.println(query);
			ResultSet rs = stat.executeQuery(query);
			//			ResultSetMetaData rsmd = rs.getMetaData();

			return rs;
			//			this.printResultSet(rs);
		} catch (Exception ex) {
			ex.printStackTrace();
			return null;
		}
	}
	
	public static void addRecordsToDescriptorDatabase(Connection conn,DescriptorData dd,boolean compressDescriptors) {
		try {

			String [] fields=getFieldsDescriptors();
			String table="Descriptors";
			
			String s = create_sql_insert(fields, table);			
			PreparedStatement prep= conn.prepareStatement(s);

			int i=1;
			prep.setString(i++, dd.ID);
//			prep.setString(i++, dd.InChi);
			prep.setString(i++, dd.InChiKey);
			
			//			prep.setString(i++, dd.to_JSON_String(dd,true));
//			long t1=System.currentTimeMillis();
			String strDescriptors="";
			
//			System.out.println(dd.CAS+"\t"+strJSON);
			
//			long t2=System.currentTimeMillis();
//			logger.info("Time to convert to json string="+(t2-t1));
			
			
			if (compressDescriptors) {
				byte [] bytesJSON=StringCompression.compress(strDescriptors,Charset.forName("UTF-8"));
				prep.setBytes(i++, bytesJSON);
			} else {				
				prep.setString(i++, strDescriptors);
			}
			
			prep.setString(i++, dd.Error);
			
//			long t3=System.currentTimeMillis();
//			logger.info("Time to compress="+(t3-t2));
			
//			prep.addBatch();
//			int [] count=prep.executeBatch();
			
			prep.executeUpdate();
//			
//			long t4=System.currentTimeMillis();
//			logger.info("Time to execute update="+(t4-t3));
			
			//			System.out.println(count.length);

		} catch (Exception ex) {
			System.out.println(ex.getMessage());

			//			ex.printStackTrace();
		}
	}

	private static String create_sql_insert(String[] fields, String table) {
		String s = "insert into "+table+" values (";

		for (int i = 1; i <= fields.length; i++) {
			s += "?";
			if (i < fields.length)
				s += ",";
		}
		s += ");";
		return s;
	}
	
	
		
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

