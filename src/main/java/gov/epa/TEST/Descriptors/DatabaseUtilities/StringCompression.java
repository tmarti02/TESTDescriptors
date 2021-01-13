package gov.epa.TEST.Descriptors.DatabaseUtilities;

import java.io.ByteArrayOutputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class StringCompression {
	public static byte[] compress(String str, Charset charset) {
	    ByteArrayOutputStream buffer = new ByteArrayOutputStream();
	    try {
	      OutputStream deflater = new GZIPOutputStream(buffer);
	      deflater.write(str.getBytes(charset));
	      deflater.close();
	    } catch (IOException e) {
	      throw new IllegalStateException(e);
	    }
	    
	    return buffer.toByteArray();
	    
	  }
	
	
	public static byte[] compressBytes(byte[] bytes) {
	    ByteArrayOutputStream buffer = new ByteArrayOutputStream();
	    try {
	      OutputStream deflater = new GZIPOutputStream(buffer);
	      deflater.write(bytes);
	      deflater.close();
	    } catch (IOException e) {
	      throw new IllegalStateException(e);
	    }
	    
	    return buffer.toByteArray();
	    
	  }
	
//	public static String compressToString(String str, Charset charset) {
//	    ByteArrayOutputStream buffer = new ByteArrayOutputStream();
//	    try {
//	      OutputStream deflater = new GZIPOutputStream(buffer);
//	      deflater.write(str.getBytes(charset));
//	      deflater.close();
//	    } catch (IOException e) {
//	      throw new IllegalStateException(e);
//	    }
//	    
//	    
//	    return buffer.toString();
//	    
//	  }

	  public static String decompress(byte[] data,
	      Charset charset) {
	    ByteArrayOutputStream buffer = new ByteArrayOutputStream();
	    ByteArrayInputStream in = new ByteArrayInputStream(data);
	    try {
	      InputStream inflater = new GZIPInputStream(in);
	      byte[] bbuf = new byte[256];
	      while (true) {
	        int r = inflater.read(bbuf);
	        if (r < 0) {
	          break;
	        }
	        buffer.write(bbuf, 0, r);
	      }
	    } catch (IOException e) {
	      throw new IllegalStateException(e);
	    }
	    return new String(buffer.toByteArray(), charset);
	  }

	  public static byte [] decompressToBytes(byte[] data) {
		    ByteArrayOutputStream buffer = new ByteArrayOutputStream();
		    ByteArrayInputStream in = new ByteArrayInputStream(data);
		    try {
		      InputStream inflater = new GZIPInputStream(in);
		      byte[] bbuf = new byte[256];
		      while (true) {
		        int r = inflater.read(bbuf);
		        if (r < 0) {
		          break;
		        }
		        buffer.write(bbuf, 0, r);
		      }
		    } catch (IOException e) {
		      throw new IllegalStateException(e);
		    }
		    return buffer.toByteArray();
		  }
	  
	  public static void main(String[] args) throws IOException {
	    StringBuilder sb = new StringBuilder();
	    while (sb.length() < 10000) {
	      sb.append("write the data here \u00A3");
	    }
	    String str = sb.toString();
	    Charset utf8 = Charset.forName("UTF-8");
	    byte[] compressed = compress(str, utf8);
	    
	    
	    System.out.println("String len=" + str.length());
	    System.out.println("Encoded len="
	        + str.getBytes(utf8).length);
	    System.out.println("Compressed len="
	        + compressed.length);

	    String decompressed = decompress(compressed, utf8);
	    System.out.println(decompressed.equals(str));
	    

	  }

}
