package gov.epa.TEST.Descriptors.DescriptorUtilities;

public class TESTConstants {
	public static String SoftwareVersion = "5.1";
	public static String SoftwareTitle = "T.E.S.T (Toxicity Estimation Software Tool)";		
	public static final String timeoutVariableName = "path.timeout";
	public static final String defaultPathGenerationTimeout = "120000"; // 2 mins
	public static int pathGenerationTimeout = 
	        Integer.valueOf(System.getProperty(timeoutVariableName, 
            defaultPathGenerationTimeout));

}
