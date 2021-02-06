

import static org.junit.Assert.assertEquals;

import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.smiles.SmilesParser;

import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.TEST.Descriptors.DescriptorFactory.DescriptorFactory;
import gov.epa.TEST.Descriptors.DescriptorUtilities.HueckelAromaticityDetector;

public class InformationContentDescriptorTest2 {

    SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    DescriptorFactory df = new DescriptorFactory(false);

    @Test
    public void testCalculateTti() throws InvalidSmilesException {

    	
        String smiles = "CCOS(=O)(=O)OCC";
        DescriptorData dd = calculateDescriptors(smiles);
        assertEquals(38.783562720275d, dd.tti, 1e-10d); 
        assertEquals(57.05122783315792d, dd.ttvi, 1e-10d); 

        smiles = "c1ccccc1";
        dd = calculateDescriptors(smiles);
        assertEquals(29.4d, dd.tti, 1e-10d); 
        assertEquals(44.1d, dd.ttvi, 1e-10d); 

        smiles = "O=n1ccccc1";// pyridine-N-oxide, neutral representation
        dd = calculateDescriptors(smiles);
        assertEquals(37.6082425418124d, dd.tti, 1e-10d); 
        assertEquals(66.0477795687537d, dd.ttvi, 1e-10d); 

        smiles="[O-][n+]1ccccc1";//pyridine-N-oxide, charge-separated
        dd = calculateDescriptors(smiles);
        assertEquals(37.6082425418124d, dd.tti, 1e-10d); 
        assertEquals(66.0477795687537d, dd.ttvi, 1e-10d); 

        smiles="O=N(=O)c(cccc1)c1"; // nitrobz
        dd = calculateDescriptors(smiles);
        assertEquals(54.31285352450502d, dd.tti, 1e-10d); 
        assertEquals(98.33409560909004d, dd.ttvi, 1e-10d); 

        smiles="O=[N+]([O-])c1ccccc1";
        dd = calculateDescriptors(smiles);
        assertEquals(54.31285352450502d, dd.tti, 1e-10d); 
        assertEquals(98.33409560909004d, dd.ttvi, 1e-10d); 

        smiles="OCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO";
        dd = calculateDescriptors(smiles);
        assertEquals(284.55394028154313d, dd.tti, 1e-10d); 
        assertEquals(296.76324261863573d, dd.ttvi, 1e-10d); 

        smiles="c1ccncc1";//pyridine
        dd = calculateDescriptors(smiles);
        assertEquals(29.4d, dd.tti, 1e-10d); 
        assertEquals(48.52783444620586d, dd.ttvi, 1e-10d); 

        smiles="c1ncncc1";
        dd = calculateDescriptors(smiles);
        assertEquals(29.4d, dd.tti, 1e-10d); 
        assertEquals(53.06942258740338d, dd.ttvi, 1e-10d); 

        smiles="CCCCCCCCCCCCc2cc5c1c(ccc3c1c2c1ccc4c2c1c3ccc2C(=O)c1ccccc41)c1ccccc1C5=O";
        dd = calculateDescriptors(smiles);
        assertEquals(17429.401862886916d, dd.tti, 1e-10d); 
        assertEquals(23726.685710471083d, dd.ttvi, 1e-10d); 

        smiles="COc7ccc(N4C(=O)c5ccc6c1ccc3c2c1c(ccc2C(=O)N(c1ccc(OC)cc1)C3=O)c1ccc(C4=O)c5c61)cc7";
        dd = calculateDescriptors(smiles);
        assertEquals(13944.381665673782d, dd.tti, 1e-10d); 
        assertEquals(19953.66627897736d, dd.ttvi, 1e-10d); 
    }
    
    private DescriptorData calculateDescriptors(String smiles) throws InvalidSmilesException {
        DescriptorFactory df = new DescriptorFactory(false);
        AtomContainer m = (AtomContainer) sp.parseSmiles(smiles);
        
        DescriptorData dd = new DescriptorData();
        
        dd.ThreeD = false;
        
        df.Calculate3DDescriptors = false;
        df.CalculateDescriptors(m, dd, true);
//        dd.WriteOut2DDescriptors(m, dd, "test-results", "|");
                
        
        return dd;
    }

}
