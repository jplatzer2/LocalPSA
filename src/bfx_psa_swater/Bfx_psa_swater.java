package bfx_psa_swater;

/**
 *
 * @author Joseph Platzer
 */
public class Bfx_psa_swater {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String q = "GCTGGAAGGCAT";
        String p = "GCAGAGCACG";
        System.out.println("Locally aligning the following strings:");
        System.out.println("String q: " + q);
        System.out.println("String p: " + p);
        
        //constructor builds matrix with passed strings
        LocalAligner aligner = new LocalAligner(q, p);
        
        int maxAlignmentScore = aligner.getMaxScore();
        System.out.println("Getting the max alignment score: " + maxAlignmentScore);
        System.out.println("This is the alignment corresponding to that score: \n");
        aligner.printMaxAlignment(maxAlignmentScore);
        System.out.println();
        
        System.out.println("Printing out the Scoring Matrix for the alignment:\n");
        aligner.printScoreMatrix();
        //print local alignment
    }
    
}
