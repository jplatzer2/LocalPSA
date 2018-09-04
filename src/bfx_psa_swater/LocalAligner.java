package bfx_psa_swater;
import java.util.Arrays;
import java.util.Collections;

/**
 *
 * @author Joseph Platzer
 */
public class LocalAligner {
    private String q;
    private String p;
    private int qlength;
    private int plength;
    private int[][] scoreMatrix;
    private char[][] pointerMatrix;
    private static final int MATCH_SCORE = 5;
    private static final int MISMATCH_SCORE = -4;
    private static final int GAP_SCORE = -4;
    //below indicate the direction to best neighbor for tracing local alignments
    private static final char DIR_DIAG = 'D';
    private static final char DIR_UP = 'U';
    private static final char DIR_LEFT = 'L';
    private static final char DIR_ZERO = 'Z';
    
    public LocalAligner(){
    }

    public LocalAligner(String q, String p) {
        this.q = q;
        this.p = p;
        this.qlength = q.length();
        this.plength = p.length();
        this.scoreMatrix = new int[qlength+1][plength+1];
        this.pointerMatrix = new char[qlength + 1][plength + 1];
        
        buildAlignmentMatrix();
    }
   
    /**
     * @return the scoreMatrix
     */
    public int[][] getScoreMatrix() {
        return scoreMatrix;
    }

    /**
     * @param scoreMatrix the scoreMatrix to set
     */
    public void setScoreMatrix(int[][] scoreMatrix) {
        this.scoreMatrix = scoreMatrix;
    }
    
    /** 
     * Build the scoring matrix with the size of the given strings.
     */
    private void buildAlignmentMatrix(){
        
        //first task is to set up matrix - so set 0's along first rows and column
        
        //row
        for(int i = 0; i <= qlength; i++){
            scoreMatrix[i][0] = 0;
            pointerMatrix[i][0] = DIR_ZERO;
        }
        //column
        for(int j = 0; j <= plength; j++){
            scoreMatrix[0][j] = 0;
            pointerMatrix[0][j] = DIR_ZERO;
        }
        
        //fill in the body of the matrix
        for(int i = 1; i <= qlength; i++){
            for(int j = 1; j <= plength; j++){
                int diagonalScore = scoreMatrix[i-1][j-1] + checkIfMatch(i,j);
                int upScore = scoreMatrix[i][j-1] + checkIfMatch(0,j);
                int leftScore = scoreMatrix[i-1][j] + checkIfMatch(i,0);
                
                //assign score to cell
                int scores[] = {0, diagonalScore, upScore, leftScore};
                int max = getMax(scores);
                scoreMatrix[i][j] = max;
                pointerMatrix[i][j] = DIR_ZERO;
                
                //assign pointer to cell
                if (diagonalScore == scoreMatrix[i][j]) {
                    pointerMatrix[i][j] = DIR_DIAG;
                }
                if (leftScore == scoreMatrix[i][j]) {
                    pointerMatrix[i][j] = DIR_LEFT;
                }
                if (upScore == scoreMatrix[i][j]) {
                    pointerMatrix[i][j] = DIR_UP;
                }
                if (0 == scoreMatrix[i][j]) {
                    pointerMatrix[i][j] = DIR_ZERO;
                }                
            }
        }
    }
    
    /**
     * Traverses the scoreMatrix to find top score.
     * @return max score
     */
    public int getMaxScore(){
        int max = 0;
        
        //loop through scoreMatrix to find the top score
        for (int i = 1; i <= qlength; i++) {
            for (int j = 1; j <= plength; j++) {
                if (scoreMatrix[i][j] > max) {
                    max = scoreMatrix[i][j];
                }
            }
        }
                
        return max;
    }
    
    /**
     * Print the scoring matrix.
     */
    public void printScoreMatrix(){
    
        //print p sequence row title
        System.out.print("   "); //gap
        for (int j=1; j <= plength; j++){
            System.out.print("   " + p.charAt(j-1));
        }

        //print out q seq and scores
        System.out.println();
        for (int i=0; i <= qlength; i++){
            if (i > 0){
                System.out.print(q.charAt(i-1) + " ");
            } else {
                System.out.print("  ");
            }

            for (int j=0; j <= plength; j++){
                System.out.print(scoreMatrix[i][j] + "   ");
            }
            System.out.println(); //essential to form a grid, not a long line...
        }
    }    
    
    /**
     * Method to check if bases in each sequence are matching.
     * @param i
     * @param j
     * @return Returns gap, match, or mismatch score depending on char at position in each string.
     */
    private int checkIfMatch(int i, int j){
        //this activates if there is a gap
        if (i == 0 || j == 0){
            return GAP_SCORE;
        }
        char qBase = q.charAt(i-1);
        char pBase = p.charAt(j-1);
        if (qBase == pBase){
            return MATCH_SCORE;
        } else {
            return MISMATCH_SCORE;
        }
        
    }
    
    /**
     * Helper function to get max of smith-waterman scoring algorithm
     * @param array
     * @return max value of scoring, always >=0
     */
    public int getMax(int[] array){ 
    int maxValue = array[0]; 
    for(int i=1;i < array.length;i++){ 
      if(array[i] > maxValue){ 
         maxValue = array[i]; 
      } 
    } 
    return maxValue; 
  }    
    
    /**
     * Moves to the max align score, and initializes printing of aligning subsequences.
     * @param maxAlignmentScore 
     */
    public void printMaxAlignment(int maxAlignmentScore){  
        //traverse to max score, and begin printing the alignment based on score and pointer matrices
        for (int i = 1; i <= qlength; i++) {
            for (int j = 1; j <= plength; j++) {
                if (scoreMatrix[i][j] == maxAlignmentScore) {
                    grabAlignments(i, j, "", "");
                }
            }
        }
    }
    
    /**
     * Method responsible for printing out sequence alignment.
     * @param i Position in matrix
     * @param j Position in matrix
     * @param align1 Recursively passed to build alignment string
     * @param align2 Recursively passed to build alignment string
     */
    private void grabAlignments(int i, int j, String align1, String align2){
        //recursive function
        //base case
        if (pointerMatrix[i][j] == DIR_ZERO){
            System.out.println(align1);
            System.out.println(align2);
            return;
        }
        
        //backtrack through matrix with pointers
        if (pointerMatrix[i][j] == DIR_LEFT) {
            grabAlignments(i-1, j, q.charAt(i-1) + align1, "_" + align2);
        }
        if (pointerMatrix[i][j] == DIR_UP) {
            grabAlignments(i, j-1, "_" + align1, p.charAt(j-1) + align2);
        }
        if (pointerMatrix[i][j] == DIR_DIAG) {
            grabAlignments(i-1, j-1, q.charAt(i-1) + align1, p.charAt(j-1) + align2);
        }
    }
}
