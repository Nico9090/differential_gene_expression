//NeedlemanWunsch.java
import java.util.Arrays;
import java.util.Scanner;

public class NeedlemanWunsch{
// Creating variables that will store the scores for when you have a match, mismatch, or gap
        private int match;
        private int mismatch;
        private int gap;
        private String seq1;
        private String seq2;
//The scores are all whole numbers
        public NeedlemanWunsch(int match, int mismatch, int gap, String seq1, String seq2){
                this.match = match;
                this.mismatch = mismatch;
                this.gap = gap;
                this.seq1 = seq1;
                this.seq2 =  seq2;
        }
// Creating methods so you are able to set and view your own scores for a match, mismatch, and gap
        public int getMatch(){
                return match;
        }
        public int getMismatch(){
                return mismatch;
        }
        public int getGap(){
                return gap;
        }

        public String getseq1(){
                return seq1;
        }

        public String getseq2(){
                return seq2;
        }

        public void setMatch(int match){
                this.match = match;
        }
        public void setMismatch(int mismatch){
                this.mismatch = mismatch;
        }
        public void setGap(int gap){
                this.gap = gap;
        }

        public void setSeq1(String seq1){
                this.seq1 = seq1;
        }

        public void setSeq2(String seq2){
                this.seq2 = seq2;
        }
//Creating a method to store provided sequences from the user
// Creating the matrix that will be used for the scoring system
       public int [][] ar1 = new int [seq1.length()][seq2.length()];
}
