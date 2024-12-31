public class NeedlemanWunsch{
    final String seq1;
    final String seq2;
    final int nrow;
    final int ncol;
    int [][] nm_matrix;
    public NeedlemanWunsch(String seq1, String seq2){
        this.seq1=seq1;
        this.seq2=seq2;
        this.nrow = seq1.length();
        this.ncol = seq2.length();
        this.nm_matrix=new int[this.nrow + 2][this.ncol + 2];
    }

    public String get_seqs(){
        return ("seq1: "+ seq1 + " "+ "seq2: "+ seq2);
        
    }

    public int[][] genMatrix(){
        for (int i = 0; i < this.nrow + 2; i++) {
            for (int j = 0; j < this.ncol + 2; j++) {
                nm_matrix[i][j] = 0;  // Example initialization, set your logic here if needed CHATGPT
        }
        }

        return nm_matrix;
        
    }
    public void showMatrix(){
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                System.out.print(nm_matrix[i][j] + " ");
        }
        
        System.out.println();
        }
    }
}
