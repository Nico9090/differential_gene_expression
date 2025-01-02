public class NeedlemanWunsch{
    final String seq1;
    final String seq2;
    final int nrow;
    final int ncol;
    int [][] nm_matrix;
    int match;
    int mismatch;
    public NeedlemanWunsch(String seq1, String seq2){
        this.seq1=seq1;
        this.seq2=seq2;
        this.nrow = seq1.length()+2;
        this.ncol = seq2.length()+2;
        this.nm_matrix=new int[this.nrow][this.ncol];
        this.match=1;
        this.mismatch=-1;
    }

    public String get_seqs(){
        return ("seq1: "+ seq1 + " "+ "seq2: "+ seq2);
        
    }

    public int[][] genMatrix(){
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                nm_matrix[i][j] = 0;  // Example initialization, set your logic here if needed CHATGPT
        }
        }

        return nm_matrix;
        
    }
    public void showMatrix(){
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                System.out.printf("%5d ",nm_matrix[i][j] );
        }
        
        System.out.println();
        }
    }
    public void doAlgorithm(){
        //Column values -1,-2,...
        for (int num=1,val=0;num< this.nrow;num++,val--){
            nm_matrix[1][num]=val;
        }
        //Row values -1,-2,...
        for (int num=1,val=0;num< this.nrow;num++,val--){
            nm_matrix[num][1]=val;
        }
    }
    public void update(){ //doesn't work yet
        for(int k=0;k<seq1.length();k++){
            for(int i=0,j=1,r=2;i<this.nrow-1 && j<this.nrow-1;i++,j++,r++){
                if(seq1.charAt(k)==seq2.charAt(i)){
                    int a =nm_matrix[1][j]+this.match;
                    int b =nm_matrix[1][j+1]+this.match;
                    int c= nm_matrix[j+1][i]+this.match;
                    nm_matrix[2][r]=Math.max(Math.max(a,b),c);
            }
        }
        }
        
    }
}
