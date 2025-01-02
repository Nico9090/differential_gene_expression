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
        nm_matrix[1][1]=0;
        nm_matrix[1][2]=0;
        for (int num=1,mismatch=-1;num< this.nrow;num++,mismatch--){
            nm_matrix[1][num]=mismatch;
        }

        for (int num=1,match=1;num< this.nrow;num++,match--){
            nm_matrix[num][1]=match;
        }
        for (int start=0,num=2; start < this.nrow && num <this.nrow;start++, num++){
            if (seq1.charAt(start)==seq2.charAt(start)){
            nm_matrix[num][num]++;
        }
        

        }
    }
    public void update(){
        if(seq1.charAt(0)==seq2.charAt(0)){
            int a=nm_matrix[1][1]+this.match;
            int b=nm_matrix[1][2]+this.match;
            int c=nm_matrix[2][1]+this.match;
            nm_matrix[2][2]=Math.max(Math.max(a,b),c);
            System.out.println("a : "+  a);
            System.out.println("b : "+  b);
            System.out.println("c : "+  c);
            System.out.println("nm_matrix[1][1] : "+  nm_matrix[1][1] );
            System.out.println("match : "+  this.match );
        }
    }
}
