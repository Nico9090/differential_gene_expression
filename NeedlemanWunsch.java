public class NeedlemanWunsch{
    final String seq1;
    final  String seq2;
   
    public NeedlemanWunsch(String seq1, String seq2){
        this.seq1=seq1;
        this.seq2=seq2;
    }

    public String get_seqs(){
        return ("seq1: "+ seq1 + "seq2: "+ seq2);
        
    }
    
}
