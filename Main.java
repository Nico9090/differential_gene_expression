public class Main{
    public static void main(String[] args) {
        NeedlemanWunsch a1 = new NeedlemanWunsch("GCATGCG","GATTACA");
        System.out.println(a1.get_seqs());
        a1.genMatrix();
        a1.showMatrix();
    }
}
