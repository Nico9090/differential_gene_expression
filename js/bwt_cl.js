export class bwt_cl{//must export because each script is run in individual scopes
    constructor(sq,tsq){
        this.sq=sq;
        this.tsq=sq+"$"


    }
    cs_lt(){
        let all_leftCS = []; //list
        let sequence=this.tsq;
        for (let i = 0; i < sequence.length; i++) {
            all_leftCS.push(sequence.slice(1)+sequence[0]);
            sequence=sequence.slice(1)+sequence[0];
}
    return all_leftCS;
    }

    table(){
        let lc_shifts=this.cs_lt();//must call method with this.
        let srt=lc_shifts.sort();
        let bwt=srt.map(chr => chr[chr.length -1]); //.map() for iteration

        let df={
            "CS":lc_shifts,
            "Sorted CS":srt,
            "BWT":bwt
        };
        return df;
    }
}