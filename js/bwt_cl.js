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
    cARR(){
        const dataframe = this.table();
        let init_list=[];
        let key=[];
        for (a_char in dataframe["BWT"]){
            init_list.push(a_char)
        };
        let ky= new Set(init_list);
        for (a_char in ky){
            key.push(a_char)
        };
        let dict_key={};
        let k=[key];
        for (let v of k){
            dict_key[v]=0;
        };
        return dict_key;
            
    }
}