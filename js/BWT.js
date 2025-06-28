//BWT constructor
const sq = "GACTATATCCTAAATACCCGCACCATTACCGACACCCGTGGCCCAAGCAG";//"let" allows variable reassignment
let mod_sq = sq + "$";
console.log(mod_sq);

function left_CS(sequence){
    let all_leftCS = []; //list
    for (let i = 0; i < mod_sq.length; i++) {
        all_leftCS.push(sequence.slice(1)+sequence[0]);
        sequence=sequence.slice(1)+sequence[0];
}
    return all_leftCS;
}

let lc_shifts=left_CS(mod_sq);
let srt=lc_shifts.sort();
let bwt=srt.map(chr => chr[chr.length -1]); //.map() for iteration

let df=[
    {"CS":lc_shifts},
    {"Sorted CS":srt},
    {"BWT":bwt}
]
console.log(df)