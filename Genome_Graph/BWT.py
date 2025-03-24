#!/usr/bin/env python3
import pandas as pd
#BW transform
seq="GACTATATCCTAAATACCCGCACCATTACCGACACCCGTGGCCCAAGCAG"
tString=seq + "$"
print(tString)

def left_circular_shift(string):
        all_left_circular_shifts=[]
        x=0
        for i in range(len(string)):
                all_left_circular_shifts.append(string[1:]+string[0])
                string=string[1:]+string[0]
                x+=1
        return all_left_circular_shifts

lc_shifts=left_circular_shift(tString)

sorted_lc_shifts=sorted(lc_shifts)
bwt = ''.join([shift[-1] for shift in sorted_lc_shifts])
#BWT table
data={
    'Circular Shifts': lc_shifts,
    'Sorted Circular Shifts': sorted_lc_shifts,
    'BWT': [bwt]*len(lc_shifts)
}
df=pd.DataFrame(data)
print(df)
