#!/usr/bin/env python3
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

print(left_circular_shift(tString))
