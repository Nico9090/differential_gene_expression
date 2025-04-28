#!/usr/bin/env python3
import pandas as pd
from collections import defaultdict
#BW transform__________________________________________
class BWT:
    def __init__(self,string):
        self.string=string + "$"
    def left_circular_shift(self):
        all_left_circular_shifts=[]
        x=0
        temp_string=self.string
        for i in range(len(self.string)):
            all_left_circular_shifts.append(temp_string[1:]+temp_string[0])
            temp_string=temp_string[1:]+temp_string[0]
            x+=1
        return all_left_circular_shifts

    def BWT_data_frame(self):
        lc_shifts=self.left_circular_shift()
        sorted_lc_shifts=sorted(lc_shifts)
        bwt = ''.join([shift[-1] for shift in sorted_lc_shifts])
        #BWT table
        data={
            'Circular Shifts': lc_shifts,
            'Sorted Circular Shifts': sorted_lc_shifts,
            'BWT': [bwt]*len(lc_shifts)
        }
        df=pd.DataFrame(data)
        return df

    def c_array(self):
        df = self.BWT_data_frame()
        chars={char:0 for char in self.string if char != "$"}
        for letter in self.string:
            if letter != "$":
                chars[letter] += 1
        chars=sorted(chars.keys())
        array_c={char:0 for char in chars}
        second_of_chars=list(chars.keys())[1]
        second_char=list(array_c.keys())[1]
        array_c[second_char]+=second_of_chars
        return chars

seq="GACTATATCCTAAATACCCGCACCATTACCGACACCCGTGGCCCAAGCAG"

data=BWT(seq)
print(data.c_array())


