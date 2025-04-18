#!/usr/bin/env python3
import pandas as pd
#BW transform__________________________________________
class BWT:
    def __init__(self,string):
        self.string=string + "$"
    def left_circular_shift(self):
        all_left_circular_shifts=[]
        x=0
        temp_string=self.string #use temp variable to avoid permanently changing self.string
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
seq="GACTATATCCTAAATACCCGCACCATTACCGACACCCGTGGCCCAAGCAG"

data=BWT(seq)
print(data.BWT_data_frame())
