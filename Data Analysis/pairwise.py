#!/usr/bin/env python3
#from collections import Counter
options= '''red orange yellow green blue indigo violet'''.split(' ')

chosen_colors=[]

i=0
j=-1

while options[i]!=options[j] and i<len(options):
        choice=input(f'{options[i]} or {options[j]}: ')
        if choice!= options[i] and choice !=options[j]:
                print('Error! Select one of the options')
                choice=input(f'{options[i]} or {options[j]}: ')
        if choice == options[i] or choice == options[j]:
                chosen_colors.append(choice)
        j-=1
        if options[j]==options[i]:
                j=-1
                i+=1

for color in options:
        print(color,chosen_colors.count(color),sep=':')


