#!/usr/bin/env python3
#from collections import Counter
options= '''red orange yellow green blue indigo violet'''.split(' ')

chosen_colors=[]
for color in options:
        for x in range(len(options)-(1+options.index(color))):
                if color!=options[x+1]:
                        choice= input(f'{color} or {options[options.index(color)+x+1]}: ')
                        if choice != color and choice !=options[x+1]:
                                print('Error! Select one of the above options')
                                continue
        #chosen_colors.append(choice)

#print(chosen_colors.count('red'))
