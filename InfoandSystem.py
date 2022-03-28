#!/usr/bin/env python
from string import ascii_uppercase

def kMap(size, equation):
    for i in range(len(equation)):
        letter = equation[i]
        if letter in ascii_uppercase:
            if(equation[i+1] == '\''):
                print(False)
            else:
                print(True)
