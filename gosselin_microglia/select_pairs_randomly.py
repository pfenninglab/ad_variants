from random import sample
import random
from itertools import chain
random.seed(17)

if __name__=="__main__":
    n = 12
    times_selected = dict()
    selected_pairs = []
    for i in range(n):
        times_selected[i] = 0
    number_range = list(range(n))
    done = False
    while not done:
        sampled_numbers = sorted(sample(number_range, 2))
        if not sampled_numbers in selected_pairs:
            if times_selected[sampled_numbers[0]]<2 and times_selected[sampled_numbers[1]]<2:
                times_selected[sampled_numbers[0]] = times_selected[sampled_numbers[0]] + 1
                times_selected[sampled_numbers[1]] = times_selected[sampled_numbers[1]] + 1        
                selected_pairs.append(sampled_numbers)
        
        done = True
        for i in range(n):
            if times_selected[i] < 2:
                done = False
    
    print(selected_pairs)
