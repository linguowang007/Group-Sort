import numpy as np
from numba import typed, uint32, types
from numba.experimental import jitclass

spec = [('counter', types.DictType(uint32, uint32))]


@jitclass(spec)
class TypeCounter:
    def __init__(self):
        self.counter = typed.Dict.empty(key_type=uint32, value_type=uint32)
    
    def from_two_numpy_array(self, num1, num2):
        for a, b in zip(num1, num2):
            if a in self.counter:
                self.counter[a] = np.add(self.counter[a], b)
            else:
                self.counter[a] = b
    
    def from_one_numpy_array(self, data1):
        num1 = data1[::2]
        num2 = data1[1::2]
        self.from_two_numpy_array(num1, num2)
    
    def to_sort_pos_len(self):
        data1 = np.empty(self.getsize(), dtype=np.uint32)
        data2 = np.empty(self.getsize(), dtype=np.uint32)
        i = 0
        for a, b in self.counter.items():
            data1[i] = a
            data2[i] = b
            i += 1
        keys = np.argsort(data1)
        data1 = data1[keys]
        data2 = data2[keys]
        return data1, data2
    
    def getsize(self):
        return len(self.counter)
    
    @staticmethod
    def add_chr_and_sort(data1, chr_id):
        num1 = data1[::2]
        num2 = data1[1::2]
        in_key = np.argsort(num1)
        num1 = num1[in_key]
        num2 = num2[in_key]
        
        num3 = np.full(len(num1), chr_id, dtype=np.uint32)
        big_data = np.empty(len(num1) * 3, dtype=np.uint32)
        i = 0
        for chrome, a, b in zip(num3, num1, num2):
            big_data[i] = chrome
            big_data[i + 1] = a
            big_data[i + 2] = b
            i += 3
        return big_data
