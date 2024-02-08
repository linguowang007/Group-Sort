import bisect
import os
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
from share_func import pairwise
from numba_type_counter import TypeCounter


class RegionMaker:
    def __init__(self, *, sam_file, bin_list, nt=16):
        self.sam_file = sam_file
        self.bin_list = bin_list  # bin-file-chrome-id, chrome-id-number
        self.nt = nt
    
    @staticmethod
    def get_break_pos(*, num_acc, chunk_size):
        current_val = 0
        sam_size = num_acc[-1]
        while current_val < sam_size:
            current_val += chunk_size
            index = bisect.bisect_left(num_acc, current_val)
            if index >= len(num_acc):
                yield len(num_acc)
                return
            yield index
            current_val = num_acc[index]
    
    @staticmethod
    def merge_index(out_name, chr_num, bin_files):
        c1 = TypeCounter()
        for file in bin_files:
            data1 = np.fromfile(file, dtype=np.uint32)
            c1.from_one_numpy_array(data1)
        pos, len1 = c1.to_sort_pos_len()
        del c1
        pos.tofile(f'{out_name}-pos')
        del pos
        len1.tofile(f'{out_name}-len')
        size = len(len1)
        del len1
        chr_data = np.full(size, chr_num, dtype=np.uint32)
        chr_data.tofile(f'{out_name}-chr')
        del chr_data
        return chr_num, f'{out_name}-chr', f'{out_name}-pos', f'{out_name}-len'
    
    def make_regions(self):
        d1 = defaultdict(list)
        for row in self.bin_list:
            d1[row[-1]].append(row[0])  # key as the chr-number, val as the bin-file
        args = [(f'{self.sam_file}-{key}.bin', key, vals) for key, vals in d1.items()]
        with Pool(self.nt) as p:
            result = p.starmap(RegionMaker.merge_index, args, chunksize=1)
            result = sorted(result)
            chr_file = [i[1] for i in result]
            pos_file = [i[2] for i in result]
            len_file = [i[3] for i in result]
            cmd1 = rf'''cat {' '.join(chr_file)} > {self.sam_file}-chr-sort.bin && rm {' '.join(chr_file)}'''
            cmd2 = rf'''cat {' '.join(pos_file)} > {self.sam_file}-pos-sort.bin && rm {' '.join(pos_file)}'''
            cmd3 = rf'''cat {' '.join(len_file)} > {self.sam_file}-len-sort.bin && rm {' '.join(len_file)}'''
            p.map(os.system, [cmd1, cmd2, cmd3], chunksize=1)

        len1 = np.fromfile(f'{self.sam_file}-len-sort.bin', dtype=np.uint32).astype(np.uint64)
        len1 = np.add.accumulate(len1)
        chunk_size = int(len1[-1] / self.nt) + 1
        break_pos = [0] + list(RegionMaker.get_break_pos(num_acc=len1, chunk_size=chunk_size))
        break_pos = [(a1, a2) for a1, a2 in pairwise(break_pos)]
        del len1

        chrome = np.fromfile(f'{self.sam_file}-chr-sort.bin', dtype=np.uint32)
        pos = np.fromfile(f'{self.sam_file}-pos-sort.bin', dtype=np.uint32)
        uniq_val = [(num, (a, b), np.unique(chrome[a:b])) for num, (a, b) in enumerate(break_pos)]
        all_pos = []
        group_order = [f'group{num}' for num, *_ in enumerate(break_pos)]
        for num, (a, b), uniq_num in uniq_val:
            region_pos = pos[a:b]
            for i in uniq_num:
                b_pos = bisect.bisect_left(chrome[a:b], i)
                e_pos = bisect.bisect_right(chrome[a:b], i)
                head = region_pos[b_pos]
                bot = region_pos[e_pos - 1]
                all_pos.append([f'group{num}', i, head, bot])
        return all_pos, group_order
