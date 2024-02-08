import itertools
import bisect
import os
from collections import Counter
from multiprocessing import Pool
from share_func import get_pos_len


class ChunkRegions:
    def __init__(self, *, sam_file, bin_list, region_list, group_order, nt):
        self.sam_file = sam_file
        self.bin_list = bin_list  # chunk-file, bin-file, chr-id
        self.nt = nt
        self.region_list = region_list
        self.group_order = group_order
    
    def group_region(self):
        bin_dict = {(a, c): b for a, b, c in self.bin_list}
        chunk_files = {i[0] for i in self.bin_list}  # keep unique chunk-file
        args = [(bin_dict, chunk_file, self.region_list, self.group_order) for chunk_file in chunk_files]
        with Pool(self.nt) as p:
            result = p.starmap(ChunkRegions.count_group_len, args, chunksize=1)
            result = list(itertools.chain.from_iterable(result))
        return result
    
    @staticmethod
    def count_group_len(bin_dict, chunk_file, region_list, group_oder):
        keep_sum = Counter()
        for row in region_list:
            g_id, chrome_id, start, end = row
            key = (chunk_file, chrome_id)
            if key in bin_dict:
                bin_file = bin_dict[key]
                d1 = get_pos_len(bin_file)
                sub_pos = d1[::2]
                sub_len = d1[1::2]
                want_start = bisect.bisect_left(sub_pos, start)
                want_end = bisect.bisect_right(sub_pos, end)
                keep_sum[g_id] += int(sub_len[want_start:want_end].sum())  # each regional sum length
        tmp_result = [(chunk_file, i, keep_sum[i]) for i in group_oder if i in keep_sum]  # re-order the group count
        tmp_acc = [0] + list(itertools.accumulate(i[-1] for i in tmp_result))
        tmp_result = [list(i) + [acc] for i, acc in zip(tmp_result, tmp_acc)]
        rm_bin = [val for key, val in bin_dict.items() if chunk_file in key]
        os.system(rf'''rm {' '.join(rm_bin)}''')  # rm bin-files associated with the chunk-file
        return tmp_result
