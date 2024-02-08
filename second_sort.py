import os
import gzip
import numpy as np
import mmap
import bisect
import itertools
from collections import defaultdict
from multiprocessing import Pool


class SecondSort:
    def __init__(self, *, sam_file, region_list, group_order, nt=16):
        self.region_list = region_list
        self.sam_file = sam_file
        self.nt = nt
        self.group_order = group_order
    
    @staticmethod
    def to_list_region(region_list):
        chunk_file = region_list[0][0]
        index = f'{chunk_file}.gzi'
        index = np.fromfile(index, dtype=np.uint64)[1:]
        
        first = np.array([0, 0], dtype=np.uint64)
        file_size = os.path.getsize(chunk_file)
        with os.popen(f'bgzip -b {index[-1]} {chunk_file} | wc -c') as h1:
            tmp = int([i.strip() for i in h1][0]) + index[-1]
            tmp = int(tmp)
            last = np.array([file_size, tmp], dtype=np.uint64)
        index = np.concatenate([first, index, last])
        compressed = index[::2]
        txt = index[1::2]
        
        keep_chunk = []
        with open(chunk_file, 'rb') as f_in:
            mmapped_in = mmap.mmap(f_in.fileno(), 0, access=mmap.ACCESS_READ)
            for row in region_list:
                chunk, group_id, query_len, start_pos = row
                left_index = bisect.bisect_left(txt, start_pos)
                right_index = bisect.bisect_right(txt, start_pos + query_len)
                out_name = f'{chunk}-{group_id}-chunk.gz'
                keep_chunk.append([group_id, out_name])
                with open(out_name, 'wb') as f_out:
                    if txt[left_index] != start_pos:
                        start = compressed[left_index - 1]
                        end = compressed[left_index]
                        contain1 = gzip.decompress(mmapped_in[start:end])[int(start_pos - txt[left_index]):]
                        contain1 = gzip.compress(contain1)
                        f_out.write(contain1)
                    
                    pos1 = compressed[left_index]
                    pos2 = compressed[right_index - 1]
                    f_out.write(mmapped_in[pos1:pos2])
                    
                    start = compressed[right_index - 1]
                    end = compressed[right_index]
                    contain1 = gzip.decompress(mmapped_in[start:end])[
                               :int(start_pos + query_len - txt[right_index - 1])]
                    contain1 = gzip.compress(contain1)
                    f_out.write(contain1)
            mmapped_in.close()
        return keep_chunk
    
    def second_sort(self, *, unmapped=None):
        gz_header = f'{self.sam_file}-header.txt.gz'
        os.system(f'samtools view -H {self.sam_file} | bgzip > {gz_header}')
        
        d1 = defaultdict(list)
        for row in self.region_list:
            chunk, group_id, query_length, start_pos = row
            d1[chunk].append(row)
        
        nt = self.nt
        if self.nt > 48:  # keep bgzip slicing using 48 or lower thread number
            nt = 48
        with Pool(nt) as p:
            result = p.map(SecondSort.to_list_region, d1.values(), chunksize=1)
            result = list(itertools.chain.from_iterable(result))
            d1 = defaultdict(list)  # key is group-id, val is split sam body file
        
        with Pool(self.nt) as p:
            for items in result:
                d1[items[0]].append(items[1])
            cmds = [rf'''cat {gz_header} {' '.join(vals)} | zcat | \
                samtools sort -T {self.sam_file}-{key}-merge.bam- -o {self.sam_file}-{key}-merge.bam
                rm {' '.join(vals)}''' for key, vals in d1.items()]
            p.map(os.system, cmds, chunksize=1)
            result = {key: f'{self.sam_file}-{key}-merge.bam' for key, _ in d1.items()}
            result = [result[i] for i in self.group_order if i in result]
            if unmapped is not None:
                result = result + [unmapped]
            bam_list = f'{self.sam_file}-2nd-bam-list.txt'
            with open(bam_list, 'wt') as h1:
                for i in result:
                    print(i, file=h1)
            os.system(rf'''samtools cat -@ {self.nt} -b {bam_list} -o {self.sam_file}-sort-by-pos.bam
                rm {bam_list} {' '.join(result)} {gz_header}''')
        return f'{self.sam_file}-sort-by-pos.bam'
