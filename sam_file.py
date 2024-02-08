import itertools
import gzip
import os
import numpy as np
from multiprocessing import Pool
from share_func import pairwise, get_keep_bin_unmap


# from chunk_file import ChunkIndex

class SamFile:
    def __init__(self, *, sam_file, nt=16, chunk_size=None):
        self.sam_file = sam_file
        self.nt = nt
        self._sam_size = None
        self._header_size = None
        self._chunk_size = chunk_size
    
    @property
    def chunk_size(self):
        if self._chunk_size is None:
            self._chunk_size = int(self.sam_size / self.nt) + 1
        return self._chunk_size
    
    @property
    def sam_size(self):
        if self._sam_size is None:
            total = np.fromfile(f'{self.sam_file}.gzi', dtype=np.uint64)[-1]
            with os.popen(f'bgzip -b {total} {self.sam_file} | wc -c') as h1:
                self._sam_size = int(int([i.strip() for i in h1][0]) + total)
        return self._sam_size
    
    @property
    def header_size(self):
        if self._header_size is None:
            total = 0
            with gzip.open(self.sam_file, 'rt') as h1:
                for line in h1:
                    if not line.startswith('@'):
                        break
                    total += len(line)
            self._header_size = total
        return self._header_size
    
    def get_break_pos(self, *, chunk_size=None, current_pos=None):
        if current_pos is None:
            current_pos = self.header_size
        if chunk_size is None:
            chunk_size = self.chunk_size
        end = self.sam_size
        yield current_pos
        while current_pos < end:
            current_pos += chunk_size  # current to end less than chunk-size
            if current_pos >= end:
                yield end
                return
            
            with os.popen(f'bgzip -b {current_pos} {self.sam_file} | head -1 | wc -c') as h1:
                left_size = int([i.strip() for i in h1][0])
                current_pos += left_size
                if current_pos >= end:
                    yield end
                else:
                    yield current_pos
    
    @staticmethod
    def first_sort(sam_file, out_name, start_pos, query_len):
        os.system(rf'''samtools view -H {sam_file} > {out_name}.header
            bgzip -b {start_pos} -s {query_len} {sam_file} | cat {out_name}.header - | \
            samtools sort -O SAM -T {out_name}- - | \
            samtools view | bgzip -i -I {out_name}.gzi > {out_name}
            rm {out_name}.header''')
        bin_list, unmapped = get_keep_bin_unmap(out_name)
        return bin_list, unmapped
    
    def split_sam(self):
        nums = [(num1, num2 - num1) for num1, num2 in pairwise(self.get_break_pos())]
        names = [(self.sam_file, f'{self.sam_file}-chunk{num}.gz', pos1, len1) for num, (pos1, len1) in enumerate(nums)]
        with Pool(self.nt) as p:
            final = p.starmap(SamFile.first_sort, names, chunksize=1)
            keep_bin = [i[0] for i in final]
            keep_bin = list(itertools.chain.from_iterable(keep_bin))
            unmapped = [i[1] for i in final if i[1]]
            return keep_bin, unmapped
