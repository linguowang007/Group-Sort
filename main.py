import sys
import os
import time
import logging
from operator import itemgetter
from sam_file import SamFile
from share_func import de_duplicates, get_chr_order
from region_maker import RegionMaker
from chunk_regions import ChunkRegions
from second_sort import SecondSort

if __name__ == '__main__':
    t1 = time.monotonic()
    
    result1 = sys.argv[1]
    nt = sys.argv[2]
    nt = int(nt)
    
    logging.basicConfig(filename=f'{result1}-sort-bgzip-slice.log', filemode='w', level=logging.INFO,
                        format='%(message)s')
    m1 = result1.endswith('.gz') and os.path.exists(f'{result1}.gzi')
    if not m1:  # not indexed bgzip sam file
        gz_sam = f'{result1}-sam.gz'
        os.system(f'samtools view -h -@ {nt} {result1} | bgzip -@ {nt} -i -I {gz_sam}.gzi > {gz_sam}')
        result1 = gz_sam
        logging.info(f'bam to bgzip sam use: {round(time.monotonic() - t1, 2)} sec.')
        t1 = time.monotonic()
    
    sam1 = SamFile(sam_file=result1, nt=nt)
    bin_list, unmapped_list = sam1.split_sam()  # chunk-file, chrome-ori, start-pos, bin-file-chrome
    
    for_chr_oder = de_duplicates(bin_list, lambda x: x[1])
    for_chr_oder = [itemgetter(0, 2)(i) for i in for_chr_oder]
    
    chrome_order = get_chr_order(sam_file=result1, chunk_pos_list=for_chr_oder)
    to_make_region = [(i[-1], chrome_order[i[1]]) for i in bin_list]
    region_mk = RegionMaker(sam_file=result1, bin_list=to_make_region, nt=nt)
    
    keep_regions, group_order = region_mk.make_regions()
    with open(f'{result1}-region.txt', 'wt') as f:
        for row in keep_regions:
            print(*row, sep='\t', file=f)
    
    to_chunk_regions = [(i[0], i[-1], chrome_order[i[1]]) for i in bin_list]  # chunk-file, bin-file, chr-id
    chunk_r1 = ChunkRegions(sam_file=result1, nt=nt, bin_list=to_chunk_regions, region_list=keep_regions,
                            group_order=group_order)
    chunk_region_list = chunk_r1.group_region()
    with open(f'{result1}-region-chunks.txt', 'wt') as f:
        for row in chunk_region_list:
            print(*row, sep='\t', file=f)
    
    unmapped = None
    unmapped_list = [i for i in unmapped_list if i]
    if unmapped_list:
        os.system(rf'''samtools view -H {result1} | bgzip > {result1}-header.gz
            cat {result1}-header.gz {' '.join(unmapped_list)} | \
            samtools view -@ {nt} -h -o {result1}-unmapped.bam -
            rm {' '.join(unmapped_list)} {result1}-header.gz''')
        unmapped = f'{result1}-unmapped.bam'
    
    sec_sort = SecondSort(sam_file=result1, region_list=chunk_region_list, group_order=group_order, nt=nt)
    final_bam = sec_sort.second_sort(unmapped=unmapped)
    rm_chunk = list({i[0] for i in bin_list})
    rm_chunk = rm_chunk + [f'{i}.gzi' for i in rm_chunk]
    os.system(rf'''rm {' '.join(rm_chunk)}''')
    logging.info(f'Sort {result1} by chr-pos use {nt} threads: {round(time.monotonic() - t1, 2)} sec.')
    os.system(f'samtools index -@ {nt} -c {final_bam}')
