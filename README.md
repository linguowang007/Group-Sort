# Dependencies for the Group-Sort Package
To use the `Group-Sort` package, please ensure the following dependencies are met and that the path of these executable tools are added to the `$PATH` environment variable:
1. python (version 3.6 or higher)
2. samtools (must support the "**samtools cat -@**" parameters)
3. bgzip

# How to Use barcode_split

1. Download the source files using "**```git clone https://github.com/linguowang007/Group-Sort.git```**".

2. excute the command in terminal: "**```python Group-Sort/main.py /where/your/input.SAM.gz number-of-threads```**"

After Group-Sort excution, the sorted BAM file will generated in the same path as the input file.
