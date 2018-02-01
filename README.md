# IsolateNS
This is a small program for dN/dS ratio calculation based on VCF file and reference gbk/gbf file.

Most of the functions in this program are mainly based on codes from other people. Thanks to Mr. Adel Qalieh who wrote 'dN/dS Calculator' which is available from https://github.com/adelq/dnds.

To use this program you need to input a VCF file to specify your SNP (only SNP is accept! no indel allow!) and a genbank file in gbk/gbf format to specify positions of genes.

The mainly usage is as following:

python WooIsolateNS.py -q file.vcf -r reference.gbf

Enter -h or --help would give following information:

-q --query     Your SNP file (in VCF format)

-r --reference Your reference file (in gbk or gbf format)

The output format would be as following:

File_name chromosome pN pS pN/pS dN dS dN/dS

Notice that if both dN and dS is 0, the ratio would be '/'; if only dN is 0, the ratio would be '-'; if only dS is 0, the ratio would be '+'.

The codes would be uploaded soon after final modifications.
