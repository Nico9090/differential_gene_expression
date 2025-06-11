#vg construct \
#  -r m_chr1.fa \
#  -v m_chr1.vcf \
#  > chr1.vg

#vg index \
#        -x chr1.xg chr1.vg

#/home/ndadzie/programs/vg find \
#       -x chr1.xg \
#       -p 1:1000000-1000100 \
#       -c 100 > subgraph.vg

vg view \
        -d subgraph.vg | dot -Tsvg \
        -o subgraph.svg

#vg view \
#       -d subgraph.vg | sfdp -Tpng \
#       -o subgraph.png
