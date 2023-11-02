import gtf_reader
annotation = gtf_reader.gtf_reader(f_file=open("/mnt/ddngs/zhangsw/database/GRCh37/gencode.v38lift37.annotation.gtf"), type="exon", s_id="gene_name", with_chr=False, add_all_isoform=True)

out = open("gencode.v38lift37.annotation.exon_total_len.tsv", "w")
out.write("Gene\tLength\n")
for gene in annotation.keys():
    for idx in annotation[gene]:
        info = annotation[gene][idx]
        total_len = 0
        for region in info.region:
            total_len += region[-1] - region[0] + 1
        out.write(gene + "\t" + str(total_len) + "\n")

out.close()

