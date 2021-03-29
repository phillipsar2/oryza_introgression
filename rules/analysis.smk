rule run_pca:
    input:
        vcf = "data/processed/filtered_snps_bpres/{REF}/oryza.{REF}.vcf.gz"
    output:
    params:
        ped_name = "data/plink/{REF}/oryza.{REF}",
        map = "data/plink/{REF}/oryza.{REF}.map",
        ped = "data/plink/{REF}/oryza.{REF}.ped",
        tmp_map = "data/pca/{REF}/{REF}.tmp.map",
        tmp_pedind = "data/pca/{REF}/{REF}.pedind",
        pedind = "data/pca/{REF}/{REF}.pedind"
    run:
        # generate plink file
        shell("plink --vcf {input.vcf} --aec --recode --out {params.ped_name}")
        # remove CM name
        shell("cat {params.map} | sed s'/CM//'g > {params.tmp_map}")
        # generate .pedind file
        shell("cat {params.ped} | cut -d ' ' -f1-5 > {params.tmp_pedind}")
        # add a population column to the file
        shell("cat {params.ped} | cut -c 1-2 | paste {params.tmp_pedind} - > {params.pedind}")
        shell("rm -i {params.tmp_pedind}")

