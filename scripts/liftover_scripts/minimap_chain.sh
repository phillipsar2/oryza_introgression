target=
query=
out=

minimap2 -ax asm5 --eqx -t 16 $target $query -o $out
samtools view --threads 16 -b -o ${out%.*}.bam ${out}
bamToPsl -dots=100 ${out%.*}.bam ${out%.*}.psl
python -m jcvi.formats.chain frompsl ${out%.*}.psl $target $query
