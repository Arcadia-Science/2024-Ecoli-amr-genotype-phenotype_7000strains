rule all:
    input:
        "output.txt",


rule create_file:
    output:
        "output.txt",
    shell:
        "echo 'Hello, Snakemake!' > {output}"
