import gzip

# ==============================
# 1. Load coordinate → rsID mapping
# ==============================
converter_37_in_1000g = {}

print('reading converter file...')
with open('/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/database/dbSNP/GRCh37_to_rsID.txt', 'r') as infile:
    for line in infile:
        items = line.strip().split('\t')
        if len(items) < 2:
            continue
        converter_37_in_1000g[items[0]] = items[1]

print(f"Loaded {len(converter_37_in_1000g)} mappings")


# ==============================
# 2. Convert function (MAIN)
# ==============================
def convert_gwas(pop, infile_path, outfile_path, samplesize):
    print(f'\nstart converting {pop}...')

    valid_alleles = {'A', 'T', 'C', 'G'}
    ambiguous = {('A','T'), ('T','A'), ('G','C'), ('C','G')}
    seen_snps = set()

    total = 0
    kept = 0

    with open(outfile_path, 'w') as outfile:
        outfile.write('SNP\tchr\tpos\teffect_allele\tother_allele\tbeta\tse\teaf\tsamplesize\n')

        with gzip.open(infile_path, 'rt') as f:
            next(f)  # skip header

            for line in f:
                total += 1
                items = line.strip().split('\t')

                try:
                    # ✅ parse chr:pos safely
                    chrom, pos = items[0].split(':')
                    chr_num = chrom.replace('chr', '')
                    chrpos = f"{chr_num}:{pos}"

                    A1 = items[1].upper()
                    A2 = items[2].upper()

                    freq = float(items[3])   # ✅ Freq1 = EAF
                    b = float(items[7])
                    se = float(items[8])

                except (IndexError, ValueError):
                    continue

                # ========================
                # QC filters
                # ========================

                # rsID mapping
                if chrpos not in converter_37_in_1000g:
                    continue
                SNP = converter_37_in_1000g[chrpos]

                # allele validity
                if A1 not in valid_alleles or A2 not in valid_alleles:
                    continue

                # ambiguous SNPs
                if (A1, A2) in ambiguous:
                    continue

                # frequency validity
                if freq <= 0 or freq >= 1:
                    continue

                # se validity
                if se <= 0:
                    continue

                # remove duplicates
                if SNP in seen_snps:
                    continue
                seen_snps.add(SNP)

                # ========================
                # Output
                # ========================
                outfile.write(
                    f"{SNP}\t{chr_num}\t{pos}\t{A1}\t{A2}\t{b}\t{se}\t{freq}\t{samplesize}\n"
                )
                kept += 1

    print(f"{pop}: total={total}, kept={kept}, removed={total-kept}")


# ==============================
# 3. Run per population
# ==============================
convert_gwas(
    'AFR',
    '/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/AFR_MetalFixed_LDSC-CORR_Results1TBL.gz',
    '02_prepare_GWAS/T2D_2023_GWAS_summary_AFR.txt',
    154160
)

convert_gwas(
    'EAS',
    '/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/EAS_MetalFixed_LDSC-CORR_Results1TBL.gz',
    '02_prepare_GWAS/T2D_2023_GWAS_summary_EAS.txt',
    427504
)

convert_gwas(
    'EUR',
    '/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/EUR_MetalFixed_LDSC-CORR_Results1TBL.gz',
    '02_prepare_GWAS/T2D_2023_GWAS_summary_EUR.txt',
    1812017
)

convert_gwas(
    'AMR',
    '/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/AMR_MetalFixed_LDSC-CORR_Results1TBL.gz',
    '02_prepare_GWAS/T2D_2023_GWAS_summary_AMR.txt',
    88743
)

convert_gwas(
    'SAS',
    '/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/SAS_MetalFixed_LDSC-CORR_Results1TBL.gz',
    '02_prepare_GWAS/T2D_2023_GWAS_summary_SAS.txt',
    50599
)