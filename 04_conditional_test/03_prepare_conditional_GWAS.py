ancestries = ["AFR", "AMR", "EAS", "EUR"]

for anc in ancestries:
    # Input/output files
    gwas_file = f"../conditional_test/01_format_gwas_summary/T2D_2023_GWAS_summary_{anc}.txt"
    cojo_file = f"02_condition_using_cojo/{anc}_T2D_COJO_cond.cma.cojo"
    output_file = f"03_prepare_conditional_GWAS/{anc}_TWAS_input.txt"

    print(f"Processing {anc}...")

    # Load GWAS allele info
    gwas_dict = {}

    with open(gwas_file) as f:
        next(f)  # skip header
        for line in f:
            items = line.strip().split()
            if len(items) < 3:
                continue

            SNP = items[0]
            A1 = items[1]
            A2 = items[2]

            gwas_dict[SNP] = (A1, A2)

    # Process COJO file
    with open(cojo_file) as f, open(output_file, "w") as out:
        header = f.readline().strip().split()
        col = {k: i for i, k in enumerate(header)}

        out.write("SNP\tA1\tA2\tZ\n")

        for line in f:
            items = line.strip().split()

            try:
                SNP = items[col["SNP"]]
                refA = items[col["refA"]]

                bC = float(items[col["bC"]])
                seC = float(items[col["bC_se"]])

                if seC == 0:
                    continue

                Z = bC / seC

                # Get original GWAS alleles
                if SNP not in gwas_dict:
                    continue

                gwas_A1, gwas_A2 = gwas_dict[SNP]

                # Align allele based on COJO refA
                if refA == gwas_A1:
                    A1 = refA
                    A2 = gwas_A2
                elif refA == gwas_A2:
                    A1 = refA
                    A2 = gwas_A1
                else:
                    continue  # allele mismatch

                out.write(f"{SNP}\t{A1}\t{A2}\t{Z}\n")

            except Exception:
                continue

    print(f"Finished {anc}: {output_file}")