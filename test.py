
def plink_ld_pairwise(lead_snp_position, pop, chrom, snp_positions, outfilename):
    # positions must be in hg19 coordinates
    # returns NaN for SNPs not in 1KG LD file; preserves order of input snp_positions
    if chrom == 'X': chrom = 23
    try:
        chrom = int(chrom)
    except:
        raise InvalidUsage("Invalid chromosome", status_code=410)
    if chrom not in np.arange(1,24):
        raise InvalidUsage("Invalid chromosome", status_code=410)
    plink_filepath = ""
    if chrom == 23:
        plink_filepath = os.path.join(MYDIR, "data", pop, "chrX")
    else:
        plink_filepath = os.path.join(MYDIR, "data", pop, f"chr{chrom}")
    # make snps file to extract:
    snps = [f"chr{str(int(chrom))}:{str(int(position))}" for position in snp_positions]
    writeList(snps, outfilename + "_snps.txt")
    plink_path = subprocess.run(args=["which","plink"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')
    plinkrun = subprocess.run(args=[
        "./plink", '--bfile', plink_filepath
        , "--chr", str(chrom)
        , "--extract", outfilename + "_snps.txt"
        , "--from-bp", str(min(snp_positions))
        , "--to-bp", str(max(snp_positions))
        , "--ld-snp", f"chr{str(int(chrom))}:{str(int(lead_snp_position))}"
        , "--r2"
        , "--ld-window-r2", "0"
        , "--ld-window", "999999"
        , "--ld-window-kb", "200000"
        , "--make-bed"
        , "--threads", "1"
        , "--out", outfilename
        ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if plinkrun.returncode != 0:
        raise InvalidUsage(plinkrun.stdout.decode('utf-8'), status_code=410)
    ld_results = pd.read_csv(outfilename + ".ld", delim_whitespace=True)
    available_r2_positions = ld_results[['BP_B', 'R2']]
    pos_df = pd.DataFrame({'pos': snp_positions})
    merged_df = pd.merge(pos_df, available_r2_positions, how='left', left_on="pos", right_on="BP_B", sort=False)[['pos', 'R2']]
    merged_df.fillna(-1, inplace=True)
    return merged_df

