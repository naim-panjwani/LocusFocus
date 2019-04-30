import subprocess

plinkrun = subprocess.run(args=[
    'plink', '--bfile', 'C:\\Users\\Naim\\Desktop\\SCS_Data_Analytics\\homework\\24-Final_project\\GWAS-QTL-Explore\\data\\EUR\\chr1'
    , "--chr", str(1)
    , "--extract", 'C:\\Users\\Naim\\Desktop\\SCS_Data_Analytics\\homework\\24-Final_project\\GWAS-QTL-Explore\\static\\session_data\\ld-10bfbeb9-24b7-45e9-8766-4169954d4071_snps.txt'
    , "--from-bp", str(205500022)
    , "--to-bp", str(205922403)
    , "--ld-snp", "chr1:205906897"
    , "--r2"
    , "--ld-window-r2", "0"
    , "--ld-window", "999999"
    , "--ld-window-kb", "200000"
    , "--make-bed"
    , "--threads", "1"
    , "--out", 'C:\\Users\\Naim\\Desktop\\SCS_Data_Analytics\\homework\\24-Final_project\\GWAS-QTL-Explore\\static\\session_data\\ld-10bfbeb9-24b7-45e9-8766-4169954d4071'
    ], stdout=subprocess.PIPE)

