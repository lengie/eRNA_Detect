### ChIPAnalysis.sh
###
###

## CBP first antibody
## filtering duplicates
macs3 filterdup -i GSM530174_CBP_Millipore_KCl_B1_E120.bedgraph --keep-dup=1 -o  GSM530174_CBP_Millipore_KCl_B1_E120.filterdup.bed -g mm --verbose 2
macs3 filterdup -i GSM530173_CBP_Millipore_un_B1_E120.bedgraph  --keep-dup=1 -o GSM530173_CBP_Millipore_un_B1_E120.filterdup.bed -g mm

# reads were already extended but let's see what the prediction is
macs3 predictd -i GSM530173_CBP_Millipore_un_B1_E120.filterdup.bed -g mm -m 5 50 

macs3 pileup -i GSM530173_CBP_Millipore_un_B1_E120.filterdup.bed -o GSM530173_CBP_Millipore_un_B1_E120.filterdup.pileup.bdg --extsize 120 -f BED
macs3 pileup -i GSM530174_CBP_Millipore_KCl_B1_E120.filterdup.bed -o GSM530174_CBP_Millipore_KCl_B1_E120.filterdup.pileup.bdg --extsize 120 -f BED

## building the local bias track from control
macs3 pileup -i GSM530173_CBP_Millipore_un_B1_E120.filterdup.bed -B --extsize 60 -o d_bg.bdg -f BED

# slocal background
macs3 pileup -i GSM530173_CBP_Millipore_un_B1_E120.filterdup.bed -B --extsize 500 -o 1k_bg.bdg -f BED
macs3 bdgopt -i 1k_bg.bdg -m multiply -p 0.12 -o 1k_bg_norm.bdg

# llocal
macs3 pileup -i GSM530173_CBP_Millipore_un_B1_E120.filterdup.bed -B --extsize 5000 -o 10k_bg.bdg -f BED
macs3 bdgopt -i 10k_bg.bdg -m multiply -p 0.012 -o 10k_bg_norm.bdg

# genome background = the_number_of_control_reads*fragment_length/genome_size = 6374811*120/1.87e9 = 0.409

## combine and generate maximum background noise
macs3 bdgcmp -m max -t 1k_bg_norm.bdg -c 10k_bg_norm.bdg -o 1k_10k_bg_norm.bdg
macs3 bdgcmp -m max -t 1k_10k_bg_norm.bdg -c d_bg.bdg -o d_1k_10k_bg_norm.bdg
macs3 bdgopt -i d_1k_10k_bg_norm.bdg -m max -p .409 -o local_bias_raw.bdg

## compare ChIP and local lambda to get scores in p-val or q-val
macs3 bdgcmp -t GSM530174_CBP_Millipore_KCl_B1_E120.filterdup.pileup.bdg -c local_lambda.bdg -m qpois -o CBP_Millipore_KCl_qvalue.bdg

## did try this as well but didn't use it
macs3 bdgcmp -t GSM530174_CBP_Millipore_KCl_B1_E120.filterdup.pileup.bdg -c local_bias_raw.bdg -m ppois -o CBP_Millipore_KCl_pvalueraw.bdg 
macs3 bdgcmp -t GSM530174_CBP_Millipore_KCl_B1_E120.filterdup.pileup.bdg -c local_lambda.bdg -m ppois -o CBP_Millipore_KCl_pvaluelambda.bdg 
macs3 bdgpeakcall -i CBP_Millipore_KCl_pvalueraw.bdg -c 1.301 -l 120 -g 75 -o CBP_Millipore_KCl_peakspvalraw.bed
macs3 bdgpeakcall -i CBP_Millipore_KCl_pvaluelambda.bdg -c 1.301 -l 120 -g 75 -o CBP_Millipore_KCl_peakspvallambda.bed


## peak calling
macs3 bdgpeakcall -i CBP_Millipore_KCl_qvalue.bdg -c 1.301 -l 120 -g 75 -o CBP_Millipore_KCl_peaksqval.bed

## CBP second antibody
## filtering duplicates
macs3 filterdup -i GSM530177_CBP_ab2832_KCl_B1_E120.bedgraph GSM530176_CBP_ab2832_KCl_B2_E120.bedgraph --keep-dup=1 -o  CBP_ab2832_KCl_E120.filterdup.bed -g mm --verbose 2
macs3 filterdup -i GSM530175_CBP_ab2832_un_B1_E120.bedgraph  --keep-dup=1 -o GSM530175_CBP_ab2832_un_B1_E120.filterdup.bed -g mm

macs3 pileup -i CBP_ab2832_KCl_E120.filterdup.bed -o CBP_ab2832_KCl_E120.filterdup.pileup.bdg --extsize 120 -f BED

## building the local bias track from control
macs3 pileup -i GSM530175_CBP_ab2832_un_B1_E120.filterdup.bed -B --extsize 60 -o d_bg.bdg -f BED

# slocal background
macs3 pileup -i GSM530175_CBP_ab2832_un_B1_E120.filterdup.bed -B --extsize 500 -o 1k_bg.bdg -f BED
macs3 bdgopt -i 1k_bg.bdg -m multiply -p 0.12 -o 1k_bg_norm.bdg

# llocal
macs3 pileup -i GSM530175_CBP_ab2832_un_B1_E120.filterdup.bed -B --extsize 5000 -o 10k_bg.bdg -f BED
macs3 bdgopt -i 10k_bg.bdg -m multiply -p 0.012 -o 10k_bg_norm.bdg

# genome background = the_number_of_control_reads*fragment_length/genome_size = 6374811*120/1.87e9 = 0.409

## combine and generate maximum background noise
macs3 bdgcmp -m max -t 1k_bg_norm.bdg -c 10k_bg_norm.bdg -o 1k_10k_bg_norm.bdg
macs3 bdgcmp -m max -t 1k_10k_bg_norm.bdg -c d_bg.bdg -o d_1k_10k_bg_norm.bdg
macs3 bdgopt -i d_1k_10k_bg_norm.bdg -m max -p .409 -o local_bias_raw.bdg

## compare ChIP and local lambda to get scores in p-val or q-val
macs3 bdgcmp -t CBP_ab2832_KCl_E120.filterdup.pileup.bdg -c local_lambda.bdg -m qpois -o CBP_ab2832_KCl_qvalue.bdg

## peak calling
macs3 bdgpeakcall -i CBP_ab2832_KCl_qvalue.bdg -c 1.301 -l 120 -g 75 -o CBP_ab2832_KCl_peaksqval.bed


