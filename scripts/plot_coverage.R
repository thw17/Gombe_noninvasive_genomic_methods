# Set working directory
# setwd("")


library("ggplot2")
library("dplyr")
library("cowplot")

# Extracted length of bed file with command awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' <file.bed>
pantro4_cds <- 30309842

# Raw reads for urine
s7365_urine_reads = 6905396
s7150_urine_reads = 43262642
s7072_urine_reads = 50270258
s7535_urine_reads = 49711762
s7650_urine_reads = 56850914
s7507_urine_reads = 41292318
s7323_urine_reads = 45626444

# Read histogram files
# "Raw" files (i.e., not downsampled)
# Urine only

cds_s7365_urine_exome <-
  mutate(read.table("s7365_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7365", Source="Urine", Baits="Exome", Reads=s7365_urine_reads)
cds_s7150_urine_exome <-
  mutate(read.table("s7150_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7150", Source="Urine", Baits="Exome", Reads=s7150_urine_reads)
cds_s7072_urine_exome <-
  mutate(read.table("s7072_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7072", Source="Urine", Baits="Exome", Reads=s7072_urine_reads)
cds_s7535_urine_exome <-
  mutate(read.table("s7535_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7535", Source="Urine", Baits="Exome", Reads=s7535_urine_reads)
cds_s7650_urine_exome <-
  mutate(read.table("s7650_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7650", Source="Urine", Baits="Exome", Reads=s7650_urine_reads)
cds_s7507_urine_exome <-
  mutate(read.table("s7507_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7507", Source="Urine", Baits="Exome", Reads=s7507_urine_reads)
cds_s7323_urine_exome <-
  mutate(read.table("s7323_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7323", Source="Urine", Baits="Exome", Reads=s7323_urine_reads)

# All raw
# cds_s7069_feces_exome <-
#   mutate(read.table("s7069_feces_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7069", Source="Feces", Baits="Exome")
# cds_s7069_feces_human_baits <-
#   mutate(read.table("s7069_feces_human_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7069", Source="Feces", Baits="Human_wgs")
# cds_s7069_feces_pts_baits <-
#   mutate(read.table("s7069_feces_pts_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7069", Source="Feces", Baits="Pts_wgs")
# cds_s7057_calculus_exome <-
#   mutate(read.table("s7057_calculus_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7057", Source="Calculus", Baits="Exome")
# cds_s7057_calculus_human_baits <-
#   mutate(read.table("s7057_calculus_human_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7057", Source="Calculus", Baits="Human_wgs")
# cds_s7057_calculus_shotgun <-
#   mutate(read.table("s7057_calculus_shotgun.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7057", Source="Calculus", Baits="Shotgun")
# cds_s7057_dentine_exome <-
#   mutate(read.table("s7057_dentine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7057", Source="Dentine", Baits="Exome")
# cds_s7365_feces_exome <-
#   mutate(read.table("s7365_feces_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7365", Source="Feces", Baits="Exome")
# cds_s7365_feces_human_baits <-
#   mutate(read.table("s7365_feces_human_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7365", Source="Feces", Baits="Human_wgs")
# cds_s7365_feces_pts_baits <-
#   mutate(read.table("s7365_feces_pts_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7365", Source="Feces", Baits="Pts_wgs")
# cds_s7365_urine_exome <-
#   mutate(read.table("s7365_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7365", Source="Urine", Baits="Exome")
# cds_s7150_feces_exome <-
#   mutate(read.table("s7150_feces_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7150", Source="Feces", Baits="Exome")
# cds_s7150_feces_human_baits <-
#   mutate(read.table("s7150_feces_human_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7150", Source="Feces", Baits="Human_wgs")
# cds_s7150_feces_pts_baits <-
#   mutate(read.table("s7150_feces_pts_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7150", Source="Feces", Baits="Pts_wgs")
# cds_s7150_urine_exome <-
#   mutate(read.table("s7150_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7150", Source="Urine", Baits="Exome")
# cds_s7072_urine_exome <-
#   mutate(read.table("s7072_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7072", Source="Urine", Baits="Exome")
# cds_s7535_urine_exome <-
#   mutate(read.table("s7535_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7535", Source="Urine", Baits="Exome")
# cds_s7433_calculus_exome <-
#   mutate(read.table("s7433_calculus_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7433", Source="Calculus", Baits="Exome")
# cds_s7433_calculus_human_baits <-
#   mutate(read.table("s7433_calculus_human_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7433", Source="Calculus", Baits="Human_wgs")
# cds_s7433_calculus_shotgun <-
#   mutate(read.table("s7433_calculus_shotgun.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7433", Source="Calculus", Baits="Shotgun")
# cds_s7650_urine_exome <-
#   mutate(read.table("s7650_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7650", Source="Urine", Baits="Exome")
# cds_s7507_feces_exome <-
#   mutate(read.table("s7507_feces_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7507", Source="Feces", Baits="Exome")
# cds_s7507_feces_human_baits <-
#   mutate(read.table("s7507_feces_human_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7507", Source="Feces", Baits="Human_wgs")
# cds_s7507_feces_pts_baits <-
#   mutate(read.table("s7507_feces_pts_baits.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7507", Source="Feces", Baits="Pts_wgs")
# cds_s7507_urine_exome <-
#   mutate(read.table("s7507_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7507", Source="Urine", Baits="Exome")
# cds_s7323_urine_exome <-
#   mutate(read.table("s7323_urine_exome.pantro4.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
#                     header=FALSE, sep="\t"), Sample="s7323", Source="Urine", Baits="Exome")

# Downsampled files (downsampled BAM to ~40 million reads)
cds_s7069_feces_exome_downsampled <-
  mutate(read.table("s7069_feces_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7069", Source="Feces", Baits="Exome")
cds_s7069_feces_human_baits_downsampled <-
  mutate(read.table("s7069_feces_human_baits.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7069", Source="Feces", Baits="Human_wgs")
cds_s7057_calculus_exome_downsampled <-
  mutate(read.table("s7057_calculus_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7057", Source="Calculus", Baits="Exome")
cds_s7057_calculus_human_baits_downsampled <-
  mutate(read.table("s7057_calculus_human_baits.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7057", Source="Calculus", Baits="Human_wgs")
cds_s7057_dentine_exome_downsampled <-
  mutate(read.table("s7057_dentine_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7057", Source="Dentine", Baits="Exome")
cds_s7365_feces_exome_downsampled <-
  mutate(read.table("s7365_feces_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7365", Source="Feces", Baits="Exome")
cds_s7365_feces_human_baits_downsampled <-
  mutate(read.table("s7365_feces_human_baits.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7365", Source="Feces", Baits="Human_wgs")
cds_s7365_feces_pts_baits_downsampled <-
  mutate(read.table("s7365_feces_pts_baits.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7365", Source="Feces", Baits="Pts_wgs")
cds_s7150_feces_exome_downsampled <-
  mutate(read.table("s7150_feces_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7150", Source="Feces", Baits="Exome")
cds_s7150_feces_human_baits_downsampled <-
  mutate(read.table("s7150_feces_human_baits.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7150", Source="Feces", Baits="Human_wgs")
cds_s7150_urine_exome_downsampled <-
  mutate(read.table("s7150_urine_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7150", Source="Urine", Baits="Exome")
cds_s7072_urine_exome_downsampled <-
  mutate(read.table("s7072_urine_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7072", Source="Urine", Baits="Exome")
cds_s7535_urine_exome_downsampled <-
  mutate(read.table("s7535_urine_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7535", Source="Urine", Baits="Exome")
cds_s7433_calculus_exome_downsampled <-
  mutate(read.table("s7433_calculus_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7433", Source="Calculus", Baits="Exome")
cds_s7433_calculus_human_baits_downsampled <-
  mutate(read.table("s7433_calculus_human_baits.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7433", Source="Calculus", Baits="Human_wgs")
cds_s7433_calculus_shotgun_downsampled <-
  mutate(read.table("s7433_calculus_shotgun.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7433", Source="Calculus", Baits="Shotgun")
cds_s7650_urine_exome_downsampled <-
  mutate(read.table("s7650_urine_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7650", Source="Urine", Baits="Exome")
cds_s7507_feces_exome_downsampled <-
  mutate(read.table("s7507_feces_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7507", Source="Feces", Baits="Exome")
cds_s7507_feces_human_baits_downsampled <-
  mutate(read.table("s7507_feces_human_baits.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7507", Source="Feces", Baits="Human_wgs")
cds_s7507_urine_exome_downsampled <-
  mutate(read.table("s7507_urine_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7507", Source="Urine", Baits="Exome")
cds_s7323_urine_exome_downsampled <-
  mutate(read.table("s7323_urine_exome.pantro4.downsampled.mapq20_noDup.genome_cov.bedopssorted.cds.hist",
                    header=FALSE, sep="\t"), Sample="s7323", Source="Urine", Baits="Exome")


# Add cumulative sum information
# "Raw", i.e., not downsampled
# Urine only
cds_s7365_urine_exome <- mutate(cds_s7365_urine_exome, percent= (sum(cds_s7365_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7150_urine_exome <- mutate(cds_s7150_urine_exome, percent= (sum(cds_s7150_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7072_urine_exome <- mutate(cds_s7072_urine_exome, percent= (sum(cds_s7072_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7535_urine_exome <- mutate(cds_s7535_urine_exome, percent= (sum(cds_s7535_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7650_urine_exome <- mutate(cds_s7650_urine_exome, percent= (sum(cds_s7650_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7507_urine_exome <- mutate(cds_s7507_urine_exome, percent= (sum(cds_s7507_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7323_urine_exome <- mutate(cds_s7323_urine_exome, percent= (sum(cds_s7323_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)

# All raw samples
# cds_s7069_feces_exome <- mutate(cds_s7069_feces_exome, percent= (sum(cds_s7069_feces_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7069_feces_human_baits <- mutate(cds_s7069_feces_human_baits, percent= (sum(cds_s7069_feces_human_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7069_feces_pts_baits <- mutate(cds_s7069_feces_pts_baits, percent= (sum(cds_s7069_feces_pts_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7057_calculus_exome <- mutate(cds_s7057_calculus_exome, percent= (sum(cds_s7057_calculus_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7057_calculus_human_baits <- mutate(cds_s7057_calculus_human_baits, percent= (sum(cds_s7057_calculus_human_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7057_calculus_shotgun <- mutate(cds_s7057_calculus_shotgun, percent= (sum(cds_s7057_calculus_shotgun$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7057_dentine_exome <- mutate(cds_s7057_dentine_exome, percent= (sum(cds_s7057_dentine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7365_feces_exome <- mutate(cds_s7365_feces_exome, percent= (sum(cds_s7365_feces_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7365_feces_human_baits <- mutate(cds_s7365_feces_human_baits, percent= (sum(cds_s7365_feces_human_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7365_feces_pts_baits <- mutate(cds_s7365_feces_pts_baits, percent= (sum(cds_s7365_feces_pts_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7365_urine_exome <- mutate(cds_s7365_urine_exome, percent= (sum(cds_s7365_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7150_feces_exome <- mutate(cds_s7150_feces_exome, percent= (sum(cds_s7150_feces_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7150_feces_human_baits <- mutate(cds_s7150_feces_human_baits, percent= (sum(cds_s7150_feces_human_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7150_feces_pts_baits <- mutate(cds_s7150_feces_pts_baits, percent= (sum(cds_s7150_feces_pts_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7150_urine_exome <- mutate(cds_s7150_urine_exome, percent= (sum(cds_s7150_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7072_urine_exome <- mutate(cds_s7072_urine_exome, percent= (sum(cds_s7072_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7535_urine_exome <- mutate(cds_s7535_urine_exome, percent= (sum(cds_s7535_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7433_calculus_exome <- mutate(cds_s7433_calculus_exome, percent= (sum(cds_s7433_calculus_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7433_calculus_human_baits <- mutate(cds_s7433_calculus_human_baits, percent= (sum(cds_s7433_calculus_human_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7433_calculus_shotgun <- mutate(cds_s7433_calculus_shotgun, percent= (sum(cds_s7433_calculus_shotgun$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7650_urine_exome <- mutate(cds_s7650_urine_exome, percent= (sum(cds_s7650_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7507_feces_exome <- mutate(cds_s7507_feces_exome, percent= (sum(cds_s7507_feces_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7507_feces_human_baits <- mutate(cds_s7507_feces_human_baits, percent= (sum(cds_s7507_feces_human_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7507_feces_pts_baits <- mutate(cds_s7507_feces_pts_baits, percent= (sum(cds_s7507_feces_pts_baits$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7507_urine_exome <- mutate(cds_s7507_urine_exome, percent= (sum(cds_s7507_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)
# cds_s7323_urine_exome <- mutate(cds_s7323_urine_exome, percent= (sum(cds_s7323_urine_exome$V2) - cumsum(V2) + V2) / pantro4_cds)

# Downsampled
cds_s7069_feces_exome_downsampled <- mutate(cds_s7069_feces_exome_downsampled, percent= (sum(cds_s7069_feces_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7069_feces_human_baits_downsampled <- mutate(cds_s7069_feces_human_baits_downsampled, percent= (sum(cds_s7069_feces_human_baits_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7057_calculus_exome_downsampled <- mutate(cds_s7057_calculus_exome_downsampled, percent= (sum(cds_s7057_calculus_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7057_calculus_human_baits_downsampled <- mutate(cds_s7057_calculus_human_baits_downsampled, percent= (sum(cds_s7057_calculus_human_baits_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7057_dentine_exome_downsampled <- mutate(cds_s7057_dentine_exome_downsampled, percent= (sum(cds_s7057_dentine_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7365_feces_exome_downsampled <- mutate(cds_s7365_feces_exome_downsampled, percent= (sum(cds_s7365_feces_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7365_feces_human_baits_downsampled <- mutate(cds_s7365_feces_human_baits_downsampled, percent= (sum(cds_s7365_feces_human_baits_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7365_feces_pts_baits_downsampled <- mutate(cds_s7365_feces_pts_baits_downsampled, percent= (sum(cds_s7365_feces_pts_baits_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7150_feces_exome_downsampled <- mutate(cds_s7150_feces_exome_downsampled, percent= (sum(cds_s7150_feces_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7150_feces_human_baits_downsampled <- mutate(cds_s7150_feces_human_baits_downsampled, percent= (sum(cds_s7150_feces_human_baits_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7150_urine_exome_downsampled <- mutate(cds_s7150_urine_exome_downsampled, percent= (sum(cds_s7150_urine_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7072_urine_exome_downsampled <- mutate(cds_s7072_urine_exome_downsampled, percent= (sum(cds_s7072_urine_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7535_urine_exome_downsampled <- mutate(cds_s7535_urine_exome_downsampled, percent= (sum(cds_s7535_urine_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7433_calculus_exome_downsampled <- mutate(cds_s7433_calculus_exome_downsampled, percent= (sum(cds_s7433_calculus_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7433_calculus_human_baits_downsampled <- mutate(cds_s7433_calculus_human_baits_downsampled, percent= (sum(cds_s7433_calculus_human_baits_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7433_calculus_shotgun_downsampled <- mutate(cds_s7433_calculus_shotgun_downsampled, percent= (sum(cds_s7433_calculus_shotgun_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7650_urine_exome_downsampled <- mutate(cds_s7650_urine_exome_downsampled, percent= (sum(cds_s7650_urine_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7507_feces_exome_downsampled <- mutate(cds_s7507_feces_exome_downsampled, percent= (sum(cds_s7507_feces_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7507_feces_human_baits_downsampled <- mutate(cds_s7507_feces_human_baits_downsampled, percent= (sum(cds_s7507_feces_human_baits_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7507_urine_exome_downsampled <- mutate(cds_s7507_urine_exome_downsampled, percent= (sum(cds_s7507_urine_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)
cds_s7323_urine_exome_downsampled <- mutate(cds_s7323_urine_exome_downsampled, percent= (sum(cds_s7323_urine_exome_downsampled$V2) - cumsum(V2) + V2) / pantro4_cds)

# Combine dataframes

cds_pantro4_calculus_dentin_downsampled <- rbind(
  cds_s7057_calculus_exome_downsampled, cds_s7057_dentine_exome_downsampled,
  cds_s7433_calculus_exome_downsampled)

cds_pantro4_feces_urine_exome_downsampled <- rbind(
  cds_s7069_feces_exome_downsampled, cds_s7365_feces_exome_downsampled,
  cds_s7150_feces_exome_downsampled, cds_s7150_urine_exome_downsampled,
  cds_s7072_urine_exome_downsampled, cds_s7535_urine_exome_downsampled,
  cds_s7650_urine_exome_downsampled, cds_s7507_feces_exome_downsampled,
  cds_s7507_urine_exome_downsampled, cds_s7323_urine_exome_downsampled)

cds_pantro4_combined_feces_downsampled <- rbind(
  cds_s7069_feces_exome_downsampled, cds_s7069_feces_human_baits_downsampled,
  cds_s7365_feces_exome_downsampled, cds_s7365_feces_human_baits_downsampled,
  cds_s7365_feces_pts_baits_downsampled, cds_s7150_feces_exome_downsampled,
  cds_s7150_feces_human_baits_downsampled, cds_s7507_feces_exome_downsampled,
  cds_s7507_feces_human_baits_downsampled)

cds_pantro4_combined_urine <- rbind(
  cds_s7365_urine_exome, cds_s7150_urine_exome,
  cds_s7072_urine_exome, cds_s7535_urine_exome,
  cds_s7650_urine_exome, cds_s7507_urine_exome, cds_s7323_urine_exome)

# Reads vs. coverage (urine--not downsampled) association
dp1_urine <- filter(cds_pantro4_combined_urine, V1 == 1)
ggplot(dp1_urine, aes(x=Reads, y=percent)) + geom_point() + geom_smooth(method='lm')

dp8_urine <- filter(cds_pantro4_combined_urine, V1 == 8)
ggplot(dp8_urine, aes(x=Reads, y=percent)) + geom_point() + geom_smooth(method='lm')

# Make plots
cds_pantro4_s7323_urine_plot <-
  ggplot(cds_s7323_urine_exome_downsampled, aes(x = V1, y= percent)) +
  geom_line(size = .5, alpha = 0.75) + scale_color_brewer(palette="Set1") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) +
  geom_vline(xintercept=4, linetype = "solid") + geom_vline(xintercept=8, linetype = "solid") +
  theme_bw() +
  xlab(label = "Depth of coverage") +
  theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) +
  ylab(label = "Proportion of cds in genome\nat X coverage or greater") +
  theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle("Callable sites in PanTro4 Ensembl CDS") +
  theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) +
  theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) +
  labs(color = "Sample (Urine)") + theme(legend.title.align=0.5)
ggsave(filename="Urine_exome_gombe_coverage_s7323only.pdf", plot=cds_pantro4_s7323_urine_plot, width=10, height=10, dpi = 800)

cds_pantro4_combined_urine_plot <-
  ggplot(cds_pantro4_combined_urine, aes(x = V1, y= percent, color= Sample)) +
  geom_line(size = .5, alpha = 0.75) + scale_color_brewer(palette="Set1") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) +
  geom_vline(xintercept=4, linetype = "solid") + geom_vline(xintercept=8, linetype = "solid") +
  theme_bw() +
  xlab(label = "Depth of coverage") +
  theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) +
  ylab(label = "Proportion of cds in genome\nat X coverage or greater") +
  theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle("Callable sites in PanTro4 Ensembl CDS") +
  theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) +
  theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) +
  labs(color = "Sample (Urine)") + theme(legend.title.align=0.5)

cds_pantro4_feces_vs_urine_plot <-
  ggplot(cds_pantro4_feces_urine_exome_downsampled, aes(x = V1, y= percent, color= Sample)) +
  geom_line(aes(linetype=Source), size = 1, alpha = 0.75) + scale_color_brewer(palette="Set1") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) +
  geom_vline(xintercept=4, linetype = "solid") + geom_vline(xintercept=8, linetype = "solid") +
  theme_bw() +
  xlab(label = "\nDepth of coverage") +
  theme(axis.title.x = element_text(margin = margin(t = 1), size = 16)) +
  ylab(label = "Proportion of cds in genome at X coverage or greater\n") +
  theme(axis.title.y = element_text(size = 16)) + theme(axis.text.x=element_text(size = 14, angle=45, hjust=1)) +
  theme(axis.text.y=element_text(size = 14)) +
  ggtitle("Callable sites in PanTro4 Ensembl CDS") +
  theme(plot.title = element_text(hjust = .5, size = 18, face="bold")) +
  theme(legend.text=element_text(size=14), legend.title=element_text(size=16, face="bold")) +
  labs(color = "Sample") + theme(legend.title.align=0.5)

cds_pantro4_calculs_vs_dentin_plot <-
  ggplot(cds_pantro4_calculus_dentin_downsampled, aes(x = V1, y= percent, color= Sample)) +
  geom_line(aes(linetype=Source), size = 1, alpha = 0.75) + scale_color_brewer(palette="Set1") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) +
  geom_vline(xintercept=4, linetype = "solid") + geom_vline(xintercept=8, linetype = "solid") +
  theme_bw() +
  xlab(label = "\nDepth of coverage") +
  theme(axis.title.x = element_text(margin = margin(t = 1), size = 16)) +
  ylab(label = "Proportion of cds in genome at X coverage or greater\n") +
  theme(axis.title.y = element_text(size = 16)) + theme(axis.text.x=element_text(size = 14, angle=45, hjust=1)) +
  theme(axis.text.y=element_text(size = 14)) +
  ggtitle("Callable sites in PanTro4 Ensembl CDS") +
  theme(plot.title = element_text(hjust = .5, size = 18, face="bold")) +
  theme(legend.text=element_text(size=14), legend.title=element_text(size=16, face="bold")) +
  labs(color = "Sample") + theme(legend.title.align=0.5)

cds_pantro4_bait_comparison_plot <-
  ggplot(cds_pantro4_combined_feces_downsampled, aes(x = V1, y= percent, color= Sample)) +
  geom_line(aes(linetype=Baits), size = 1, alpha = 0.75) + scale_color_brewer(palette="Set1") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 12), breaks = seq(0,12,2), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) +
  theme_bw() +
  xlab(label = "\nDepth of coverage") +
  theme(axis.title.x = element_text(margin = margin(t = 1), size = 16)) +
  ylab(label = "Proportion of cds in genome at X coverage or greater\n") +
  theme(axis.title.y = element_text(size = 16)) + theme(axis.text.x=element_text(size = 14, angle=45, hjust=1)) +
  theme(axis.text.y=element_text(size = 14)) +
  ggtitle("Callable sites in PanTro4 Ensembl CDS") +
  theme(plot.title = element_text(hjust = .5, size = 18, face="bold")) +
  theme(legend.text=element_text(size=16),legend.title=element_text(size=14, face="bold")) +
  labs(color = "Sample") + theme(legend.title.align=0.5)

# Save plots
ggsave(filename="Calculus_vs_dentin_exome_gombe_coverage_comparison.pdf", plot=cds_pantro4_calculs_vs_dentin_plot, width=10, height=10, dpi = 800)
ggsave(filename="Urine_exome_gombe_coverage_comparison.pdf", plot=cds_pantro4_combined_urine_plot, width=10, height=10, dpi = 800)
ggsave(filename="Feces_vs_urine_exome_gombe_coverage_comparison.pdf", plot=cds_pantro4_feces_vs_urine_plot, width=10, height=10, dpi = 800)
ggsave(filename="Bait_comparison_gombe_coverage_comparison.pdf", plot=cds_pantro4_bait_comparison_plot, width=10, height=10, dpi = 800)


#########################################
# Currently not using these
# cds_pantro4_combined <- rbind(
#   cds_s7069_feces_exome, cds_s7069_feces_human_baits, cds_s7069_feces_pts_baits,
#   cds_s7057_calculus_exome, cds_s7057_calculus_human_baits, cds_s7057_calculus_shotgun,
#   cds_s7057_dentine_exome, cds_s7365_feces_exome, cds_s7365_feces_human_baits,
#   cds_s7365_feces_pts_baits, cds_s7365_urine_exome, cds_s7150_feces_exome,
#   cds_s7150_feces_human_baits, cds_s7150_feces_pts_baits, cds_s7150_urine_exome,
#   cds_s7072_urine_exome, cds_s7535_urine_exome, cds_s7433_calculus_exome,
#   cds_s7433_calculus_human_baits, cds_s7433_calculus_shotgun, cds_s7650_urine_exome,
#   cds_s7507_feces_exome, cds_s7507_feces_human_baits, cds_s7507_feces_pts_baits,
#   cds_s7507_urine_exome, cds_s7323_urine_exome)
#
# cds_pantro4_combined_urine <- rbind(
#   cds_s7365_urine_exome, cds_s7150_urine_exome,
#   cds_s7072_urine_exome, cds_s7535_urine_exome,
#   cds_s7650_urine_exome, cds_s7507_urine_exome, cds_s7323_urine_exome)
#
# cds_pantro4_combined_feces <- rbind(
#   cds_s7069_feces_exome, cds_s7069_feces_human_baits, cds_s7069_feces_pts_baits,
#   cds_s7365_feces_exome, cds_s7365_feces_human_baits, cds_s7365_feces_pts_baits,
#   cds_s7150_feces_exome, cds_s7150_feces_human_baits, cds_s7150_feces_pts_baits,
#   cds_s7507_feces_exome, cds_s7507_feces_human_baits, cds_s7507_feces_pts_baits)
################################################

# # Plot
# cds_pantro4_combined_plot <-
#   ggplot(cds_pantro4_combined, aes(x = V1, y= percent, color= Sample)) +
#   geom_line(size = .5, alpha = 0.75) +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +
#   scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) +
#   geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") +
#   geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.958) + theme_bw() +
#   xlab(label = "Depth of coverage") +
#   theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) +
#   ylab(label = "Proportion of cds in genome\nat X coverage or greater") +
#   theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) +
#   ggtitle("Callable sites in PanTro4 CDS") +
#   theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) +
#   theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) +
#   labs(color = "Sample") + theme(legend.title.align=0.5)
#
# cds_pantro4_combined_urine_plot <-
#   ggplot(cds_pantro4_combined_urine, aes(x = V1, y= percent, color= Sample)) +
#   geom_line(size = .5, alpha = 0.75) +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +
#   scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) +
#   geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") +
#   geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.958) + theme_bw() +
#   xlab(label = "Depth of coverage") +
#   theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) +
#   ylab(label = "Proportion of cds in genome\nat X coverage or greater") +
#   theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) +
#   ggtitle("Callable sites in PanTro4 Ensembl CDS") +
#   theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) +
#   theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) +
#   labs(color = "Sample (Urine)") + theme(legend.title.align=0.5)
#
# cds_pantro4_combined_feces_plot <-
#   ggplot(cds_pantro4_combined_feces, aes(x = V1, y= percent, color= Sample)) +
#   geom_line(size = .5, alpha = 0.75) +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +
#   scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) +
#   geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") +
#   geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.958) + theme_bw() +
#   xlab(label = "Depth of coverage") +
#   theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) +
#   ylab(label = "Proportion of cds in genome\nat X coverage or greater") +
#   theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) +
#   ggtitle("Callable sites in PanTro4 Ensembl CDS") +
#   theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) +
#   theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) +
#   labs(color = "Sample (Feces)") + theme(legend.title.align=0.5)
#
# ggsave(filename="Urine_exome_gombe_coverage_comparison.png", plot=cds_pantro4_combined_urine_plot, width=10, height=10, dpi = 800)
# ggsave(filename="Feces_exome_gombe_coverage_comparison.png", plot=cds_pantro4_combined_feces_plot, width=10, height=10, dpi = 800)
