{\rtf1\ansi\ansicpg1252\deff0\nouicompat\deflang1033{\fonttbl{\f0\fnil\fcharset0 Calibri;}}
{\*\generator Riched20 10.0.17134}\viewkind4\uc1 
\pard\sa200\sl276\slmult1\f0\fs22\lang9 import pysam\par
import sys\par
import copy\par
\par
tmpfilename                   = "/mnt/MiSeq/Jolly/quantgene/Synthetic_datasets/chk_all3.sam"\par
hotspot_filename              = "/mnt/MiSeq/Jolly/quantgene/Synthetic_datasets/Quantgene_Hotspots_V3_filtered.txt"\par
SNP_reference_fasta_filename  = "/mnt/MiSeq/Jolly/quantgene/Synthetic_datasets/SNPs_ref_checked.fasta"\par
SNP_alternate_fasta_filename  = "/mnt/MiSeq/Jolly/quantgene/Synthetic_datasets/SNPs_alt_checked.fasta"\par
read_name_file                = "/mnt/MiSeq/Jolly/quantgene/Synthetic_datasets/read_names"\par
\par
\par
header = \{ 'HD': \{'VN': '1.0'\},\par
            'SQ': [\{'LN': 249250621 , 'SN': 'chr1'\},\par
                   \{'LN': 135534747, 'SN': 'chr10'\},\par
                   \{'LN': 135006516, 'SN': 'chr11'\},\par
                   \{'LN': 133851895, 'SN': 'chr12'\},\par
                   \{'LN': 107349540, 'SN': 'chr14'\},\par
                   \{'LN': 102531392, 'SN': 'chr15'\},\par
                   \{'LN': 90354753, 'SN': 'chr16'\},\par
                   \{'LN': 81195210, 'SN': 'chr17'\},\par
                   \{'LN': 78077248, 'SN': 'chr18'\},\par
                   \{'LN': 59128983, 'SN': 'chr19'\},\par
                   \{'LN': 63025520, 'SN': 'chr20'\},\par
                   \{'LN': 51304566, 'SN': 'chr22'\},\par
                   \{'LN': 198022430, 'SN': 'chr3'\},\par
                   \{'LN': 180915260, 'SN': 'chr5'\},\par
                   \{'LN': 171115067, 'SN': 'chr6'\},\par
                   \{'LN': 159138663, 'SN': 'chr7'\},\par
                   \{'LN': 146364022, 'SN': 'chr8'\},\par
                   \{'LN': 141213431, 'SN': 'chr9'\},\par
                   \{'LN': 155270560, 'SN': 'chrX'\}] \}\par
\par
SNP_ref_mutation_name = []\par
SNP_ref_fasta         = []\par
SNP_alt_mutation_name = []\par
SNP_alt_fasta         = []\par
read_names            = []\par
\par
#read the fasta file with reference SNP sequence\par
with open(SNP_reference_fasta_filename) as ref_fasta:\par
 for f in ref_fasta:\par
   f       = f.strip('\\n')\par
   content = f.split("\\t")\par
   SNP_ref_mutation_name.append(content[0])\par
   SNP_ref_fasta.append(content[1])\par
\par
#read the fasta file with alternate SNP sequence\par
with open(SNP_alternate_fasta_filename) as alt_fasta:\par
 for f in alt_fasta:\par
   f       = f.strip('\\n')\par
   content = f.split("\\t")\par
   SNP_alt_mutation_name.append(content[0])\par
   SNP_alt_fasta.append(content[1])\par
\par
\par
#read the read name file to insert as read name in bam file so that we can do umi analysis\par
with open(read_name_file) as reads:\par
 for f in reads:\par
   f       = f.strip('\\n')\par
   read_names.append(f)\par
\par
with pysam.AlignmentFile(tmpfilename, "wb", header = header) as outf:\par
# open the reference file for all SNP hotspots\par
 with open(hotspot_filename) as hotspots:\par
    for line in hotspots:\par
     line = line.strip('\\n')\par
     location = line.split("\\t")\par
     reference_id = 0\par
     print(location[0], ":", location[1], "-", location[2])\par
     # get the chromosome location from hotspot file and extract the header dictionary for the index\par
     for kk in header:\par
         if kk.title() == "Sq":\par
             list_of_chromosomes = header[kk]\par
             for chr in range(len(list_of_chromosomes)):\par
                 for key in list_of_chromosomes[chr]:\par
                     if list_of_chromosomes[chr][key] == location[0]:\par
                         reference_id = chr\par
                         break\par
\par
     # get the reference sequence based on mutation we are interested\par
     for name in SNP_ref_mutation_name:\par
         if name in location[4]:\par
             j = SNP_ref_mutation_name.index(name)\par
             query_sequence_ref = SNP_ref_fasta[j]\par
             length1 = len(query_sequence_ref)\par
             break\par
\par
\par
\par
      # get the alternate sequence based on mutation we are interested\par
     for name in SNP_alt_mutation_name:\par
         if name in location[4]:\par
            j = SNP_alt_mutation_name.index(name)\par
            query_sequence_alt = SNP_alt_fasta[j]\par
            length2 = len(query_sequence_alt)\par
\par
     if location[4].find('del') != -1:\par
             print("Deletion found")\par
             deletion_length = int(location[2]) - int(location[1])\par
             if location[1] > location[3]:\par
                 match_length1 = int(location[1]) - int(location[3]) + 1\par
                 match_length2 = length2 - match_length1\par
                 query_start = int(location[3])\par
             elif location[1] < location[3]:\par
                 match_length2 = int(location[3]) - (int(location[1]) + 1)\par
                 match_length1 = length2 - match_length2\par
                 query_start = int(location[3]) - length2\par
             cig_temp1 = str(match_length2) + "M" + str(deletion_length) + "D" + str(match_length1) + "M"\par
\par
             # for reference Deletion processing with 15 reads matching all criteria\par
             for x in range(5):\par
                 query_name1 = read_names.pop()\par
                 for y in range(3):\par
                     a = pysam.AlignedSegment()\par
                     a.query_name = query_name1\par
                     a.query_sequence = query_sequence_ref\par
                     a.flag = 99\par
                     a.reference_id = reference_id\par
                     a.reference_start = query_start\par
                     a.mapping_quality = 20\par
                     cig_temp = str(length1) + "M"\par
                     a.cigarstring = (cig_temp)\par
                     a.next_reference_id = reference_id\par
                     a.next_reference_start = int(location[3]) + 20\par
                     a.template_length = length1\par
                     quality = "<" * length1\par
                     a.query_qualities = pysam.qualitystring_to_array(quality)\par
                     a.tags = (("NM", 0), ("RG", "L1"))\par
                     outf.write(a)\par
\par
             for x in range(5):\par
                query_name2 = read_names.pop()\par
                for y in range(3):\par
                     b = pysam.AlignedSegment()\par
                     b.query_name = query_name2\par
                     b.query_sequence = query_sequence_alt\par
                     b.flag = 99\par
                     b.reference_id = reference_id\par
                     b.reference_start = query_start\par
                     b.mapping_quality = 60\par
                     b.cigarstring = (cig_temp1)\par
                     b.next_reference_id = reference_id\par
                     b.next_reference_start = int(location[3]) + 20\par
                     b.template_length = length2\par
                     quality = "<" * length2\par
                     b.query_qualities = pysam.qualitystring_to_array(quality)\par
                     b.tags = (("NM", 1), ("RG", "L1"))\par
                   #  print(b.query_name,"\\t",b.reference_start, "\\t", b.cigarstring, "\\t", b.query_sequence)\par
                     outf.write(b)\par
\par
     elif location[4].find('ins') != -1:\par
             insertion_length = int(location[2]) - int(location[1])\par
             if location[1] > location[3]:\par
                 match_length1 = int(location[1]) - int(location[3])\par
                 match_length2 = length2 - match_length1 - insertion_length\par
                 query_start = int(location[3])\par
             elif location[1] < location[3]:\par
                 match_length2 = int(location[3]) - (int(location[1]) + 1)\par
                 match_length1 = length2 - match_length2 - insertion_length\par
                 query_start = int(location[1]) - (int(location[3]) - int(location[1]))\par
             cig_temp1 = str(match_length2) + "M" + str(insertion_length) + "I" + str(match_length1) + "M"\par
\par
             # for reference insertion processing\par
             for x in range(5):\par
                 query_name1 = read_names.pop()\par
                 for y in range(3):\par
                     a = pysam.AlignedSegment()\par
                     a.query_name = query_name1\par
                     a.query_sequence = query_sequence_ref\par
                     a.flag = 99\par
                     a.reference_id = reference_id\par
                     a.reference_start = query_start\par
                     a.mapping_quality = 20\par
                     cig_temp = str(length1) + "M"\par
                     a.cigarstring = (cig_temp)\par
                     a.next_reference_id = reference_id\par
                     a.next_reference_start = int(location[3]) + 20\par
                     a.template_length = length1\par
                     quality = "<" * length1\par
                     a.query_qualities = pysam.qualitystring_to_array(quality)\par
                     a.tags = (("NM", 0), ("RG", "L1"))\par
                     outf.write(a)\par
\par
             for x in range(5):\par
                query_name2 = read_names.pop()\par
                for y in range(3):\par
                     b = pysam.AlignedSegment()\par
                     b.query_name = query_name2\par
                     b.query_sequence = query_sequence_alt\par
                     b.flag = 99\par
                     b.reference_id = reference_id\par
                     b.reference_start = query_start\par
                     b.mapping_quality = 60\par
                     b.cigarstring = (cig_temp1)\par
                     b.next_reference_id = reference_id\par
                     b.next_reference_start = int(location[3]) + 20\par
                     b.template_length = length2\par
                     quality = "<" * length2\par
                     b.query_qualities = pysam.qualitystring_to_array(quality)\par
                     b.tags = (("NM", 1), ("RG", "L1"))\par
                 #    print(b.query_name, "\\t", b.reference_start, "\\t", b.cigarstring, "\\t", b.query_sequence)\par
                     outf.write(b)\par
\par
     elif location[4] not in('del','ins'):\par
          print("SNP processing \\t")\par
          if location[1] > location[3]:\par
              query_start  =  int(location[3])\par
              match_length1 = int(location[1]) - int(location[3]) -1\par
              match_length2 = length2 - match_length1 -1\par
              query_start1 =  int(location[3]) - 20\par
              #quality calculation\par
              quality1 = "<"*match_length1+"!"+">"*match_length2\par
              print(quality1)\par
          elif location[1] < location[3]:\par
              query_start = int(location[3]) - length1\par
              query_start1 = int(location[3]) + 20\par
              match_length2 = int(location[3]) - int(location[2])\par
              match_length1 = length2 - match_length1 - 1\par
              # quality calculation\par
              quality1 = "<" * match_length1 + "!" + ">" * match_length2\par
              #bad primer location\par
              query_start1 = int(location[3]) - length1 - 20\par
\par
          #for normal SNP processing\par
          for x in range(5):  # for reference reads with good primer start, mapping quality, base quality and no sec/suppl alignments\par
           query_name1 = read_names.pop()\par
           for y in range(3):\par
             a = pysam.AlignedSegment()\par
             a.query_name = query_name1\par
             a.query_sequence=query_sequence_ref\par
             a.flag = 99\par
             a.reference_id = reference_id\par
             a.reference_start = query_start\par
             a.mapping_quality = 40\par
             cig_temp=str(length1)+"M"\par
             a.cigarstring = (cig_temp)\par
             a.next_reference_id = reference_id\par
             a.next_reference_start=int(location[3]) + 20\par
             a.template_length= length1\par
             quality = "<"*length1\par
             a.query_qualities = pysam.qualitystring_to_array(quality)\par
             a.tags = (("NM", 0),("RG", "L1"))\par
             outf.write(a)\par
\par
          # for alternate sequence with right primer, base quality, mapping quality, no sec or suppl alignments\par
          for x in range(5):\par
            query_name1 = read_names.pop()\par
            for y in range(3):\par
             b = pysam.AlignedSegment()\par
             b.query_name = query_name1\par
             b.query_sequence = query_sequence_alt\par
             b.flag = 99\par
             b.reference_id = reference_id\par
             b.reference_start = query_start\par
             b.mapping_quality = 60\par
             cig_temp1 = str(length2) + "M"\par
             b.cigarstring = (cig_temp1)\par
             b.next_reference_id = reference_id\par
             b.next_reference_start = int(location[3]) + 20\par
             b.template_length = length2\par
             quality = "<" * length2\par
             b.query_qualities = pysam.qualitystring_to_array(quality)\par
             b.tags = (("NM", 1), ("RG", "L1"))\par
             outf.write(b)\par
\par
          # for alternate sequence with incorrect primer start but good base quality, mapping quality, no sec or suppl alignments\par
          for x in range(5):\par
             query_name1 = read_names.pop()\par
             b = pysam.AlignedSegment()\par
             b.query_name = query_name1\par
             b.query_sequence = query_sequence_ref\par
             b.flag = 99\par
             b.reference_id = reference_id\par
             b.reference_start = query_start1\par
             b.mapping_quality = 60\par
             cig_temp1 = str(length1) + "M"\par
             b.cigarstring = (cig_temp1)\par
             b.next_reference_id = reference_id\par
             b.next_reference_start = int(location[3]) + 20\par
             b.template_length = length1\par
             quality = "<" * length1\par
             b.query_qualities = pysam.qualitystring_to_array(quality)\par
             b.tags = (("NM", 1), ("RG", "L1"))\par
             outf.write(b)\par
\par
         # for alternate sequence with correct primer start, good base quality, bad mapping quality, no sec or suppl alignments\par
          for x in range(5):\par
             query_name1 = read_names.pop()\par
             b = pysam.AlignedSegment()\par
             b.query_name = query_name1\par
             b.query_sequence = query_sequence_ref\par
             b.flag = 99\par
             b.reference_id = reference_id\par
             b.reference_start = query_start\par
             b.mapping_quality = 12\par
             cig_temp1 = str(length2) + "M"\par
             b.cigarstring = (cig_temp1)\par
             b.next_reference_id = reference_id\par
             b.next_reference_start = int(location[3]) + 20\par
             b.template_length = length1\par
             quality = "<" * length1\par
             b.query_qualities = pysam.qualitystring_to_array(quality)\par
             b.tags = (("NM", 1), ("RG", "L1"))\par
             outf.write(b)\par
\par
        # for alternate sequence with correct primer start, good base quality, good mapping quality, but sec alignments\par
          for x in range(5):\par
             query_name1 = read_names.pop()\par
             b = pysam.AlignedSegment()\par
             b.query_name = query_name1\par
             b.query_sequence = query_sequence_ref\par
             b.flag = 323\par
             b.reference_id = reference_id\par
             b.reference_start = query_start\par
             b.mapping_quality = 60\par
             cig_temp1 = str(length2) + "M"\par
             b.cigarstring = (cig_temp1)\par
             b.next_reference_id = reference_id\par
             b.next_reference_start = int(location[3]) + 20\par
             b.template_length = length2\par
             quality = "<" * length2\par
             b.query_qualities = pysam.qualitystring_to_array(quality)\par
             b.tags = (("NM", 1), ("RG", "L1"))\par
             outf.write(b)\par
\par
        # for alternate sequence with correct primer start, good base quality, good mapping quality, but suppl alignments\par
          for x in range(5):\par
             query_name1 = read_names.pop()\par
             b = pysam.AlignedSegment()\par
             b.query_name = query_name1\par
             b.query_sequence = query_sequence_ref\par
             b.flag = 2115\par
             b.reference_id = reference_id\par
             b.reference_start = query_start\par
             b.mapping_quality = 60\par
             cig_temp1 = str(length2) + "M"\par
             b.cigarstring = (cig_temp1)\par
             b.next_reference_id = reference_id\par
             b.next_reference_start = int(location[3]) + 20\par
             b.template_length = length2\par
             quality = "<" * length2\par
             b.query_qualities = pysam.qualitystring_to_array(quality)\par
             b.tags = (("NM", 1), ("RG", "L1"))\par
             outf.write(b)\par
\par
        # for alternate sequence with soft_clipping but correct primer start, good base quality, good mapping quality, sec or suppl alignments\par
          for x in range(5):\par
             soft_clip = 30\par
             match_length = length2 - soft_clip\par
             cig_temp1 = str(match_length) + "M" + str(soft_clip) + "S"\par
             query_name1 = read_names.pop()\par
             b = pysam.AlignedSegment()\par
             b.query_name = query_name1\par
             b.query_sequence = query_sequence_ref\par
             b.flag = 99\par
             b.reference_id = reference_id\par
             b.reference_start = query_start\par
             b.mapping_quality = 60\par
             b.cigarstring = (cig_temp1)\par
             b.next_reference_id = reference_id\par
             b.next_reference_start = int(location[3]) + 20\par
             b.template_length = length2\par
             quality = "<" * length2\par
             b.query_qualities = pysam.qualitystring_to_array(quality)\par
             b.tags = (("NM", 1), ("RG", "L1"))\par
             outf.write(b)\par
\par
}
 