#####################################
#                                   #
# Jayadev Joshi, Cleveland Clinic.  #
#                                   #
#####################################


import os, glob, sys
import argparse

parser = argparse.ArgumentParser(description='Description of your program')

parser.add_argument('-I','--InPutDir', help='Enter the path of BASE_DIRECTORY', required=True)
parser.add_argument('-O','--OutPutDir', help='Enter the path for Out_put files', required=True)
parser.add_argument('-P','--DataType', help='if single end SE or if paried end PE', required=False, default='SE')
parser.add_argument('-G','--GenomeDir', help='if single end SE or if paried end PE', required=False)
parser.add_argument('-g','--sjdbGTFfile', help='if single end SE or if paried end PE', required=False)

args = vars(parser.parse_args())

path_to_raw_files = args['InPutDir']
path_to_out_files = args['OutPutDir']
Seq_data_type = args['DataType']
Path_to_Ref_Genome = args['GenomeDir']
Path_to_ref_GTF_file = args['sjdbGTFfile']

files = glob.glob(path_to_raw_files+'/*.fastq')

if not os.path.exists(path_to_out_files+'/Result_files'):
    os.makedirs(path_to_out_files+'/Result_files')

print """
 #############################################################
 #                  Step 1   RUNNING FASTP                   #
 #############################################################
"""

if Seq_data_type == 'SE':

    if not os.path.exists(path_to_out_files+'/Result_files'+'/FASTP_Out_files'):
        os.makedirs(path_to_out_files+'/Result_files'+'/FASTP_Out_files')

    for file in files:

        os.environ["In_file"] = file
        os.environ["Out_file"] = path_to_out_files+'/Result_files'+'/FASTP_Out_files/'+file.split('/')[len(file.split('/'))-1].split('.')[len(file.split('/')[len(file.split('/'))-1].split('.'))-2]+'_Out.fastq'
        os.system('fastp -w 22 -i $In_file -o $Out_file')

else:

    if not os.path.exists(path_to_out_files+'/Result_files'+'/FASTP_Out_files'):
        os.makedirs(path_to_out_files+'/Result_files'+'/FASTP_Out_files')

    R1s = glob.glob(path_to_raw_files+'/*_1.fastq')
    R2s = glob.glob(path_to_raw_files+'/*_2.fastq')

    for i,n in enumerate(R1s):
        n.split('/')[len(n.split('/'))-1].split('_')[len(n.split('/')[len(n.split('/'))-1].split('_'))-2]

        for r in R2s:
            if n.split('/')[len(n.split('/'))-1].split('_')[len(n.split('/')[len(n.split('/'))-1].split('_'))-2] in r:


                os.environ['R1'] = n
                os.environ['R2'] = r

                os.environ['R1_out'] = path_to_out_files+'/Result_files'+'/FASTP_Out_files/'+n.split('/')[len(n.split('/'))-1].strip('.fastq')+'_Out.fastq'
                os.environ['R2_out'] = path_to_out_files+'/Result_files'+'/FASTP_Out_files/'+r.split('/')[len(n.split('/'))-1].strip('.fastq')+'_Out.fastq'

                os.system('fastp -w 22 --in1 $R1 --in2 $R2 --out1 $R1_out --out2 $R2_out')

print """
 #############################################################
 #                  Step 2   RUNNING STAR                    #
 #############################################################
"""


if Seq_data_type == 'SE':


    files = glob.glob(path_to_out_files+'/Result_files/FASTP_Out_files/*.fastq')

    print files
    
    if not os.path.exists(path_to_out_files+'/Result_files'+'/STAR_Out_files'):
        os.makedirs(path_to_out_files+'/Result_files'+'/STAR_Out_files')

    for file in files:

        os.environ["in_file"] = file 
        os.environ["prefix"] = path_to_out_files+'/Result_files'+'/STAR_Out_files/'+file.split('/')[len(file.split('/'))-1].split('_Out')[0]+'_Out'
        os.environ["Path_to_Ref_Genome"] = Path_to_Ref_Genome 
        os.environ["Path_to_ref_GTF_file"] = Path_to_ref_GTF_file

        os.system('STAR  --runMode alignReads --genomeDir $Path_to_Ref_Genome --genomeLoad NoSharedMemory --readFilesIn $in_file --readFilesCommand "zcat -fc" --outStd SAM --runThreadN 18 --outFilterMultimapNmax 10 --outSAMmode Full --outSAMattributes Standard --outSAMstrandField intronMotif --outFileNamePrefix $prefix --outReadsUnmapped Fastx --outFilterScoreMinOverLread 0.9 --outFilterMismatchNoverLmax 0.05 --outFilterMismatchNmax 4 --sjdbGTFfile $Path_to_ref_GTF_file --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate  --runDirPerm All_RWX')
        os.system('rm $in_file')

else:   

    R1_out = glob.glob(path_to_out_files+'/Result_files'+'/FASTP_Out_files/*1_Out.fastq')
    R2_out = glob.glob(path_to_out_files+'/Result_files'+'/FASTP_Out_files/*2_Out.fastq')

    if not os.path.exists(path_to_out_files+'/Result_files'+'/STAR_Out_files'):
        os.makedirs(path_to_out_files+'/Result_files'+'/STAR_Out_files')

    for i,x in enumerate(R1_out):
        if R1_out[i].split('/')[len(R1_out[i].split('/'))-1].split('_')[0] in R2_out[i]:

            os.environ["prefix"] = path_to_out_files+'/Result_files'+'/STAR_Out_files/'+R1_out[i].split('/')[len(R1_out[i].split('/'))-1].split('_')[0]+'_Out'
            os.environ["R1_out"] = R1_out[i]
            os.environ["R2_out"] = R2_out[i]
            os.environ["Path_to_Ref_Genome"] = Path_to_Ref_Genome 
            os.environ["Path_to_ref_GTF_file"] = Path_to_ref_GTF_file

            #print os.environ["prefix"] 
            #print os.environ["R1_out"] 
            #print os.environ["R2_out"] 
            #print os.environ["Path_to_Ref_Genome"]  
            #print os.environ["Path_to_ref_GTF_file"] 

            os.system('STAR  --runMode alignReads --genomeDir $Path_to_Ref_Genome --genomeLoad NoSharedMemory --readFilesIn $R1_out $R2_out --readFilesCommand "zcat -fc" --outStd SAM --runThreadN 18 --outFilterMultimapNmax 10 --outSAMmode Full --outSAMattributes Standard --outSAMstrandField intronMotif --outFileNamePrefix $prefix --outReadsUnmapped Fastx --outFilterScoreMinOverLread 0.9 --outFilterMismatchNoverLmax 0.05 --outFilterMismatchNmax 4 --sjdbGTFfile $Path_to_ref_GTF_file --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate  --runDirPerm All_RWX')


print """
 #############################################################
 #                  Step 3   RUNNING STRINGTIE               #
 #############################################################
"""

if not os.path.exists(path_to_out_files+'/Result_files'+'/STRINGTIE_Out_files'):
    os.makedirs(path_to_out_files+'/Result_files'+'/STRINGTIE_Out_files')


files  = glob.glob(path_to_out_files+'/Result_files'+'/STAR_Out_files/*.bam')

for f in files:

    print f.split('/')[len(f.split('/'))-1].split('_')[0]+'_out.gtf'

    os.environ["Str_out"] = path_to_out_files+'/Result_files'+'/STRINGTIE_Out_files/'+f.split('/')[len(f.split('/'))-1].split('_')[0]+'_out.gtf'
    os.environ["Str_in_file"]  = f
    os.environ["Str_prefix"]  = f.split('/')[len(f.split('/'))-1].split('_')[0]

    os.system('stringtie -p 20 -G $Path_to_ref_GTF_file -o $Str_out -l $Str_prefix $Str_in_file')


print """
 #############################################################
 #                  Step 4   creating merge_file             #
 #############################################################
"""

print "creating 'merge_file_list.txt' file and generating 'stringtie_merged.gtf'"

files = glob.glob(path_to_out_files+'/Result_files'+'/STRINGTIE_Out_files/*.gtf')


merge_file = open(path_to_out_files+'/Result_files/merge_file_list.txt','w')

for f in files:
    merge_file.write(f+'\n')

merge_file.close()


os.environ["merge_list"] = path_to_out_files+'/Result_files/merge_file_list.txt'

os.environ["Str_merge_out"] = path_to_out_files+'/Result_files/stringtie_merged.gtf'

os.system('stringtie --merge -p 20 -G $Path_to_ref_GTF_file  -o $Str_merge_out $merge_list')


print """
 #############################################################
 #                  Step 5   STRINGTIE with Merge Data       #
 #############################################################
"""

files  = glob.glob(path_to_out_files+'/Result_files'+'/STAR_Out_files/*.bam')


if not os.path.exists(path_to_out_files+'/Result_files'+'/Merged_STRINGTIE_Out_files'):
    os.makedirs(path_to_out_files+'/Result_files'+'/Merged_STRINGTIE_Out_files')

for f in  files:

    os.environ["Merge_Str_out"] = path_to_out_files+'/Result_files'+'/Merged_STRINGTIE_Out_files/'+f.split('/')[len(f.split('/'))-1].split('_')[0]+'_out.gtf'
    os.environ["Merge_Str_in_file"]  = f

    Out_file_dir  = f.split('/')[len(f.split('/'))-1].split('_')[0]

    os.environ["Merge_Str_out"] = path_to_out_files+'/Result_files'+'/Merged_STRINGTIE_Out_files/'+Out_file_dir+'/'+f.split('/')[len(f.split('/'))-1].split('_')[0]+'_out.gtf'

    if not os.path.exists(path_to_out_files+'/Result_files'+'/Merged_STRINGTIE_Out_files/'+Out_file_dir):
        os.makedirs(path_to_out_files+'/Result_files'+'/Merged_STRINGTIE_Out_files/'+Out_file_dir)
    
    os.system('stringtie -e -B -p 22 -G $Str_merge_out -o $Merge_Str_out $Merge_Str_in_file') 
