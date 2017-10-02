#! /bin/csh
#$ -V
#$ -S /bin/csh
#$ -o /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -e /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -cwd
#$ -t 1-74
#$ -l h_vmem=8G
#$ -tc 50
#$ -l long
set samples=(e772ff12-17b1-473e-b581-d1303398e09b 2f8baf51-8d2a-4d62-a92a-a62f1fd5cd77 776c994e-b68a-4f52-90c8-02c40e0e80b2 9d5005b5-e663-483d-95d4-be5a8f33d9a6 0ee3d529-a972-4ddb-9fa1-b8c8f0b26684 8e16aad5-78ec-4246-8013-0b760f0aea2b ce92e89d-ef33-45cf-a7a7-0e509f8c6c29 18c27e7a-2a5f-438e-a4c0-39140e3d95ed fd51273d-6b69-4e29-a96e-53839d3c6f77 78b46f0d-4e59-4ea4-981c-0a3f6c59019f 551032b2-6327-428c-83f8-55d312a000c6 ceffbb44-6e4f-44f8-9afa-407b9bc1a35c 26fa0e90-90f6-4643-9ab7-8b6ffda9b635 b82ecb1c-6c7b-49f6-86b4-5362ff6629ea 6679198d-74a0-4d98-9cf7-ad76581eca66 9713549a-2f0b-47ef-a4f5-2f0b8959e9e0 c8894a92-26e1-4a44-9311-b0e4300dc0e0 4ccd57c5-bbac-4144-b159-9f80b99cce61 c2b24763-5491-4741-b38a-267e4d22af0d cd8945cc-c07c-4ddd-9f15-d63349ceb0da 73c224ae-95b0-48a9-b583-18104301aca5 0a993224-ab18-4013-9655-8128980305bb d611a9d9-9c71-4a32-beb9-cc8058a72f4b 6d7bc15b-c2e9-42d9-9fd1-004ad3a13e17 ac52fde5-c591-43df-9294-3b97fdba5a7a a7f1561a-dada-4d30-a4a8-76318de080b9 25f7e41a-a33c-43d9-a242-a80c921bb53b 21fc93b7-e01a-4942-ba6b-c9a5028c4e60 6f9f4fd3-0812-47ce-9236-4457079b3d03 1a8e4541-eaec-48f0-b570-290e75f57323 f5529875-cf5f-4cde-80de-3ebb6191dfb6 02f11ff5-366b-46a4-97ef-a592b1414549 31c3d889-9960-496a-8eda-853b998bba5f 51b1ee75-a4de-48a9-aea2-4f46223d4dfd b65d5780-cec1-4566-b421-aa536c01482f f8ac0b85-2b14-4145-9290-c5d9e9fa94b8 f7ba182b-fc71-414d-99b7-ebd116596020 10dca98f-036f-4e10-8b79-e1c24748acfb da8bf409-a572-46b1-ae2c-6501dfe6292e 7b9357e2-13eb-46a0-a39f-c8e16922de4d 332505bd-b865-4f3c-a56d-baf01e5ff95c 1afcb3a8-a593-48ee-926c-1ce2d0c5af90 60da9169-68be-4cd8-8b9e-b747ecfdacf4 318667b4-2ca3-4f83-a62e-ef51ccbb6fd8 95fde270-407c-4b76-a670-3a6fcf693b0c d54da8bc-6655-42a4-b3b2-1ef4ca28db9e f7b154ff-61c6-4faa-9a89-2c083626c84d 662dfb77-cbc6-4128-aef9-1255b5b3a2a0 c18b05e7-584e-4788-a9b7-ba08292aa433 9d4193e0-fa0a-41cd-98d9-7be5118b9420 47ad98e1-c9a4-4ebc-be0d-50c4714cf4b3 1571b866-6da0-48a3-875c-21cbaef4e7cb 946ebe91-b00c-4c96-be31-2ef7adc8566b e19d592f-1825-48d9-bfb8-4534f74f7b6d ec5a5c1d-0e91-4459-95fa-f77a6bcb4e9c b0c5a9d5-c06d-4b0e-91d0-22134c95f642 e527fe1d-ccb6-40e2-b450-1b167254e615 a5aaa6f7-9e10-49df-a017-1ee4fbf3837a 657515ff-92cc-42f2-af9d-f1ad55ad4193 bea9fc86-b591-43c8-897e-d96fa45ff977 fab18fd0-f436-4a4a-9f84-838b12c04ec4 05c2820e-fce9-4b90-8469-d28e3a1ee54f 660b5854-2c6b-4724-83b1-aa74b9e3af8f b4ea2516-aa5f-45ae-b8be-23a2e8155ab5 901f8907-4330-4b96-bb83-6b935a9e8e37 35fe03dd-d7f1-46ee-b9d1-ff0136154bef 4a073fb2-25a6-4916-ab68-9bbac1370ded cf7bf0f0-7ee2-4502-b3d2-02cc543b8f86 b5419f58-7833-4008-b6e6-9480b4f96593 075061d7-4c2f-4a19-a9fd-219d2ffc882b 70d48af1-e1c5-4cfb-85ae-345775c0b641 043b3105-e669-4ebd-84d2-98818e7604a8 2bb32af7-248d-49dd-9f7d-442a5124c3f3 3d843b5e-d6aa-428e-bbe3-db1a92dd353c)
set barcodes=(TCGA-DV-5576 TCGA-13-0891 TCGA-AG-3602 TCGA-AA-3877 TCGA-AG-3574 TCGA-60-2721 TCGA-60-2720 TCGA-10-0926 TCGA-13-1495 TCGA-21-1078 TCGA-AG-3600 TCGA-13-0762 TCGA-13-0757 TCGA-13-0906 TCGA-13-1485 TCGA-13-0727 TCGA-13-0921 TCGA-13-0730 TCGA-DV-5574 TCGA-13-0766 TCGA-17-Z035 TCGA-DV-5573 TCGA-13-0714 TCGA-AA-3870 TCGA-13-0924 TCGA-13-0764 TCGA-AA-3930 TCGA-04-1332 TCGA-DD-A116 TCGA-DV-5575 TCGA-AA-3872 TCGA-13-0911 TCGA-AA-3869 TCGA-DD-A1EE TCGA-AA-3562 TCGA-17-Z010 TCGA-13-0916 TCGA-A8-A0A4 TCGA-13-0768 TCGA-13-0794 TCGA-13-0904 TCGA-24-1467 TCGA-AA-3939 TCGA-DV-5569 TCGA-13-0726 TCGA-17-Z014 TCGA-CW-5589 TCGA-13-0797 TCGA-85-8351 TCGA-13-0720 TCGA-06-0649 TCGA-A2-A0SV TCGA-AG-3587 TCGA-13-0883 TCGA-13-0803 TCGA-25-1312 TCGA-17-Z000 TCGA-09-0364 TCGA-13-0919 TCGA-D1-A0ZU TCGA-10-0928 TCGA-CW-5588 TCGA-13-0889 TCGA-AA-3941 TCGA-17-Z018 TCGA-AG-3598 TCGA-13-0758 TCGA-18-4086 TCGA-AB-2864 TCGA-24-0966 TCGA-10-0937 TCGA-DD-A1EA TCGA-13-0717 TCGA-AA-3548)
set outs=(/nrnb/users/ramarty/TCGA/exomes/TCGA-DV-5576 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0891 /nrnb/users/ramarty/TCGA/exomes/TCGA-AG-3602 /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3877 /nrnb/users/ramarty/TCGA/exomes/TCGA-AG-3574 /nrnb/users/ramarty/TCGA/exomes/TCGA-60-2721 /nrnb/users/ramarty/TCGA/exomes/TCGA-60-2720 /nrnb/users/ramarty/TCGA/exomes/TCGA-10-0926 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-1495 /nrnb/users/ramarty/TCGA/exomes/TCGA-21-1078 /nrnb/users/ramarty/TCGA/exomes/TCGA-AG-3600 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0762 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0757 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0906 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-1485 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0727 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0921 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0730 /nrnb/users/ramarty/TCGA/exomes/TCGA-DV-5574 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0766 /nrnb/users/ramarty/TCGA/exomes/TCGA-17-Z035 /nrnb/users/ramarty/TCGA/exomes/TCGA-DV-5573 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0714 /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3870 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0924 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0764 /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3930 /nrnb/users/ramarty/TCGA/exomes/TCGA-04-1332 /nrnb/users/ramarty/TCGA/exomes/TCGA-DD-A116 /nrnb/users/ramarty/TCGA/exomes/TCGA-DV-5575 /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3872 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0911 /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3869 /nrnb/users/ramarty/TCGA/exomes/TCGA-DD-A1EE /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3562 /nrnb/users/ramarty/TCGA/exomes/TCGA-17-Z010 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0916 /nrnb/users/ramarty/TCGA/exomes/TCGA-A8-A0A4 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0768 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0794 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0904 /nrnb/users/ramarty/TCGA/exomes/TCGA-24-1467 /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3939 /nrnb/users/ramarty/TCGA/exomes/TCGA-DV-5569 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0726 /nrnb/users/ramarty/TCGA/exomes/TCGA-17-Z014 /nrnb/users/ramarty/TCGA/exomes/TCGA-CW-5589 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0797 /nrnb/users/ramarty/TCGA/exomes/TCGA-85-8351 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0720 /nrnb/users/ramarty/TCGA/exomes/TCGA-06-0649 /nrnb/users/ramarty/TCGA/exomes/TCGA-A2-A0SV /nrnb/users/ramarty/TCGA/exomes/TCGA-AG-3587 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0883 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0803 /nrnb/users/ramarty/TCGA/exomes/TCGA-25-1312 /nrnb/users/ramarty/TCGA/exomes/TCGA-17-Z000 /nrnb/users/ramarty/TCGA/exomes/TCGA-09-0364 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0919 /nrnb/users/ramarty/TCGA/exomes/TCGA-D1-A0ZU /nrnb/users/ramarty/TCGA/exomes/TCGA-10-0928 /nrnb/users/ramarty/TCGA/exomes/TCGA-CW-5588 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0889 /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3941 /nrnb/users/ramarty/TCGA/exomes/TCGA-17-Z018 /nrnb/users/ramarty/TCGA/exomes/TCGA-AG-3598 /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0758 /nrnb/users/ramarty/TCGA/exomes/TCGA-18-4086 /nrnb/users/ramarty/TCGA/exomes/TCGA-AB-2864 /nrnb/users/ramarty/TCGA/exomes/TCGA-24-0966 /nrnb/users/ramarty/TCGA/exomes/TCGA-10-0937 /nrnb/users/ramarty/TCGA/exomes/TCGA-DD-A1EA /nrnb/users/ramarty/TCGA/exomes/TCGA-13-0717 /nrnb/users/ramarty/TCGA/exomes/TCGA-AA-3548)

set sample=$samples[$SGE_TASK_ID]
set barcode=$barcodes[$SGE_TASK_ID]
set out=$outs[$SGE_TASK_ID]

date
hostname
mkdir $out
mkdir $out/hlaHD
mkdir /tmp/$barcode

bash /cellar/users/ramarty/Projects/kir/KIR_development/data_gathering/bin/GDC.exome.sh $sample /tmp/$barcode/full_exome.bam
echo Bam downloaded.

/cellar/users/ramarty/programs/samtools-1.3.1/samtools index /tmp/$barcode/full_exome.bam
echo index.

/cellar/users/ramarty/programs/samtools-1.3.1/samtools view -b /tmp/$barcode/full_exome.bam "chr6:22361036-47685581" > $out/chr6.bam

/cellar/users/ramarty/programs/samtools-1.3.1/samtools view -b /tmp/$barcode/full_exome.bam "chr19:54051467-55502666" > $out/chr19.bam

/cellar/users/ramarty/programs/samtools-1.3.1/samtools view -b -f 4 /tmp/$barcode/full_exome.bam > $out/unmapped.bam
echo sliced.

/cellar/users/ramarty/programs/samtools-1.3.1/samtools merge /tmp/$barcode/merged.bam $out/chr6.bam $out/unmapped.bam
echo merged.

python /cellar/users/ramarty/Projects/kir/KIR_development/data_gathering/bin/convert_to_fastq2.py /tmp/$barcode/merged.bam /tmp/$barcode/merged_sorted /tmp/$barcode/merged_1.fastq /tmp/$barcode/merged_2.fastq cellar
echo Fastq stripped.

hlahd.sh -t 8 -m 70 -f ~/programs/hlahd.1.0.0/freq_data/ /tmp/$barcode/merged_1.fastq /tmp/$barcode/merged_2.fastq ~/programs/hlahd.1.0.0/HLA_gene.split.txt ~/programs/hlahd.1.0.0/dictionary/ sampleID $out/hlaHD
echo HLA-HD completed.

rm $out/merged*
rm $out/full_exome*
rm -r /tmp/$barcode/*
date
