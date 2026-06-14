#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH -J Q96EB6_CIT
#SBATCH --mem=0
#SBATCH --output=Q96EB6_CIT-%A.out
#SBATCH -D /home/rpearson/research/23ValidationSet/CIT

cd /home/rpearson/research/23ValidationSet/CIT

srun --ntasks=1 --nodes=1 --cpus-per-task=28 \
bash -c "echo Making Q96EB6 directory.; mkdir -p Q96EB6 ; \
cp ./fastas/Q96EB6.fasta ./Q96EB6 ; echo Running Q96EB6 ; \
cd /home/rpearson/research/23ValidationSet/CIT/Q96EB6 ; mv Q96EB6.fasta seq.fasta ; \
perl /home/rpearson/Structure_Prediction_Tools/C-I-TASSER-1.0/I-TASSERmod/runI-TASSER.pl \
-pkgdir /home/rpearson/Structure_Prediction_Tools/C-I-TASSER-1.0 \
-libdir /home/rpearson/Structure_Prediction_Tools/CIT_Lib \
-seqname Q96EB6.fasta -datadir /home/rpearson/research/23ValidationSet/CIT/Q96EB6 \
-outdir /home/rpearson/research/23ValidationSet/CIT/Q96EB6 \
-runstyle parallel -homoflag benchmark -idcut 0.3 -light true -nmodel 5 -hours 5 \
-LBS false -EC false -GO false -java_home /usr -cit true; \
echo Q96EB6 complete."
