echo "Fitness is:"
grep -n "FITNESS" archipelago/G0/done/*
npop=(`ls archipelago/G0/prog | wc -l`)
nbatch=(`squeue -u cns7ae | wc -l`)
nbatch=$((nbatch-2))
echo "*************************"
echo "**** PROG:$npop SQUEUE:$nbatch ****"
echo "*************************"

