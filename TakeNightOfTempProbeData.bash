port=/dev/ttyACM0
int_time=60
var_time=10
n_minutes=660

save_freq=20

total_time=$(($n_minutes*60))
n_exps=$(($total_time/$int_time))
start=1

for (( c=$start; c<=$n_exps; c++ ))
do
  save_all=$(($c%$save_freq))
  if [ "$save_all" == 1 ]
  then
      save_all=1
  else
      save_all=0
  fi
  echo int_time var_time save_all
  echo $int_time $var_time $save_all
	python3 ReadDataFromJimArduinoCircuit.py $port $int_time $var_time $save_all
  sleep 0.1
done
