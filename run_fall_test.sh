n=100
#n=10

rm results.txt

for ((i=0;i<$n;i++))
do
    ./coordconv > out_meas.txt
    tail -n 1 out_meas.txt >> results.txt
done

