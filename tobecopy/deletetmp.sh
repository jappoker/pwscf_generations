for i in 400.0 500.0 600.0 700.0 800.0 900.0 1000.0 1100.0 1200.0 1300.0 1400.0 1500.0 ; do
cd scf_$i
rm -rf tmp
sleep 3
cd ..
done
wait
