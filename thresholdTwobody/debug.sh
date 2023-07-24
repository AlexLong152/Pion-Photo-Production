rm *.tmp
make clean
make
chmod +x run.twobodyvia2Ndensity
gdb --args twobodyvia2Ndensity input.dat
