clear
rm out
gfortran -o  out main.f
#gfortran -o -ffixed-line-length-132 out main.f
./out
