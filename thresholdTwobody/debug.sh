rm *.tmp
clear;
make clean
cp Makefile makedebug
#replace -O optimizeation with no optimization at all for better debugging
sed -i 's/FFLAGS= -O -/FFLAGS= -O0 -/' makedebug
make -f makedebug
chmod +x run.twobodyvia2Ndensity
echo ""
echo "Debuging with gdb debugger"
echo "------------------------------------------------------------------"
echo "Set breakpoints with:------------------b filename.f:linenumber"
echo "Run code until breakpoint with:--------r or run"
echo "Print with:----------------------------p:variablename"
echo "Go to next line in file with:----------n"
echo "Step into function with:---------------s or step"
echo "https://cs.baylor.edu/~donahoo/tools/gdb/tutorial.html"
echo ""
gdb --args twobodyvia2Ndensity input.dat
