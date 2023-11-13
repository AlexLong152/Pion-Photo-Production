# no arguments ---------------------------- runs input.dat
# argument "test", or "t", just passes ---- runs testInput.dat
# argument "physical" or "p" -------------- runs physicalInput.dat
clear
make clean
make
chmod +x run.twobodyvia2Ndensity

if [ $# -eq 0 ]; then
# if no arguments then run input.dat
./run.twobodyvia2Ndensity input.dat
fi

if [[ "$1" == "test" || "$1" == "t" ]]
then # if "test" is passed then run testInput.dat
./run.twobodyvia2Ndensity testInput.dat
fi

if [[ "$1" == "nosymtest" || "$1" == "nosymt" ]]
then # if "test" is passed then run testInput.dat
./run.twobodyvia2Ndensity nosym-testInput.dat
fi
if [[ "$1" == "p" || "$1" == "physical" ]]
then # if "physical" is passed then run physicalInput.dat
./run.twobodyvia2Ndensity physicalInput.dat
fi
