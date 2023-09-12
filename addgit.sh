cp ../get_file.py .
git add get_file.py
git add thresholdOnebody/*
git add thresholdTwobody/*
git add .testing/mySigma/main.f
git add .testing/mySigma/run.sh
git add .testing/sigma_plus_minus_sigma/main.f
git add .testing/sigma_plus_minus_sigma/run.sh
git add common-densities/*

git add common-densities/density-modules/*
git add README.md
git commit 

rm get_file.py

echo "To push to github run: git push -u origin main"
echo "or: git push -uf origin main"
echo "and log in with credentials"

