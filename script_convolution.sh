
# You can run this file in Unix shell to convolute all files that matches the signature
# just go to the directory where those files are located and run 
# $ bash script_convolution.sh
# I have used the forllowing signature as an example

signature=SitePercolationL1__entropy-order_L*
files=`ls $signature`
echo $files

for file in $files:
do
    echo $file
    ./convolution --in $file --with 1,2,3,4 --without 0 --threads 6 --precision 15

done
