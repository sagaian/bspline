bspline
=======

//Download Git
sudo apt-get install git

//Clone Project. This is your BSpline folder.
git clone https://github.com/sagaian/bspline.git

//Download LAPACK and BLAS
sudo apt-get install liblapack-dev

//Download F2C
sudo apt-get install libf2c2-dev

//Levmar Optimizer
//Download code from  http://www.ics.forth.gr/~lourakis/levmar/
//Extract into a folder of your choice, cd into folder. In terminal, type

make
sudo cp liblevmar.a /usr/local/lib/

//Download gsl libraries
sudo apt-get install gsl

// Run BSpline
// cd into bspline folder. In terminal type

make clean
make
./spline


//Using Git

// get most up to date version
git pull

// add files to commit
git add .

// commit
git commit -m "some comment"

// push
git push
