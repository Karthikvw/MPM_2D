 # Packages for Linux/WSL

The document provides the necessary packages to run the API in a Linux/WSL environment along with installation guidlines

## Anaconda

Install anaconda by copying the download link from the webportal and follow the below instructions

>wget "link"

>chmod +x ***.sh

>***.sh

>conda update -all

## gcc and pip

>sudo apt install python3-pip

## cmake

>sudo apt install cmake

## [pybind11](https://pybind11.readthedocs.io/en/stable/basics.html)

>pip install pybind11

## git

>sudo apt install git-all

# Setting up git repository

## Cloning the repository:

>git clone "***repository link***"  

## Creating and changing to build directory

>mkdir build  

>cd build

## Generating the executables

>cmake ..  

>make  

>make install

