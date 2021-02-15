//C++ header
#include<memory>
#include<iostream>
#include<string>
#include<vector>
#include<array>
#include<iomanip>

//Pybind11 header
#include<pybind11/pybind11.h>
#include<pybind11/iostream.h>
#include<pybind11/numpy.h>
#include<pybind11/stl.h>

namespace py = pybind11;

#include<MPM_Grid.hpp>         //Header of the class to bind

void init_PythonMPMGrid(py::module &m)
{
    //Registering a Python type
    auto  MPMGridClass = py::class_<MPM_Grid, std::shared_ptr<MPM_Grid>>(m, "MPM_Grid");

    //Registering Constructor
    MPMGridClass.def(
        py::init<array<double,2>, array<size_t,2>, array<double,2>>()
        ,"Grid Constructor docstring"
        ,py::arg("L")
        ,py::arg("N")
        ,py::arg("O") = array<double,2> ({0,0})
        ,py::call_guard<py::scoped_ostream_redirect>()
    );

    //Register member properties
    MPMGridClass.def_property_readonly(
        "NoNodes"
        ,[](MPM_Grid &m){return m.NoNodes;}
        ,"Number of nodes"
    );

    //Registering a member function
    MPMGridClass.def(
            "VTKGrid"
            ,&MPM_Grid::VTKGrid
            ,"Creates VTK files"
            ,py::arg("output_file_name")
            ,py::arg("step")
    );
}
