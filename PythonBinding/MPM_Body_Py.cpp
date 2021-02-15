//c++ header
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

//Header of the classes to bind
#include<MPM_Body.hpp>
#include<MPM_STVKSolid.hpp>
#include<MPM_NHSolid.hpp>
#include<MPM_Fluid.hpp>

//Helper function - to view as a numpy array
template<typename T, size_t size>
inline py::array_t<double> numpy_view(std::vector<std::array<T, size>> &Matrix)
{
  return py::array_t<double>(
            { Matrix.size(), size },
            { sizeof(std::array<T, size>), sizeof(double) },
            (double*)Matrix[0].data(),
            py::capsule((double*)Matrix[0].data(), [](void *f){}));
}

void init_PythonMPMBody(py::module &m)
{
    //Registering a Python type
    auto  MPMBodyClass = pybind11::class_<MPM_Body, std::shared_ptr<MPM_Body>>(m, "MPM_Body");

    //Register member properties
    MPMBodyClass.def_property_readonly(
        "NoMP"
        ,[](MPM_Body &m){return m.NoMP;}
        ,"Number of material points"
    );

    //Register member functions
    MPMBodyClass.def(
        "setvelocity"
        ,&MPM_Body::setvelocity
        ,"Set the velocity for the whole body"
        ,py::arg("vx")
        ,py::arg("vy")
        ,py::call_guard<py::scoped_ostream_redirect>()
    );
    
    MPMBodyClass.def(
        "KE"
        ,&MPM_Body::KE
        ,"Calculates the Kinetic energy for the whole body"
        ,py::arg("rho")
        ,py::call_guard<py::scoped_ostream_redirect>()
    );    
        
    MPMBodyClass.def(
        "SE"
        ,&MPM_Body::SE
        ,"Calculates the Strain energy for the whole body"
        ,py::call_guard<py::scoped_ostream_redirect>()
    );

    MPMBodyClass.def(
            "VTKMaterialpoints"
            ,&MPM_Body::VTKMaterialpoints
            ,"Creates VTK files"
            ,py::arg("output_file_name")
            ,py::arg("step")
    );

    //Register a numpy view on an internal data field (to make it modifyable from python)
    MPMBodyClass.def(
        "MPC"
        ,[](MPM_Body &Body)
        {
           return numpy_view(Body.MPC);
        }
        ,"Numpy access to the bodies MPC"
        ,py::call_guard<py::scoped_ostream_redirect>()
    );
};


void init_PythonMPMSTVKSolid(py::module &m)
{
    //Registering a Python type
    auto  MPMSTVKSolidClass = pybind11::class_<MPM_STVKSolid, MPM_Body, std::shared_ptr<MPM_STVKSolid>>(m, "MPM_STVKSolid");

    //Registering Constructor
    MPMSTVKSolidClass.def(
        py::init<double, double, double, vector<array<double,3>> &>()
        ,"Constructor for MPMSolid Class"
        ,py::arg("Emod")
        ,py::arg("Density")
        ,py::arg("Poisson")
        ,py::arg("MPV")
        ,py::call_guard<py::scoped_ostream_redirect>()
    );

};

void init_PythonMPMNHSolid(py::module &m)
{
    //Registering a Python type
    auto  MPMNHSolidClass = pybind11::class_<MPM_NHSolid, MPM_Body, std::shared_ptr<MPM_NHSolid>>(m, "MPM_NHSolid");

    //Registering Constructor
    MPMNHSolidClass.def(
        py::init<double, double, double, vector<array<double,3>> &>()
        ,"Constructor for MPMNHSolid Class"
        ,py::arg("Emod")
        ,py::arg("Density")
        ,py::arg("Poisson")
        ,py::arg("MPV")
        ,py::call_guard<py::scoped_ostream_redirect>()
    );

};

void init_PythonMPMFluid(py::module &m)
{
    //Registering a Python type
    auto  MPMFluidClass = pybind11::class_<MPM_Fluid, MPM_Body, std::shared_ptr<MPM_Fluid>>(m, "MPM_Fluid");

    //Registering Constructor
    MPMFluidClass.def(
        py::init<double, double, double, vector<array<double,3>> &, double>()
        ,"Constructor for MPMFluid Class"
        ,py::arg("visc")
        ,py::arg("Density")
        ,py::arg("bulk_mod")
        ,py::arg("MPV")
        ,py::arg("gamma") = 7
        ,py::call_guard<py::scoped_ostream_redirect>()
    );

};
