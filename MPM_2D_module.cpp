//Include C++ header
#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <iomanip>

//Include the Pybind11 header
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

//Include cpp files
#include <PythonBinding/MPM_Body_Py.cpp>
#include <PythonBinding/MPM_Grid_Py.cpp>
#include <PythonBinding/MPM_Solver_Py.cpp>

namespace py = pybind11;

//Declaration for bindings defined in another library
void init_PythonMPMGrid(py::module &m);
void init_PythonMPMBody(py::module &m);
void init_PythonMPMSTVKSolid(py::module &m);
void init_PythonMPMNHSolid(py::module &m);
void init_PythonMPMFluid(py::module &m);
void init_PythonMPMSolver(py::module &m);

//Generating pybind module
PYBIND11_MODULE(MPM_2D, m) 
{
    m.doc() = "The MPM2D Module";

    init_PythonMPMGrid(m);          //Registering the MPM_Grid class to Python module
    init_PythonMPMBody(m);          //Registering the MPM_Body class to Python module
    init_PythonMPMSTVKSolid(m);         //Registering the MPM_Solid class to Python module
    init_PythonMPMNHSolid(m);
    init_PythonMPMFluid(m);         //Registering the MPM_Fluid class to Python module
    init_PythonMPMSolver(m);        //Registering the MPM_Solver class to Python module   

}
