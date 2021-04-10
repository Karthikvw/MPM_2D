//C++ headder
#include<memory>
#include<iostream>
#include<string>
#include<vector>
#include<array>
#include<iomanip>
#include<cmath>
#include<functional>
#include<cfloat>
#include<utility>

//Pybind11 headder
#include<pybind11/pybind11.h>
#include<pybind11/iostream.h>
#include<pybind11/numpy.h>
#include<pybind11/stl.h>
#include<pybind11/functional.h>

namespace py = pybind11;

//Header of the class to bind
#include <MPM_Solver.hpp>

void init_PythonMPMSolver(py::module &m)
{
    //Registering Python type
    auto  MPMSolverClass = pybind11::class_<MPM_Solver, std::shared_ptr<MPM_Solver>>(m, "MPM_Solver");

    //Registering the constructor
    MPMSolverClass.def(
        py::init<>()
        ,"Solver class"
        ,py::call_guard<py::scoped_ostream_redirect>()
    );

    //Register member function
    MPMSolverClass.def(
        "Solve"
        ,&MPM_Solver::Solve
        ,py::arg("Grid")
        ,py::arg("Bodies")
        ,py::arg("NoS")
        ,py::arg("dt")
        ,py::arg("alpha")
    );

    //Adding Boundary Conditions
    MPMSolverClass.def(
        "addBC"
        ,[](MPM_Solver &Solver, const std::function<bool(std::array<double,2>)> &Criterion, const double vx, const double vy)
          {
            MPM_Solver::BC newBC(Criterion);
            if (vx!=-DBL_MAX) {newBC.active_vx=true; newBC.value_vx=vx;}
            if (vy!=-DBL_MAX) {newBC.active_vy=true; newBC.value_vy=vy;}
            Solver.BoundaryConditions.push_back(newBC);
            return &Solver;
          }
          ,"Add boundary conditions."
          ,pybind11::arg("SelectionFunction")
          ,pybind11::arg("vx") = -DBL_MAX
          ,pybind11::arg("vy") = -DBL_MAX
          ,pybind11::call_guard<pybind11::scoped_ostream_redirect>()
    );

};
