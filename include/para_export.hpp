#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <iomanip>

inline void vtk_export(
    std::string output_file_name, 
    size_t number,
    std::vector<std::array<double,2>> &Points,
    std::vector<std::array<size_t,4>> &Cells,
    std::vector<std::string> &ScalarPointDataNames,
    std::vector< std::vector<double> > &ScalarPointData,
    std::vector<std::string> &VectorPointDataNames,
    std::vector< std::vector<std::array<double,2>> > &VectorPointData
)
{
    // input sanity checks
    if (ScalarPointDataNames.size() != ScalarPointData.size())
        throw std::runtime_error("Error! Number of scalar data names does not match provided data.");
    if (VectorPointDataNames.size() != VectorPointData.size())
        throw std::runtime_error("Error! Number of vector data names does not match provided data.");
    for (size_t i=0; i<ScalarPointDataNames.size(); i++) if (ScalarPointData[i].size() != Points.size())
        throw std::runtime_error("Error! Number of provided "+ScalarPointDataNames[i]+" does not match number of nodes.");
    for (size_t i=0; i<VectorPointDataNames.size(); i++) if (VectorPointData[i].size() != Points.size())
        throw std::runtime_error("Error! Number of provided "+VectorPointDataNames[i]+" does not match number of nodes.");

    // Prepar output name
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(6) << number;
    std::ofstream VTKOut(output_file_name +oss.str() +".vtk");

    // write header
    VTKOut << "# vtk DataFile Version 2.0" << std::endl;
    VTKOut << "VTKOUTPUT FROM MPM2D" << std::endl;
    VTKOut << "ASCII" << std::endl;
    VTKOut << "DATASET POLYDATA" << std::endl;
    
    // write points
    VTKOut << "POINTS " << Points.size() <<" float" << std::endl;
    for (auto &pt : Points)
        VTKOut << pt[0] << " " << pt[1] << " " << 0 << std::endl;
    VTKOut << std::endl;

    // write cells
    VTKOut << "POLYGONS " << Cells.size() << " " << 5*Cells.size() << std::endl;
    for (auto &cell : Cells)
        VTKOut << "4 " << cell[0] << " " << cell[1] << " " << cell[2] << " " << cell[3] << std::endl;
    VTKOut << std::endl;

    // write point data
    VTKOut << "POINT_DATA " << Points.size() << std::endl;
    for (size_t i=0; i<ScalarPointDataNames.size(); i++)
    {
        VTKOut << "SCALARS " << ScalarPointDataNames[i] <<" float 1" << std::endl;
        VTKOut << "LOOKUP_TABLE default " << std::endl;
        for (auto &d : ScalarPointData[i]) VTKOut << d << std::endl;
        VTKOut << std::endl;
    }
    for (size_t i=0; i<VectorPointDataNames.size(); i++)
    {
        VTKOut << "VECTORS " << VectorPointDataNames[i] <<" float" << std::endl;
        for (auto &d : VectorPointData[i]) VTKOut << d[0] << " " << d[1] << " " << 0 << std::endl;
        VTKOut << std::endl;
    }


    VTKOut.close();

}

// overload to have only points which will appear in paraview 
// as cells
inline void vtk_export(
    std::string output_file_name, 
    size_t number,
    std::vector<std::array<double,2>> &Points,
    std::vector<std::string> &ScalarPointDataNames,
    std::vector< std::vector<double> > &ScalarPointData,
    std::vector<std::string> &VectorPointDataNames,
    std::vector< std::vector<std::array<double,2>> > &VectorPointData
)
{
    // input sanity checks
    if (ScalarPointDataNames.size() != ScalarPointData.size())
        throw std::runtime_error("Error! Number of scalar data names does not match provided data.");
    if (VectorPointDataNames.size() != VectorPointData.size())
        throw std::runtime_error("Error! Number of vector data names does not match provided data.");
    for (size_t i=0; i<ScalarPointDataNames.size(); i++) if (ScalarPointData[i].size() != Points.size())
        throw std::runtime_error("Error! Number of provided "+ScalarPointDataNames[i]+" does not match number of nodes.");
    for (size_t i=0; i<VectorPointDataNames.size(); i++) if (VectorPointData[i].size() != Points.size())
        throw std::runtime_error("Error! Number of provided "+VectorPointDataNames[i]+" does not match number of nodes.");

    // Prepar output name
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(6) << number;
    std::ofstream VTKOut(output_file_name +oss.str() +".vtk");

    // write headder
    VTKOut << "# vtk DataFile Version 2.0" << std::endl;
    VTKOut << "VTKOUTPUT FROM MPM2D" << std::endl;
    VTKOut << "ASCII" << std::endl;
    VTKOut << "DATASET POLYDATA" << std::endl;
    
    // write points
    VTKOut << "POINTS " << Points.size() <<" float" << std::endl;
    for (auto &pt : Points)
        VTKOut << pt[0] << " " << pt[1] << " " << 0 << std::endl;
    VTKOut << std::endl;

    // write cells
    VTKOut << "VERTICES " << Points.size() << " " << 2*Points.size() << std::endl;
    for (size_t n=0; n<Points.size(); n++)
        VTKOut << "1 " << n << std::endl;
    VTKOut << std::endl;

    // write point data
    VTKOut << "CELL_DATA " << Points.size() << std::endl;
    for (size_t i=0; i<ScalarPointDataNames.size(); i++)
    {
        VTKOut << "SCALARS " << ScalarPointDataNames[i] <<" float 1" << std::endl;
        VTKOut << "LOOKUP_TABLE default " << std::endl;
        for (auto &d : ScalarPointData[i]) VTKOut << d << std::endl;
        VTKOut << std::endl;
    }
    for (size_t i=0; i<VectorPointDataNames.size(); i++)
    {
        VTKOut << "VECTORS " << VectorPointDataNames[i] <<" float" << std::endl;
        for (auto &d : VectorPointData[i]) VTKOut << d[0] << " " << d[1] << " " << 0 << std::endl;
        VTKOut << std::endl;
    }


    VTKOut.close();

}