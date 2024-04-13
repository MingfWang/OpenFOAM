/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-----------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    TrajectoryToVTK

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include <string>

#include <vtk-9.3/vtkPolyData.h>
#include <vtk-9.3/vtkSmartPointer.h>
#include <vtk-9.3/vtkPoints.h>
#include <vtk-9.3/vtkXMLPolyDataWriter.h>
#include <vtk-9.3/vtkPolyLine.h>
#include <vtk-9.3/vtkCellArray.h>
#include <vtk-9.3/vtkIntArray.h>
#include <vtk-9.3/vtkPointData.h>
#include <vtk-9.3/vtkCellData.h>

void points2vtk(DynamicList<DynamicList<vector>> &path, fvMesh &mesh, fileName &VTKDir);

// ************************************************************************* //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // 粒子一般在入口处出发，这里就是指定了粒子出发的边界
    IOdictionary trajectoryDict(
            IOobject
            (
                    "trajectoryDict",
                    runTime.system(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                    )
            );
    word boundaryName; // 给定的边界名称
    trajectoryDict.lookup("boundaryName") >> boundaryName;
    label indexDensity; // 选定粒子的间隔
    trajectoryDict.lookup("indexDensity") >> indexDensity;
    scalar maxCoordinate; // 给定计算域流向的长度，用于筛选中断的情况。中断是因为碰到了noSlip边界条件。
    trajectoryDict.lookup("maxCoordinate") >> maxCoordinate;
    word flowDirection;
    trajectoryDict.lookup("flowDirection") >> flowDirection;

    // 创建一个文件夹用于存放VTK文件
    fileName VTKDir = mesh.time().path() / "Trajectory"; // time()函数返回的网格注册到哪个时间
    mkDir(VTKDir);

    Info << VTKDir << endl;
    // 创建文件指针，写出每隔颗粒的移动距离
    autoPtr<OFstream> distancePtr;
    distancePtr.reset(new OFstream(VTKDir / "distance.csv"));

    // 写出位置文件
    autoPtr<OFstream> positionPtr;
    positionPtr.reset(new OFstream(VTKDir / "position.csv"));

    instantList times = runTime.times(); // 得到所有时间
    runTime.setTime(times.last(), 0); // 设置当前时可为最后一个时间点
    // 主要是为了读取最后一个时间步长的速度
    Info << nl << "使用的数据时刻为： " << runTime.timeName() << endl;

    // 读取速度场
    Info << nl << "读取速度场：" << endl;
    volVectorField U
            (
                    IOobject
                    (
                            "U",
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                            ),
                    mesh
                    );

    // 创建一个一维vector数组，用于存放每个起点
    DynamicList<vector> initialParticlePositions; // 用于存放每个颗粒的起始位置
    DynamicList<DynamicList<vector>> finalParticlePositions; // 用于存放每个颗粒的轨迹坐标
    DynamicList<scalar> distance; // 用于存放每个颗粒的运行长度

    label index = mesh.boundaryMesh().findIndex(boundaryName);
    // index是inlet边界的索引号
    // 当没有要查找的边界时，会返回-1
    if ( index == -1)
    {
        FatalErrorInFunction << "没有该边界！" << abort(FatalError);
    }

    // 选择边界面上的若干点
    forAll(mesh.boundaryMesh()[index], iP)
    {
        // 每隔indexDensity选择一个起始点，但是不包括0，因为第0个面网格位置太偏僻
        if (iP % indexDensity == 0 && iP != 0) {
            DynamicList<vector> temp;
            temp.append(mesh.boundaryMesh()[index].faceCentres()[iP]);
            initialParticlePositions.append(temp);
        }
    }

    forAll(initialParticlePositions, index)
    {
        Info << "开始追踪第" << index << "个颗粒！" << endl;
        // 从initialParticlePositions中起始颗粒的位置坐标
        vector initialParticleInner = initialParticlePositions[index];

        // 准备存储颗粒的动态列表Dynamic List
        DynamicList<vector> particlePositionsInner;
        scalar timeTaken(0.0);
        scalar distanceTraveled(0.0);

        // 将初始位置的坐标append到动态列表dynamicList中
        particlePositionsInner.append(initialParticleInner);

        // 找到颗粒所在的网格编号，注意是网格不是边界面编号
        label cellID = mesh.findCell(initialParticleInner);
        Info << nl << "初始颗粒位于网格 " << cellID << "中！" << endl;

        vector currentPosition = initialParticleInner; // 当前位置
        vector newPosition(0, 0, 0); // 移动到新的位置
        label iterCount(1); // 计数器
        const label maxIters(10000000); // 颗粒的最大移动次数

        Info << "开始追踪颗粒" << endl;
        while (cellID != -1)
        {
            vector velocity = U[cellID]; // 得到网格的速度
            scalar charLen = Foam::cbrt(mesh.V()[cellID]); // 得到网格的特征尺寸
            scalar dt = charLen / mag(velocity); // 得到特征时间步长
            newPosition = currentPosition + velocity * dt; // 计算得到新的颗粒位置

            scalar dist = mag(newPosition - currentPosition); // 得到颗粒在当前时间步的位移

            // 输出信息
            Info << nl << "Lagrangian time step: " << iterCount << nl
                 << tab << "current position = " << currentPosition << nl
                 << tab << "new position = " << newPosition << nl
                 << tab << "local distance traveled = " << dist << nl
                 << tab << "local time taken = " << dt << nl
                 << tab << "currently cell no." << cellID << endl;

            // 将数据添加到列表中
            distanceTraveled += dist;
            timeTaken += dt;
            particlePositionsInner.append(newPosition);
            currentPosition = newPosition;
            iterCount++;

            // 当前所在的网格
            cellID = mesh.findCell(currentPosition);

            if (iterCount > maxIters)
            {
                FatalErrorInFunction << "达到最大的跟踪时间步数！" << abort(FatalError);
            }

            if (cellID == -1) // 颗粒离开计算区域
            {
                Info << nl
                     << "颗粒离开计算区域！" << nl
                     << "颗粒运行的总距离： " << distanceTraveled << nl
                     << "颗粒停留时间： " << timeTaken << endl;
            }

        }
        Info << "第" << index << "个颗粒追踪结束！" << nl << endl;

        // 增加筛选，如果大于长度就写出否则不写出，tortuosity大于1
        if ((flowDirection == "x" && currentPosition.x() > maxCoordinate) ||
            (flowDirection == "y" && currentPosition.y() > maxCoordinate) ||
            (flowDirection == "z" && currentPosition.z() > maxCoordinate))
        {
            finalParticlePositions.append(particlePositionsInner); // 将每个颗粒的轨迹存放于finalParticlePositions这个二维动态列表中
            // 写出每个颗粒的移动总距离
            distancePtr() << index << ", " << distanceTraveled << endl;
            distance.append(distanceTraveled); //用于存放所有选中颗粒的运行长度
        }
    }

    forAll(finalParticlePositions, index)
    {
        Info << "开始写入第" << index << "个颗粒！" << endl;
        forAll(finalParticlePositions[index], indexInner)
        {
            positionPtr() << index << ", "
                          << finalParticlePositions[index][indexInner].x() << ", "
                          << finalParticlePositions[index][indexInner].y() << ", "
                          << finalParticlePositions[index][indexInner].z() << endl;
        }
        positionPtr() << endl;
        Info << "第" << index << "个颗粒写入完成！" << endl;
    }
    // 写出vtk文件中
    Info << endl;
    points2vtk(finalParticlePositions, mesh, VTKDir);

    Info << nl;
    runTime.printExecutionTime(Info);
    Info << "End" << endl;

    return 0;
}

void points2vtk(DynamicList<DynamicList<vector>> &path, fvMesh &mesh, fileName &VTKDir)
{
    Info << "开始将颗粒的路径写入到VTK文件中！" << endl;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    // vtkCellArray的主要作用是管理拓扑关系

    // 创建存储线编号的标量数组
    vtkSmartPointer<vtkIntArray> lineIds = vtkSmartPointer<vtkIntArray>::New();
    lineIds->SetName("LineIds");
    lineIds->SetNumberOfComponents(1);
    lineIds->SetNumberOfTuples(path.size());


    forAll(path, index)
    {
        vtkSmartPointer<vtkPolyLine> polyline = vtkSmartPointer<vtkPolyLine>::New();
        polyline->GetPointIds()->SetNumberOfIds(path[index].size()); // 设置本条线有多少个点
        // 将点写入到线中
        forAll(path[index], indexInner)
        {
            // pointId是点的全局编号
            vtkIdType pointId = points->InsertNextPoint(path[index][indexInner].x(),
                                                        path[index][indexInner].y(),
                                                        path[index][indexInner].z());
            polyline->GetPointIds()->SetId(indexInner, pointId);
            // SetId中的第一个参数是点在本线中编号，第二个参数是点在全局中的编号
        }
        lines->InsertNextCell(polyline);

        // 为每条线的标量数组添加线编号
        lineIds->SetValue(index, index);
    }

    // vtkPolyData的主要作用是存储所有数据
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetLines(lines);
    polydata->GetCellData()->AddArray(lineIds); // 将标量数组添加到Cell Data中

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    fileName VTKName = "/trajectory.vtp";
    fileName fullDir = VTKDir + VTKName; // 也可以用“+”将两个字符串链接起来

    writer->SetFileName(fullDir.c_str()); //
    writer->SetInputData(polydata);
    writer->Write();

    Info << "所有颗粒的路径写入完成！" << endl;
}