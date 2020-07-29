#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPointPicker.h>
#include <vtkSphereSource.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkConeSource.h>
#include <vtkProperty.h>
#include <vtkPLYWriter.h>
#include"rendercone.h"
#include <stdio.h>
#include<iostream>
#include<sstream>
#include"pch.h"
//#include"algorithm.h"
#include "optimal_nonrigid_icp.h"
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkRenderingOpenGL2);

#include <vtkNew.h>
#include <vtkConeSource.h>
#include <vtkRenderWindow.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballActor.h>

//#include <ImportLibVTK8_2_0.h>
//#include <AutoInitVTKModules.h>

int main()
{
	// 创建圆锥数据
	vtkNew<vtkConeSource> cone;
	cone->SetHeight(20);         //设置圆锥的高度和半径  
	cone->SetRadius(20);
	cone->SetResolution(100);    //设置分辨率，值越大，越趋近于圆锥体
	//创建对应的mapper 
	vtkNew<vtkPolyDataMapper> coneMapper;
	coneMapper->SetInputConnection(cone->GetOutputPort());
	// 创建对应的Actor
	vtkNew<vtkActor> coneActor;
	coneActor->SetMapper(coneMapper);
	// 针对每个actor创建一个vtkRender，添加actor 
	vtkNew<vtkRenderer> coneRenderer;
	coneRenderer->AddActor(coneActor);
	//下面的代码就是设置视口的，将屏幕分为4个视口  
	//前两个参数是视口左下角点的坐标(xmin,ymin)，后两个参数是右上角的坐标(xmax,ymax)  
	coneRenderer->SetViewport(0, 0.5, 0.5, 1);
	coneRenderer->SetBackground(0.1, 0.2, 0.4);

	// 创建立方体
	vtkNew<vtkCubeSource> cube;
	cube->SetXLength(40);     //设置立方体的长宽高  
	cube->SetYLength(50);
	cube->SetZLength(60);
	vtkNew<vtkPolyDataMapper> cubeMapper;
	cubeMapper->SetInputConnection(cube->GetOutputPort());
	vtkNew<vtkActor> cubeActor;
	cubeActor->SetMapper(cubeMapper);
	vtkNew<vtkRenderer> cubeRenderer;
	cubeRenderer->AddActor(cubeActor);
	cubeRenderer->SetViewport(0.5, 0.5, 1, 1);
	cubeRenderer->SetBackground(0, 1, 0);

	//创建圆柱体  
	vtkNew<vtkCylinderSource> cylinder;
	vtkNew<vtkPolyDataMapper> cylinderMapper;
	cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
	vtkNew<vtkActor> cylinderActor;
	cylinderActor->SetMapper(cylinderMapper);
	vtkNew<vtkRenderer> cylinderRenderer;
	cylinderRenderer->AddActor(cylinderActor);
	cylinderRenderer->SetViewport(0, 0, 0.5, 0.5);

	//创建平面  
	vtkNew<vtkPlaneSource> plane;
	vtkNew<vtkPolyDataMapper> planeMapper;
	planeMapper->SetInputConnection(plane->GetOutputPort());
	vtkNew<vtkActor> planeActor;
	planeActor->SetMapper(planeMapper);
	vtkNew<vtkRenderer> planeRenderer;
	planeRenderer->AddActor(planeActor);
	planeRenderer->SetViewport(0.5, 0, 1, 0.5);
	planeRenderer->SetBackground(0, 0, 1);

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(coneRenderer);
	renderWindow->AddRenderer(cubeRenderer);
	renderWindow->AddRenderer(cylinderRenderer);
	renderWindow->AddRenderer(planeRenderer);
	renderWindow->SetSize(800, 600);
	renderWindow->Render();

	//添加鼠标交互  
	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow(renderWindow);
	vtkNew<vtkInteractorStyleTrackballActor> style;
	interactor->SetInteractorStyle(style);

	//初始化交互器 并开始执行事件循环  
	interactor->Initialize();
	interactor->Start();

}