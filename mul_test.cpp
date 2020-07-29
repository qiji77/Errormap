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
	// ����Բ׶����
	vtkNew<vtkConeSource> cone;
	cone->SetHeight(20);         //����Բ׶�ĸ߶ȺͰ뾶  
	cone->SetRadius(20);
	cone->SetResolution(100);    //���÷ֱ��ʣ�ֵԽ��Խ������Բ׶��
	//������Ӧ��mapper 
	vtkNew<vtkPolyDataMapper> coneMapper;
	coneMapper->SetInputConnection(cone->GetOutputPort());
	// ������Ӧ��Actor
	vtkNew<vtkActor> coneActor;
	coneActor->SetMapper(coneMapper);
	// ���ÿ��actor����һ��vtkRender�����actor 
	vtkNew<vtkRenderer> coneRenderer;
	coneRenderer->AddActor(coneActor);
	//����Ĵ�����������ӿڵģ�����Ļ��Ϊ4���ӿ�  
	//ǰ�����������ӿ����½ǵ������(xmin,ymin)�����������������Ͻǵ�����(xmax,ymax)  
	coneRenderer->SetViewport(0, 0.5, 0.5, 1);
	coneRenderer->SetBackground(0.1, 0.2, 0.4);

	// ����������
	vtkNew<vtkCubeSource> cube;
	cube->SetXLength(40);     //����������ĳ����  
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

	//����Բ����  
	vtkNew<vtkCylinderSource> cylinder;
	vtkNew<vtkPolyDataMapper> cylinderMapper;
	cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
	vtkNew<vtkActor> cylinderActor;
	cylinderActor->SetMapper(cylinderMapper);
	vtkNew<vtkRenderer> cylinderRenderer;
	cylinderRenderer->AddActor(cylinderActor);
	cylinderRenderer->SetViewport(0, 0, 0.5, 0.5);

	//����ƽ��  
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

	//�����꽻��  
	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow(renderWindow);
	vtkNew<vtkInteractorStyleTrackballActor> style;
	interactor->SetInteractorStyle(style);

	//��ʼ�������� ����ʼִ���¼�ѭ��  
	interactor->Initialize();
	interactor->Start();

}