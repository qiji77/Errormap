#pragma once
#include <vtkCellLocator.h>
#include <vtkSmartPointer.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyData.h>
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
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkCellLocator.h>
class renderCone
{
public:
	renderCone()
	{
		renderer =
			vtkSmartPointer<vtkRenderer>::New();
	}
	
	vtkSmartPointer<vtkRenderWindowInteractor> Getrenderwindow() {
		return this->renderWindowInteractor;
	}
	void SetChild(renderCone* child)
	{
		next = child;
	}
	void SetColor(double r, double g, double b)
	{
		renderer->SetBackground(r, g, b); // Background color dark red
	}
	void SetActor(vtkSmartPointer<vtkActor>& getactor)//设置actor
	{
		actor = getactor;
	}
	void SetPos(double* inputpos) {
		double pos[3];
		pos[0] = inputpos[0];
		pos[1] = inputpos[1];
		pos[2] = inputpos[2];
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->Update();

		vtkSmartPointer<vtkPolyDataMapper> mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(sphereSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->SetPosition(pos);
		actor->SetScale(0.2);//设置显示点的大小
		actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		renderer->AddActor(actor);
	}
	void SetInteractor() {
		renderWindowInteractor->SetPicker(pointPicker);
	}
	void Render()
	{
		//Create a cone
		coneSource =
			vtkSmartPointer<vtkConeSource>::New();
		coneSource->Update();

		//Create a mapper and actor
		//mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
		//mapper->SetInputConnection(coneSource->GetOutputPort());

		//actor =
		//	vtkSmartPointer<vtkActor>::New();
		//actor->SetMapper(mapper);


		renderWindow =
			vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);

		//renderWindowInteractor =
		//	vtkSmartPointer<vtkRenderWindowInteractor>::New();

		renderWindowInteractor->SetRenderWindow(renderWindow);

		//Add the actors to the scene
		renderer->AddActor(actor);
		//vtkSmartPointer<PointPickerInteractorStyle> interactorStyle =
		//	vtkSmartPointer<PointPickerInteractorStyle>::New();
		//renderWindowInteractor->SetInteractorStyle(interactorStyle);
		//vtkSmartPointer<vtkInteractorStyleRubberBandPick> interactorStyle = vtkSmartPointer<vtkInteractorStyleRubberBandPick>::New();
		//renderWindowInteractor->SetInteractorStyle(interactorStyle);
		//Render and interact
		renderWindow->Render();
		if (next == nullptr)
		{
			renderWindowInteractor->Start();
		}
	}

public:
	vtkSmartPointer<vtkConeSource>		coneSource;
	vtkSmartPointer<vtkPolyDataMapper>	mapper;
	vtkSmartPointer<vtkActor>			actor;
	vtkSmartPointer<vtkRenderer>		renderer;
	vtkSmartPointer<vtkRenderWindow>	renderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor>	renderWindowInteractor= vtkSmartPointer<vtkRenderWindowInteractor>::New();;
	vtkSmartPointer<vtkPointPicker> pointPicker;
	renderCone* next{ nullptr };
};