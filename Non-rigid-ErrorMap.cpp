// Non-rigid-win32.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <string>
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
#include"rendercone.h"
#include <vtkKdTree.h>
#include <stdio.h>
#include<iostream>
#include<sstream>
#include<vtkLODActor.h>
//#include"algorithm.h"
#include<vtkScalarBarActor.h>
//#include "optimal_nonrigid_icp.h"
#include <vtkAutoInit.h>
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCamera.h"
#include<vtkLineSource.h>
#include <vtkPolyLine.h>
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkRenderingOpenGL2);
#include <vtkLookupTable.h>
typedef std::vector<double> T;
std::vector<T> pickpoints;
std::vector<int> PickId;
std::vector<int> PickId_For_range;
std::vector<int> PickID_Range;
//std::vector<double*> test;
bool ShowTargetCor = false;
vtkSmartPointer<vtkPoints> templet_list= vtkSmartPointer<vtkPoints>::New();
vtkSmartPointer<vtkPoints> target_list = vtkSmartPointer<vtkPoints>::New();
vtkSmartPointer<vtkPoints> source_list = vtkSmartPointer<vtkPoints>::New();
std::vector<double> error_list;
std::vector<double> stiffness_error_list;
std::vector<std::vector<double>> corrent_num;
std::vector<std::vector<double>> X_Matrix_list;
std::vector<std::vector<int>> Neighbor_list;
static renderCone cone;
static renderCone child;
static renderCone source;
void MakeLine(vtkRenderer * renderer);
void Make_range_list() {
	for (int i = 0; i < PickId_For_range.size(); i++) {
		double rxyz[3];
		templet_list->GetPoint(PickId_For_range[i], rxyz);
		for (int j = 0; j < templet_list->GetNumberOfPoints(); j++) {
			double xyz[3];
			templet_list->GetPoint(j, xyz);
			double dis = std::sqrt((rxyz[0] - xyz[0])*(rxyz[0] - xyz[0]) + (rxyz[1] - xyz[1])*(rxyz[1] - xyz[1]) + (rxyz[2] - xyz[2])*(rxyz[2] - xyz[2]));
			if (dis < 0.7) {
				PickID_Range.push_back(j);
			}
		}
	}
}
bool IsIn(std::vector<int> list, int num) {
	for (int i = 0; i < list.size(); i++) {
		if (num == list[i]) {
			return true;
		}
	}
	return false;
}
class PointPickerInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static PointPickerInteractorStyle* New();
	vtkTypeMacro(PointPickerInteractorStyle, vtkInteractorStyleTrackballCamera);

	virtual void OnRightButtonDown()
	{
		std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " " << this->Interactor->GetEventPosition()[1] << std::endl;
		/*vtkSmartPointer < vtkPointPicker> pickertest = vtkSmartPointer < vtkPointPicker>::New();
		pickertest->Pick(this->Interactor->GetEventPosition()[0],
			this->Interactor->GetEventPosition()[1],
			0,  // always zero.
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double pickpoint[3];
		pickertest->GetPickPosition(pickpoint);
		int num = pickertest->GetPointId();
		double xyz[3];
		templet_list->GetPoint(num, xyz);
		std::cout << "Pick id is " << num << std::endl;
		std::cout << "Picked value: " << pickpoint[0] << " " << pickpoint[1] << " " << pickpoint[2] << std::endl;
		std::cout << "Picked value: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
		std::cout << "Distance error is" << error_list[num] << std::endl;
		PickId_For_range.push_back(num);
		if (PickId_For_range.size() == 3) {
			MakeLine(this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		}
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->Update();

		vtkSmartPointer<vtkPolyDataMapper> mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(sphereSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->SetPosition(pickpoint);
		actor->SetScale(0.5);//设置显示点的大小
		actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);
		vtkInteractorStyleTrackballCamera::OnRightButtonDown();
		*/
	
		//this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
		//	this->Interactor->GetEventPosition()[1],
		//	0,  // always zero.
		//	this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		//this->Interactor->GetPicker()->GetPickPosition(picked);
		//int cur_id = 0;
		//cur_id=this->Interactor->GetPicker()->GetPickFromList();
		vtkSmartPointer < vtkPointPicker> pickertest = vtkSmartPointer < vtkPointPicker>::New();
			pickertest->Pick(cone.renderWindowInteractor->GetEventPosition()[0],
				cone.renderWindowInteractor->GetEventPosition()[1],
				0,  // always zero.
				cone.renderWindowInteractor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
			double pickpoint[3];
			pickertest->GetPickPosition(picked);
			int num = pickertest->GetPointId();
			std::cout << "Pick id is " << num << std::endl;
			std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;
			std::cout << "Distance error is" << error_list[num] << std::endl;
			double pos[3];
			pos[0] = corrent_num[num][0];
			pos[1] = corrent_num[num][1];
			pos[2] = corrent_num[num][2];
			//target_list->GetPoint(corrent_num[num],pos);
			vtkSmartPointer<vtkSphereSource> sphereSource =
				vtkSmartPointer<vtkSphereSource>::New();
			sphereSource->Update();

			vtkSmartPointer<vtkPolyDataMapper> mapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputConnection(sphereSource->GetOutputPort());
			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->SetPosition(pos);
			actor->SetScale(0.5);//设置显示点的大小
			actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
			//child.SetActor(actor);
			child.SetPos(pos);

			//源模型上的点
			double source_pos[3];
			source_list->GetPoint(num, source_pos);
			vtkSmartPointer<vtkSphereSource> sphereSource3 =
				vtkSmartPointer<vtkSphereSource>::New();
			sphereSource3->Update();

			vtkSmartPointer<vtkPolyDataMapper> mapper3 =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper3->SetInputConnection(sphereSource3->GetOutputPort());
			vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
			actor3->SetMapper(mapper);
			actor3->SetPosition(pos);
			actor3->SetScale(0.5);//设置显示点的大小
			actor3->GetProperty()->SetColor(1.0, 0.0, 0.0);
			//child.SetActor(actor);
			source.SetPos(source_pos);
			//源模型的点（上部）
			//child->AddActor(actor);
			vtkSmartPointer<vtkSphereSource> sphereSource2 =
				vtkSmartPointer<vtkSphereSource>::New();
			sphereSource2->Update();

			vtkSmartPointer<vtkPolyDataMapper> mapper2 =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper2->SetInputConnection(sphereSource2->GetOutputPort());
			vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
			actor2->SetMapper(mapper2);
			actor2->SetPosition(picked);
			actor2->SetScale(0.1);//设置显示点的大小
			actor2->GetProperty()->SetColor(1.0, 0.0, 0.0);

			cone.renderWindowInteractor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);

			vtkInteractorStyleTrackballCamera::OnRightButtonDown();
		
	}
};

vtkStandardNewMacro(PointPickerInteractorStyle);
void InitLineActor(vtkPoints* points, vtkActor* lactor, double color[3]) {
	vtkSmartPointer<vtkPolyLine> polyLine =vtkSmartPointer<vtkPolyLine>::New();
	
	polyLine->GetPointIds()->SetNumberOfIds(points->GetNumberOfPoints());
	for (unsigned int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		polyLine->GetPointIds()->SetId(i, i);
	}

	vtkSmartPointer<vtkCellArray> cells =
		vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(polyLine);
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(points);
	polyData->SetLines(cells);
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(polyData);

	lactor->SetMapper(mapper);
	lactor->GetProperty()->SetLineWidth(1);
	lactor->GetProperty()->SetColor(color);
}
void MakeLine(vtkRenderer * renderer) {
	Make_range_list();
	for (int i = 0; i < PickID_Range.size(); i++) {
		vtkPoints* points = vtkPoints::New();
		double xyz[3];
		int num = PickID_Range[i];
		templet_list->GetPoint(num, xyz);
		points->InsertNextPoint(xyz[0], xyz[1], xyz[2]);
		points->InsertNextPoint(corrent_num[num][0] + 30, corrent_num[num][1], corrent_num[num][2]);
		vtkSmartPointer<vtkActor> temp_actor = vtkSmartPointer<vtkActor>::New();
		double color[3] = { 1,0,0 };
		InitLineActor(points, temp_actor, color);
	    renderer->AddActor(temp_actor);
	}
	for (int i = 0; i < PickID_Range.size(); i++) {
		vtkPoints* points = vtkPoints::New();
		double xyz[3];
		int num = PickID_Range[i];
		templet_list->GetPoint(num, xyz);
		points->InsertNextPoint(xyz[0], xyz[1], xyz[2]);
		double sxyz[3];
		source_list->GetPoint(num, sxyz);
		points->InsertNextPoint(sxyz[0] - 30, sxyz[1], sxyz[2]);
		vtkSmartPointer<vtkActor> temp_actor = vtkSmartPointer<vtkActor>::New();
		double color[3] = { 0,1,0 };
		InitLineActor(points, temp_actor, color);
		renderer->AddActor(temp_actor);
	}
}
int main(int argc, char**argv)
{
	    std::string source_filename = "BaseModel.ply";
		std::string template_filename = "BaseModel_choose_mouth_rigid_no_cal_199.ply";
		std::string target_filename = "ScreamNew.ply";

		vtkSmartPointer<vtkPLYReader> template_reader = vtkSmartPointer<vtkPLYReader>::New();
		template_reader->SetFileName(template_filename.c_str());
		template_reader->Update();

		vtkSmartPointer<vtkPLYReader> target_reader = vtkSmartPointer<vtkPLYReader>::New();
		target_reader->SetFileName(target_filename.c_str());
		target_reader->Update();

		vtkSmartPointer<vtkPolyData> template_polyData = template_reader->GetOutput();
		vtkSmartPointer<vtkPolyData> target_polyData = target_reader->GetOutput();
		templet_list = template_polyData->GetPoints();
		target_list = target_polyData->GetPoints();
		vtkSmartPointer<vtkPLYReader> source_reader = vtkSmartPointer<vtkPLYReader>::New();
		source_reader->SetFileName(source_filename.c_str());
		source_reader->Update();
		vtkSmartPointer<vtkPolyData> source_Polydata = source_reader->GetOutput();
		vtkSmartPointer<vtkPolyDataMapper> source_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		source_mapper->SetInputData(source_Polydata);
		vtkSmartPointer<vtkActor> source_actor = vtkSmartPointer<vtkActor>::New();
		source_actor->SetMapper(source_mapper);
		source_actor->GetProperty()->SetColor(0, 0.8, 0.2);
		source_actor->GetProperty()->SetOpacity(0.7);
		source_actor->SetPosition(-30, 0, 0);
		vtkSmartPointer<vtkInteractorStyleRubberBandPick> interactorStyle3 =
			vtkSmartPointer<vtkInteractorStyleRubberBandPick>::New();
		source_list = source_Polydata->GetPoints();
		//vtkSmartPointer<vtkActor> template_actor = vtkSmartPointer<vtkActor>::New();
		//template_actor->SetMapper(template_mapper);
	//	InitLineActor(points, template_actor, color);
		ifstream cor_text;
		cor_text.open("199BaseModelcor_no_cal.txt");
		std::string str_cor;
		while (std::getline(cor_text, str_cor)) {
			std::istringstream input(str_cor);
			std::vector<double> a;
			double temp;
			while (input >> temp) {
				a.push_back(temp);
			}
			corrent_num.push_back(a);
		}
		//GetVertices(templet_list, corrent_list);
		double max_error=0;
		double max_sec_error = 0;
		double min_error=100;
		
		for (int i = 0; i < templet_list->GetNumberOfPoints(); i++) {
			double xyz[3];
			templet_list->GetPoint(i, xyz);
			double txyz[3];
			txyz[0]=corrent_num[i][0];
			txyz[1] = corrent_num[i][1];
			txyz[2] = corrent_num[i][2];
			double error = std::sqrt((xyz[0] - txyz[0])*(xyz[0] - txyz[0]) + (xyz[1] - txyz[1])*(xyz[1] - txyz[1])
				+ (xyz[2] - txyz[2])*(xyz[2] - txyz[2]));
			error = error;
			error_list.push_back(error);
			if (max_error < error) {
				max_sec_error = max_error;
				max_error = error;

			}
			if (min_error > error) {
				min_error = error;
			}
		}
		//计算刚性loss
		ifstream mat_text;
		mat_text.open("199BaseModelX_Matrix_no_cal.txt");
		std::string str;
		while (std::getline(mat_text, str)) {
			std::istringstream input(str);
			std::vector<double> a;
			double temp;
			while (input >> temp) {
				a.push_back(temp);
			}
			X_Matrix_list.push_back(a);
		}
		//计算邻居
		int num = template_polyData->GetNumberOfCells();
		std::vector<std::vector<int>> Topl_list(template_polyData->GetNumberOfPoints());
		for (int i = 0; i < template_polyData->GetNumberOfCells(); ++i)//遍历模型中的单元，找到所有边，存储每个点的和该点相连的点的ID
		{
			vtkCell* cell = template_polyData->GetCell(i);
			int a = cell->GetPointId(0);
			int b = cell->GetPointId(1);
			int c = cell->GetPointId(2);
			if (IsIn(Topl_list[a], b) == false) {//判断该点是否已经在周围点列表中
				Topl_list[a].push_back(b);
			}
			if (IsIn(Topl_list[a], c) == false) {
				Topl_list[a].push_back(c);
			}
			if (IsIn(Topl_list[b], a) == false) {
				Topl_list[b].push_back(a);
			}
			if (IsIn(Topl_list[b], c) == false) {
				Topl_list[b].push_back(c);
			}
			if (IsIn(Topl_list[c], b) == false) {
				Topl_list[c].push_back(b);
			}
			if (IsIn(Topl_list[c],a) == false) {
				Topl_list[c].push_back(a);
			}
		

		}
		//设置最大误差和最小误差，在映射误差图颜色时用
		double max_stiff_error = 0;
		double min_stiff_error = 100;
		for (int i = 0; i < template_polyData->GetNumberOfPoints(); i++) {//遍历原始模型的点
			double testPoint[3];
			template_polyData->GetPoint(i, testPoint);//得到原始模型上的第i个点
			vtkSmartPointer<vtkIdList> result =
				vtkSmartPointer<vtkIdList>::New();
			//pointTree->FindClosestNPoints(k, testPoint, result);//找到4个最近点，结果保存在result中
			double cur_offset[3];//当前点的变换的位置（这里的X_Matrix里保存的XT，是论文中用来变换坐标的矩阵，这里将点扩充为齐次坐标表示）
			cur_offset[0] = X_Matrix_list[3 * i][0] * testPoint[0] + X_Matrix_list[3 * i][1] * testPoint[1] +
				X_Matrix_list[3 * i][2] * testPoint[2] + X_Matrix_list[3 * i][3];
			cur_offset[1] = X_Matrix_list[3 * i+1][0] * testPoint[0] + X_Matrix_list[3 * i+1][1] * testPoint[1] +
				X_Matrix_list[3 * i+1][2] * testPoint[2] + X_Matrix_list[3 * i+1][3];
			cur_offset[2] = X_Matrix_list[3 * i + 2][0] * testPoint[0] + X_Matrix_list[3 * i + 2][1] * testPoint[1] +
				X_Matrix_list[3 * i + 2][2] * testPoint[2] + X_Matrix_list[3 * i + 2][3];
			double near_list[10][3];
			double mat_err_temp=0;//定义error
			for (int j = 0; j < Topl_list[i].size(); j++) {
				int num = Topl_list[i][j];//得到最近点的id
				near_list[j][0]= X_Matrix_list[3 * num][0] * testPoint[0] + X_Matrix_list[3 * num][1] * testPoint[1] +
					X_Matrix_list[3 * num][2] * testPoint[2] + X_Matrix_list[3 * num][3];
				near_list[j][1] = X_Matrix_list[3 * num+1][0] * testPoint[0] +
					X_Matrix_list[3 * num+1][1] * testPoint[1] + X_Matrix_list[3 * num+1][2] * testPoint[2] + X_Matrix_list[3 * num+1][3];
				near_list[j][2] = X_Matrix_list[3 * num+2][0] * testPoint[0] +
					X_Matrix_list[3 * num+2][1] * testPoint[1] + X_Matrix_list[3 * num+2][2] * testPoint[2] + X_Matrix_list[3 * num+2][3];
			}
			for (int er_num = 0; er_num < Topl_list[i].size(); er_num++) {//计算loss
				mat_err_temp += std::sqrt((near_list[er_num][0] - cur_offset[0])*(near_list[er_num][0] - cur_offset[0]) +
					(near_list[er_num][1] - cur_offset[1])*(near_list[er_num][1] - cur_offset[1]) +
					(near_list[er_num][2] - cur_offset[2])*(near_list[er_num][2] - cur_offset[2]));
			}
			mat_err_temp = mat_err_temp / Topl_list[i].size();
			stiffness_error_list.push_back(mat_err_temp);
			if (max_stiff_error < mat_err_temp) {
				max_stiff_error = mat_err_temp;
			}
			if (min_stiff_error > mat_err_temp) {
				min_stiff_error = mat_err_temp;
			}
		}
		
		
		
		vtkSmartPointer<vtkFloatArray>scalars = vtkSmartPointer<vtkFloatArray>::New();
		for (int i = 0; i < template_polyData->GetNumberOfPoints(); i++){
			//double temp = error_list[i];
			double temp = error_list[i];
			scalars->InsertTuple1(i, temp);
		}

		template_polyData->GetPointData()->SetScalars(scalars);
		vtkSmartPointer<vtkPolyDataMapper> template_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		template_mapper->SetInputData(template_polyData);
		vtkSmartPointer<vtkLookupTable>pColorTable = vtkSmartPointer<vtkLookupTable>::New();
		int numofpoints = template_polyData->GetNumberOfPoints();
		pColorTable->SetHueRange(0.6, 0.0);
		pColorTable->SetAlphaRange(1.0, 1.0);
		pColorTable->SetValueRange(1.0, 1.0);
		pColorTable->SetSaturationRange(1.0, 1.0);
		pColorTable->SetNumberOfTableValues(25600);
		pColorTable->Build();
		template_mapper->SetScalarRange(min_stiff_error,max_stiff_error);
		template_mapper->SetLookupTable(pColorTable);
		vtkSmartPointer<vtkScalarBarActor> scalarBar =
			vtkSmartPointer<vtkScalarBarActor>::New();
		scalarBar->SetLookupTable(template_mapper->GetLookupTable());
		scalarBar->SetTitle("Error Bar");
		scalarBar->SetNumberOfLabels(5); //设置5个标签  
		vtkSmartPointer<vtkPolyDataMapper> target_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		target_mapper->SetInputData(target_polyData);
		vtkSmartPointer<vtkActor> template_actor = vtkSmartPointer<vtkActor>::New();
		template_actor->SetMapper(template_mapper);
		template_actor->SetPosition(0, 0, 0);
		//template_actor->GetProperty()->SetOpacity(0.7);
		//template_actor->GetProperty()->SetColor(0.8, 0.2, 0);


		vtkSmartPointer<vtkActor> target_actor = vtkSmartPointer<vtkActor>::New();
		target_actor->SetMapper(target_mapper);
		target_actor->SetPosition(30, 0, 0);
		target_actor->GetProperty()->SetColor(0, 0.2, 0.8);
		target_actor->GetProperty()->SetOpacity(0.7);

		vtkCamera *camera = vtkCamera::New();
		camera->SetPosition(1, 1, 1);
		camera->SetFocalPoint(0, 0, 0);

		vtkRenderer *renderer = vtkRenderer::New();
		vtkRenderWindow *renWin = vtkRenderWindow::New();
		renWin->AddRenderer(renderer);
		vtkSmartPointer<PointPickerInteractorStyle> interactorStyle_figer =
			vtkSmartPointer<PointPickerInteractorStyle>::New();
		//vtkSmartPointer<vtkInteractorStyleRubberBandPick> interactorStyle4 =
		//	vtkSmartPointer<vtkInteractorStyleRubberBandPick>::New();
		vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
		iren->SetRenderWindow(renWin);
		
		//renderer->AddActor(source_actor);
		renderer->AddActor(template_actor);
		//renderer->AddActor(target_actor);
		renderer->AddActor2D(scalarBar);
		renderer->SetActiveCamera(camera);
		renderer->ResetCamera();
		renderer->SetBackground(1, 1, 1);
		iren->SetInteractorStyle(interactorStyle_figer);
		renWin->SetSize(400, 400);
		renWin->Render();
		iren->Start();
		
	    vtkSmartPointer<PointPickerInteractorStyle> interactorStyle =
			vtkSmartPointer<PointPickerInteractorStyle>::New();
		cone.Getrenderwindow()->SetInteractorStyle(interactorStyle);

		//cone.SetColor(0.5, 0, 0);
		cone.SetActor(template_actor);
		cone.SetInteractor();

		vtkSmartPointer<PointPickerInteractorStyle> interactorStyle2 =
			vtkSmartPointer<PointPickerInteractorStyle>::New();
		child.Getrenderwindow()->SetInteractorStyle(interactorStyle2);
		child.SetColor(0, 0.5, 0);
		child.SetActor(target_actor);
		cone.SetChild(&child);
	
		source.Getrenderwindow()->SetInteractorStyle(interactorStyle3);
		source.SetActor(source_actor);
		source.SetColor(0, 0, 5);
		child.SetChild(&source);

		cone.Render();
		child.Render();
		source.Render();
		
	return 0;
}
// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
