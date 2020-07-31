#pragma once
#ifndef OPTIMAL_NONRIGID_ICP_H
#define OPTIMAL_NONRIGID_ICP_H

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkCellLocator.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <vector>

#include <boost/shared_ptr.hpp>


typedef vtkSmartPointer<vtkPolyData> Mesh;

typedef Mesh Template;
typedef Mesh Target;
typedef std::pair<int, int> Edge;
typedef boost::shared_ptr< std::vector<Edge> > Edges;
typedef vtkSmartPointer<vtkPoints> Vertices;
typedef boost::shared_ptr< std::vector<float> > Weights;

class OptimalNonrigidICP
{
public:
	OptimalNonrigidICP(Template _template, Target _target, std::vector<std::vector<double>>& pickpoints,std::vector<int>& picknum) :_template(_template), _target(_target), pickpoints(pickpoints),PickNum(picknum){}//构造函数
	~OptimalNonrigidICP() {}

	void init()//首先调用该函数
	{
		edgesInit();//初始化边
		verticesInit();//初始化顶点
		nearestSearchInit();//初始化最邻近搜索
	}

	void initCompute()//第二步调用，初始化计算
	{
		//ComputeDis();//计算偏移距离
		correspondencesInit();//初始化计算相关点
		weightsInit();//权重计算初始化
		InitRangNum();//计算

	}

	void edgesInit()
	{


		if (_edges == NULL) _edges.reset(new std::vector<Edge>);
		int num = _template->GetNumberOfCells();
		for (int i = 0; i < _template->GetNumberOfCells(); ++i)//GetNumberOfCells()函数用于获取模型的单元数，一个单元有三个边
		{
			vtkCell* cell = _template->GetCell(i);
			int a = cell->GetPointId(0);
			int b = cell->GetPointId(1);
			int c = cell->GetPointId(2);

			_edges->push_back(Edge(a, b));
			_edges->push_back(Edge(b, c));
			_edges->push_back(Edge(c, a));
			
		}
		for (int i = 0; i < _edges->size(); i++) {
			Edge edge = (*_edges)[i];
			int a = edge.first;
			int b = edge.second;
			if (GetInRangeNum(a) || GetInRangeNum(b)) {
				addaloha.push_back(true);
			}
			else {
				addaloha.push_back(false);
			}
		}
		
	}

	void nearestSearchInit()
	{
		_cellLocator = vtkSmartPointer<vtkCellLocator>::New();
		_cellLocator->SetDataSet(_target);
		_cellLocator->BuildLocator();//构建定位器		
	}

	void verticesInit()
	{
		_vertices = _template->GetPoints();//顶点数
		target_vertices = _target->GetPoints();
	}
	int GetIdinList(int num) {
		for (int i = 0; i < PickNum.size()-20; i++) {
			if (num == PickNum[i]) {
				return i;
			}
		}
		return -1;
	}
	int computeintemplelist(double* xyz, std::vector<std::vector<double>> pickpoints, int num_start, int num_end,double dis=2) {//计算是否在列表中，如果在则返回对应的队列位置，否则返回-1；
		for (int num_pick = num_start; num_pick < num_end; num_pick++) {
			double dis = std::sqrt((xyz[0] - pickpoints[num_pick][0])*(xyz[0] - pickpoints[num_pick][0]) +
				(xyz[1] - pickpoints[num_pick][1])*(xyz[1] - pickpoints[num_pick][1]) +
				(xyz[2] - pickpoints[num_pick][2])*(xyz[2] - pickpoints[num_pick][2]));
			if (dis < 0.05) {
				
				return num_pick;
			}
			if (dis<2) {
				PickRangeNum.push_back(num_pick);
			}
		}
		return -1;
	}
	void InitRangNum() {
		for (int i = 0; i < _vertices->GetNumberOfPoints(); i++) {
			double xyz[3];
			_vertices->GetPoint(i, xyz);
			for (int num_pick = 0; num_pick < 20; num_pick++) {
				double dis = std::sqrt((xyz[0] - pickpoints[num_pick][0])*(xyz[0] - pickpoints[num_pick][0]) +
					(xyz[1] - pickpoints[num_pick][1])*(xyz[1] - pickpoints[num_pick][1]) +
					(xyz[2] - pickpoints[num_pick][2])*(xyz[2] - pickpoints[num_pick][2]));
				if (dis < 2) {
					PickRangeNum.push_back(num_pick);
				}
			}
		}
	}
	bool GetInRangeNum(int num) {
		for (int i = 0; i < PickRangeNum.size(); i++) {
			if (num == PickRangeNum[i]) {
				return true;
			}
		}
		return false;
	}
	int CosmputeRangelist(double* xyz, std::vector<std::vector<double>> pickpoints, int num_start, int num_end, double dis = 2) {//计算是否在列表中，如果在则返回对应的队列位置，否则返回-1；
		for (int num_pick = num_start; num_pick < num_end; num_pick++) {
			double dis = std::sqrt((xyz[0] - pickpoints[num_pick][0])*(xyz[0] - pickpoints[num_pick][0]) +
				(xyz[1] - pickpoints[num_pick][1])*(xyz[1] - pickpoints[num_pick][1]) +
				(xyz[2] - pickpoints[num_pick][2])*(xyz[2] - pickpoints[num_pick][2]));
			if (dis < dis) {
				return num_pick;
			}
		}
		return -1;
	}
	void correspondencesInit()//设置最初的相关点
	{
		if (_correspondences == NULL) _correspondences = Vertices::New();

		_correspondences->SetNumberOfPoints(_vertices->GetNumberOfPoints());
		for (int i = 0; i < _vertices->GetNumberOfPoints(); i++)
		{
			double testPoint[3];
			_vertices->GetPoint(i, testPoint);
			double closestPoint[3];
			double closestPointDist2;
			vtkIdType cellId;
			int subId;
			/*if (computeintemplelist(testPoint, pickpoints, 0, 20) > -1) {
				double targetpos[3];
				targetpos[0] = pickpoints[computeintemplelist(testPoint, pickpoints, 0, 20) + 20][0];
				targetpos[1] = pickpoints[computeintemplelist(testPoint, pickpoints, 0, 20) + 20][1];
				targetpos[2] = pickpoints[computeintemplelist(testPoint, pickpoints, 0, 20) + 20][2];
				_correspondences->SetPoint(i, targetpos);
			}*/
			int num = GetIdinList(i);
			if (num>-1) {
				double targetpos[3];
				int cor_id = num+ 20;
				targetpos[0] = pickpoints[cor_id][0];
				targetpos[1] = pickpoints[cor_id][1];
				targetpos[2] = pickpoints[cor_id][2];
				_correspondences->SetPoint(i, targetpos);
			}
			else {
				_cellLocator->FindClosestPoint(testPoint, closestPoint, cellId, subId, closestPointDist2);
				_correspondences->SetPoint(i, closestPoint);
			}
		}
	}
	void GetVertices(vtkSmartPointer<vtkPoints>&  templet_list, vtkSmartPointer<vtkPoints>& corrent_list) {
		templet_list->SetNumberOfPoints(_vertices->GetNumberOfPoints()); 
		corrent_list->SetNumberOfPoints(_correspondences->GetNumberOfPoints());
		for (int i = 0; i < _vertices->GetNumberOfPoints(); i++) {
			double xyz[3];
			_vertices->GetPoint(i, xyz);
			templet_list->SetPoint(i, xyz);
		}
		for (int i = 0; i < _correspondences->GetNumberOfPoints(); i++) {
			double xyz[3];
			_correspondences->GetPoint(i, xyz);
			corrent_list->SetPoint(i, xyz);
		}
	}
	
	void weightsInit()//初始化权重，都赋值1
	{
		if (_weights == NULL) _weights.reset(new std::vector<float>());
		_weights->resize(_vertices->GetNumberOfPoints());
		for (int i = 0; i < _vertices->GetNumberOfPoints(); i++) {

			(*_weights)[i] = 1.0f;
			double xyz[3];
			_vertices->GetPoint(i, xyz);
			if (Addweight(xyz, pickpoints, 0, 1)) {
				(*_weights)[i] = 2.0f;
			}
		}
	}
	bool Addweight(double* xyz, std::vector<std::vector<double>> pickpoints, int num_start, int num_end) {
		for (int num_pick = num_start; num_pick < num_end; num_pick++) {
			double dis = std::sqrt((xyz[0] - pickpoints[num_pick][0])*(xyz[0] - pickpoints[num_pick][0]) +
				(xyz[1] - pickpoints[num_pick][1])*(xyz[1] - pickpoints[num_pick][1]) +
				(xyz[2] - pickpoints[num_pick][2])*(xyz[2] - pickpoints[num_pick][2]));
			if (dis < 0.05) {
				return true;
			}
		}
		return false;
	}
	int compute(float alpha, float beta, float gamma)
	{

		stiffness = 1;
		lastmatrix.setZero();
		//To do nonrigid icp registration
		//while (stiffness > minstiffness) {

		int loop = 1;
		int firstnum = 1;
		while (loop) {
			int n = _vertices->GetNumberOfPoints();
			int m = _edges->size();

			Eigen::SparseMatrix<float> A(4 * m + n, 4 * n);//A矩阵

			std::vector< Eigen::Triplet<float> > alpha_M_G;
			for (int i = 0; i < m; ++i)
			{
				Edge edge = (*_edges)[i];
				int a = edge.first;
				int b = edge.second;
				/*if (addaloha[i]) {
					for (int j = 0; j < 3; j++)
						alpha_M_G.push_back(Eigen::Triplet<float>(i * 4 + j, a * 4 + j, alpha*2));
					alpha_M_G.push_back(Eigen::Triplet<float>(i * 4 + 3, a * 4 + 3, alpha * 2 * gamma));//构建（Es）M_G矩阵

					for (int j = 0; j < 3; j++)
						alpha_M_G.push_back(Eigen::Triplet<float>(i * 4 + j, b * 4 + j, -alpha * 2));
					alpha_M_G.push_back(Eigen::Triplet<float>(i * 4 + 3, b * 4 + 3, -alpha * 2 * gamma));
				}*/
			
					for (int j = 0; j < 3; j++)
						alpha_M_G.push_back(Eigen::Triplet<float>(i * 4 + j, a * 4 + j, alpha));
					alpha_M_G.push_back(Eigen::Triplet<float>(i * 4 + 3, a * 4 + 3, alpha * gamma));//构建（Es）M_G矩阵

					for (int j = 0; j < 3; j++)
						alpha_M_G.push_back(Eigen::Triplet<float>(i * 4 + j, b * 4 + j, -alpha));
					alpha_M_G.push_back(Eigen::Triplet<float>(i * 4 + 3, b * 4 + 3, -alpha * gamma));
				
				
			}
			std::cout << "alpha_M_G calculated!" << std::endl;

			std::vector< Eigen::Triplet<float> > W_D;
			for (int i = 0; i < n; ++i)
			{
				double xyz[3];
				_vertices->GetPoint(i, xyz);

				float weight = (*_weights)[i];
				for (int j = 0; j < 3; ++j) W_D.push_back(Eigen::Triplet<float>(4 * m + i, i * 4 + j, weight * xyz[j]));//Triplet是三元组，这里应为点的坐标
				W_D.push_back(Eigen::Triplet<float>(4 * m + i, i * 4 + 3, weight));
			}
			std::cout << "W_D calculated!" << std::endl;

			std::vector< Eigen::Triplet<float> > _A = alpha_M_G;
			_A.insert(_A.end(), W_D.begin(), W_D.end());
			std::cout << "_A calculated!" << std::endl;

			A.setFromTriplets(_A.begin(), _A.end());
			std::cout << "A calculated!" << std::endl;

			Eigen::MatrixX3f B = Eigen::MatrixX3f::Zero(4 * m + n, 3);
			for (int i = 0; i < n; ++i)
			{
				double xyz[3];
				_correspondences->GetPoint(i, xyz);
				float weight = (*_weights)[i];
				for (int j = 0; j < 3; j++)
					B(4 * m + i, j) = weight * xyz[j];
			}
			std::cout << "B calculated!" << std::endl;
			
			Eigen::SparseMatrix<float> ATA = Eigen::SparseMatrix<float>(A.transpose()) * A;
			std::cout << "ATA calculated!" << std::endl;
			Eigen::MatrixX3f ATB = Eigen::SparseMatrix<float>(A.transpose()) * B;
			std::cout << "ATB calculated!" << std::endl;

			Eigen::ConjugateGradient< Eigen::SparseMatrix<float> > solver;
			solver.compute(ATA);
			std::cout << "solver computed ATA!" << std::endl;
			if (solver.info() != Eigen::Success)
			{
				std::cerr << "Decomposition failed" << std::endl;
				return 1;
			}
			Eigen::MatrixX3f X = solver.solve(ATB);
			std::cout << "X calculated!" << std::endl;

			Eigen::Matrix3Xf XT = X.transpose();
			for (int i = 0; i < n; ++i)
			{

				double xyz[3];
				_vertices->GetPoint(i, xyz);
				Eigen::Vector4f point(xyz[0], xyz[1], xyz[2], 1.0f);
				Eigen::Vector3f point_transformed = XT.block<3, 4>(0, 4 * i) * point;//从(0,4i)开始的3*4变形矩阵
				_vertices->SetPoint(i, point_transformed[0], point_transformed[1], point_transformed[2]);//根据计算结果变换点

				double txyz[3];
				_correspondences->GetPoint(i, txyz);//得到相关点
				//if ( i < 10) std::cout << XT.block<3, 4>(0, 4*i) << std::endl;
				if (i < 10) std::cout << xyz[0] << "," << xyz[1] << "," << xyz[2] << " -> "
					<< point_transformed[0] << " " << point_transformed[1] << " " << point_transformed[2] << " -> "
					<< txyz[0] << "," << txyz[1] << "," << txyz[2] << std::endl;

			}
			double norm = epsilon + 1;
			if (firstnum > 1) {
				norm = (X - lastmatrix).norm();
			}
			//std::cout << "The norm value in this round (inner loop): " << norm << std::endl << std::endl;
			if (norm < epsilon || norm != norm) {
				loop = 0;
				if (norm != norm) {
					cout << "Bad Matrix Error!" << endl;
					//errorInt = 1;
				}
			}
			firstnum += 1;
			lastmatrix = X;
		}

		stiffness -= 0.1;
		//}

		return 0;
	}

protected:
	Template _template;
	Target _target;

	Edges _edges;
	Vertices _vertices;
	Vertices target_vertices;
	Vertices _correspondences;
	Weights _weights;
	float epsilon = 2;
	float stiffness = 1;
	float minstiffness = 0.1;
	double distance[3] = { 0, 0, 0 };
	Eigen::MatrixX3f lastmatrix;
	std::vector<std::vector<double>> pickpoints;//源模型中选中点的列表
	std::vector<int> PickRangeNum;//选中点的范围内的点
	std::vector<int> PickNum;//选中的点的ID
	std::vector<bool> addaloha;//是否增加alpha
private:
	vtkSmartPointer<vtkCellLocator> _cellLocator;
};

#endif


