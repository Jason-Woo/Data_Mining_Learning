/*
datasets\\Aggregation_cluster=7.txt
datasets\\flame_cluster=2.txt
datasets\\Jain_cluster=2.txt
datasets\\Pathbased_cluster=3.txt
datasets\\Spiral_cluster=3.txt
*/
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>
#include<math.h>
#include<Eigen/Dense> 
#include<stdlib.h>
#include<time.h>

using namespace std;
using namespace Eigen; 

struct point
{
    double x;
    double y;
    int flag;
};

int k;
vector<vector<double> > y;

struct cluster
{
    vector<double> center;
    vector<int> id;
    void update_center()
    {
        double sum_x = 0;
        double sum_y = 0;
        for(int i = 0 ; i < k ; i ++)
        {
        	double tmp_sum = 0 ;
        	for(int j = 0 ;j < id.size(); j ++)
        	{
        		tmp_sum += y[i][id[j]];
			}
			tmp_sum /= id.size();
			center[k] = tmp_sum;
		}
    };
};

vector<double> eigen_val,D;
vector<vector<double> > W, L, L0, eigen_vec;
double sigma;
int num_cluster;
vector<point> P;
string input_file, output_file;
vector<cluster> C0;
int iteration;

double cal_dis(point p1, point p2)
{
    double tmp = pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2);
    double tmp2 = tmp / (-1 * 2 * sigma * sigma);
    return exp(tmp2);
}

void Init()
{
    k = 15;
    iteration = 2000;
    num_cluster = 7;
	sigma = 5;
    input_file = "datasets\\Aggregation_cluster=7.txt";
    output_file = "result\\spectral_Aggregation_cluster=7.txt";
    C0.resize(num_cluster);
    for(int i = 0; i < num_cluster; i ++) C0[i].center.resize(k);
}

void Input(string file_name)
{
    ifstream fin(file_name);
    if(! fin)
    {
        cout << "OPEN FILE ERROR" << endl;
        return;
    }

    string buffer;
    int input_cnt = 0;
    while(getline(fin, buffer))
    {
        point p;
		stringstream ss(buffer);
        string tmp;
        ss >> tmp;
        double tmp_num = stod(tmp);
        p.x = tmp_num;
        ss >> tmp;
        tmp_num = stod(tmp);
        p.y = tmp_num;
        p.flag = -1;
        P.push_back(p);
        input_cnt ++;
    }
}

void Generate()
{
    int size = P.size();
    D.resize(size);
    W.resize(size);
    for(int i = 0; i < size; i ++) W[i].resize(size);
    L.resize(size);
    for(int i = 0; i < size; i ++) L[i].resize(size);
    L0.resize(size);
    for(int i = 0; i < size; i ++) L0[i].resize(size);
    for(int i = 0; i < size; i ++)
    {
        W[i][i] = 0;
        for(int j = i + 1; j < size; j ++)
        {
            W[i][j] = cal_dis(P[i], P[j]);
            W[j][i] = W[i][j];
        }
    }
    for(int i = 0; i < size; i ++)
    {
    	double tmp_w = 0;
    	for(int j = 0; j < size; j ++)
    	{
    		tmp_w += W[i][j];
		}
		D[i] = tmp_w;
	}
    for(int i = 0; i < size; i ++)
    {
        L[i][i] = D[i] - W[i][i];
        for(int j = i + 1; j < size; j ++)
        {
            L[i][j] = 0 - W[i][j];
            L[j][i] = L[i][j];
        }
    }
    for(int i = 0; i < size; i ++)
    {
        for(int j = 0; j < size; j ++) L0[i][j] = L[i][j] / sqrt(D[i] * D[j]);
    }
    cout << "Running Matrix!" << endl;
    Matrix<double, Dynamic, Dynamic> U;
    U.resize(size, size);
    for(int i = 0; i < size; i ++)
    {
        for(int j = 0; j < size; j ++) U(i, j) = L0[i][j];
    }
    cout << "Matrix Assignment Success!" << endl;
    EigenSolver<MatrixXd> es(U);
    eigen_val.resize(size);
	eigen_vec.resize(size);
	for(int i = 0; i < size; i ++) eigen_vec[i].resize(size);
	y.resize(k);
	for(int i = 0; i < k; i ++) y[i].resize(size);
	auto au = es.eigenvalues();
	for(int i = 0; i < size; i ++)
	{
		eigen_val[i] = au[i].real();
	} 
	cout << "Matrix Eigen Values Success!" << endl;
	Matrix<complex<double>, Dynamic, Dynamic> m = es.eigenvectors(); 
	for(int i = 0; i < size; i ++)
	{
		for(int j = 0; j < size ; j ++)
		{
			eigen_vec[i][j] = m.row(i)[j].real();
		} 
	}
	cout << "Matrix Eigen Vector Success!" << endl;
	vector<double> eigen_val_0 = eigen_val;
	sort(eigen_val.begin(), eigen_val.end());
	for(int i = 0; i < k ; i ++)
	{
		vector<double>::iterator it = find(eigen_val_0.begin(),eigen_val_0.end(),eigen_val[i]);
		if(it == eigen_val_0.end())
		{
			cout << "ERROR0!" << endl;
			return;
		}

		for(int j = 0; j < size; j ++)
		{
			y[i][j] = eigen_vec[it - eigen_val_0.begin()][j];
		}		
	}
	for(int i = 0; i < size ; i ++)
	{
		double tmp_sum = 0;
		for(int j = 0 ; j < k; j ++)
		{
			tmp_sum += pow(y[j][i], 2);
		}
		tmp_sum = sqrt(tmp_sum);
		for(int j = 0 ; j < k; j ++)
		{
			y[j][i] /= tmp_sum;
		}
	} 
}

double cal_dis2(int id1, int id2)
{
	double tmp_sum = 0 ;
	for(int i = 0; i < k ; i ++) tmp_sum  += pow((y[i][id1] - y[i][id2]) ,2);
	return tmp_sum;
}

void KM_Rand_Init()
{
    vector<int> seed;
    
    srand(time(NULL));
    for(int i = 0; i < num_cluster; i ++)
	{
		int tmp = rand() % P.size();
    	seed.push_back(tmp);
	} 
//    int tmp = rand() % P.size();
//    seed.push_back(tmp);
//
//    for(int i = 1; i < num_cluster; i ++)
//    {
//        double total_sum = 0;
//        vector<double> dis;
//        for(int j = 0; j < P.size(); j ++)
//        {
//            double min_dis = 999999999;
//            for(int l = 0; l < i; l++)
//            {
//                double tmp_dis = cal_dis2(j, seed[l]);
//                if(tmp_dis < min_dis) min_dis = tmp_dis;
//            }
//            total_sum += min_dis;
//            dis.push_back(min_dis);
//        }
//        srand(time(NULL));
//        double Rand_num =  ((rand() % 10 + 1) * total_sum) / 10;
//        int seed_num = 0;
//        while(Rand_num > 0)
//        {
//            Rand_num -= dis[seed_num];
//            seed_num ++;
//        }
//        bool flag = true;
//        while(flag)
//        {
//            vector<int>::iterator it = find(seed.begin(), seed.end(), seed_num);
//            if(it == seed.end() && seed_num < P.size() && seed_num >= 0) flag = false;
//            else if(seed_num >= P.size()) seed_num --;
//            else seed_num ++;
//        }
//        
//        seed.push_back(seed_num);
//    }
    for(int i = 0; i < num_cluster; i ++)
    {
        for(int j = 0; j < k ; j ++) C0[i].center[j] = y[j][seed[i]];
        C0[i].id.push_back(seed[i]);
        P[seed[i]].flag = i;
    }

}

void KMeans()
{
    for(int i = 0; i < C0.size(); i ++) C0[i].id.clear();
	for(int i = 0; i < P.size(); i ++)
    {
        int n_point = -1;
        double min_dis = 999999999;
        for(int j = 0; j < C0.size(); j ++)
        {
            double tmp_dis = 0 ;
			for(int l = 0; l < k ; l ++) tmp_dis  += pow((C0[j].center[l] - y[l][i]) ,2);
            if(tmp_dis < min_dis) 
            {
                min_dis = tmp_dis;
                n_point = j;
            }
        }
        if(n_point == -1) cout << "ERRROR" << endl;
        C0[n_point].id.push_back(i);
        P[i].flag = n_point;
    }
    for(int i = 0; i < C0.size() ; i ++)
    {
        C0[i].update_center();
    }
}

void Output(string file_name)
{
    ofstream fout(file_name);
	if(!fout)
	{
		cout<<"WRITE DATA ERROR!"<<endl;
		return;
	}
    for(int i = 0; i < P.size(); i ++) fout << P[i].flag << endl;
    fout.close();
}

int main()
{
	Init();
	Input(input_file);
	Generate();
	KM_Rand_Init();
	for(int i = 0; i < iteration; i ++)
	{
		KMeans();
		if(i % 100 == 0) cout << i << endl;
	} 
	Output(output_file);
}
