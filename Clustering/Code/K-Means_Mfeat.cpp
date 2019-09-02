#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<time.h>
#include<algorithm>
#include<cmath>

using namespace std;

struct point
{
    vector<double> ch;
    int flag;
};

struct cluster
{
    vector<double> center;
    vector<point> V;
    void update_center()
    {
        for(int i = 0 ; i < 649; i ++)
        {
        	double sum;
        	for(int j = 0; j < V.size(); j ++)
	        {
	            sum += V[j].ch[i];
	        }
	        sum  /= 649;
	        center[i] = sum;
		}
    };
};

int iteration, num_cluster;
vector<point> V0;
vector<cluster> C0;
string input_file, output_file;

double cal_dis(point p1, point p2)
{
    double dis = 0;
    for(int i = 0; i < 649; i ++)
    {
    	dis += pow((p1.ch[i] - p2.ch[i]), 2);
	}
	dis = sqrt(dis);
    return dis;
}

void Init()
{
    iteration = 1000;
    num_cluster = 10;
    input_file = "datasets\\Mfeat.txt";
    output_file = "result\\KMeans_Mfeat.txt";
    C0.resize(num_cluster);
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
    int test_num = 0;
    while(getline(fin, buffer))
    {
        point p;
		stringstream ss(buffer);
        string tmp;
        for(int i = 0; i < 649; i ++)
        {
        	ss >> tmp;
        	double tmp_num = stod(tmp);
        	p.ch.push_back(tmp_num);
		}
        p.flag = -1;

        V0.push_back(p);
    }
}

void norm()
{
	double max,min;
	for(int i = 0; i < 649; i ++)
	{
		max = V0[0].ch[i];
		min = max;
		for(int j = 0 ; j < 2000; j ++)
		{
			if(V0[j].ch[i] < min) min = V0[j].ch[i];
			if(V0[j].ch[i] > max) max = V0[j].ch[i];
		}
		double tmp_num  = max - min;
		for(int j = 0 ; j < 2000; j ++)
		{
			V0[j].ch[i] = (V0[j].ch[i] - min) * 100 / tmp_num;
		}
	}
}

void KM_Rand_Init()
{
    vector<int> seed;
    
    srand(time(NULL));
    for(int i = 0; i < num_cluster; i ++)
	{
		int tmp = rand() % V0.size();
    	seed.push_back(tmp);
	} 
    for(int i = 0; i < num_cluster; i ++)
    {
        C0[i].center.resize(649);
		for(int j = 0; j < 649; j ++)
        {
        	C0[i].center[j] = V0[seed[i]].ch[j];
		}
        C0[i].V.push_back(V0[seed[i]]);
        V0[seed[i]].flag = i;
    }

}

void KMeans()
{
    for(int i = 0; i < C0.size() ; i ++)
    {
        C0[i].V.clear();
    }
	for(int i = 0; i < V0.size(); i ++)
    {
        int n_point = -1;
        double min_dis = 999999;
        for(int j = 0; j < C0.size(); j ++)
        {
            point center;
			center.ch.resize(649);
            for(int k = 0; k < 649; k++)
            {
            	center.ch[k] = C0[j].center[k];
			}
            double tmp_dis = cal_dis(center, V0[i]);
            if(tmp_dis < min_dis) 
            {
                min_dis = tmp_dis;
                n_point = j;
            }
        }
        if(n_point == -1) cout << "ERRROR" << endl;
        C0[n_point].V.push_back(V0[i]);
        V0[i].flag = n_point;
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
    for(int i = 0; i < V0.size(); i ++) fout << V0[i].flag << endl;
    fout.close();
}

int main()
{
	Init();
    Input(input_file);
    norm();
    KM_Rand_Init();
    for(int i = 0; i < iteration; i ++)
	{
		KMeans();
		if(i % 100 == 0) cout << i << endl;
	} 
    Output(output_file);
}


