/*
datasets\\Aggregation_cluster=7.txt
datasets\\flame_cluster=2.txt
datasets\\Jain_cluster=2.txt
datasets\\Pathbased_cluster=3.txt
datasets\\Spiral_cluster=3.txt
*/
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<time.h>
#include<algorithm>

using namespace std;

struct point
{
    double x;
    double y;
    int flag;
};

struct cluster
{
    double center_x;
    double center_y;
    vector<point> V;
    void update_center()
    {
        double sum_x = 0;
        double sum_y = 0;
        for(int i = 0; i < V.size(); i ++)
        {
            sum_x += V[i].x;
            sum_y += V[i].y;
        }
        sum_x /= V.size();
        sum_y /= V.size();
        center_x = sum_x;
        center_y = sum_y;
    };
};

int iteration, num_cluster;
double max_x, max_y, min_x, min_y;
vector<point> V0;
vector<cluster> C0;
string input_file, output_file;

double cal_dis(point p1, point p2)
{
    double x1 = p1.x;
    double x2 = p2.x;
    double y1 = p1.y;
    double y2 = p2.y;
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

void Init()
{
    max_x = 0;
    max_y = 0;
    min_x = 999999;
    min_y = 999999;
    iteration = 2000;
    num_cluster = 3;
    input_file = "datasets\\Spiral_cluster=3.txt";
    output_file = "result\\KMeans_Spiral_cluster=3.txt";
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
        ss >> tmp;
        double tmp_num = stod(tmp);

        p.x = tmp_num;
        if(tmp_num < min_x) min_x = tmp_num;
        if(tmp_num > max_x) max_x = tmp_num;
        ss >> tmp;
        tmp_num = stod(tmp);
        p.y = tmp_num;
        if(tmp_num < min_y) min_y = tmp_num;
        if(tmp_num > max_y) max_y = tmp_num;
        p.flag = -1;

        V0.push_back(p);
    }
}

void KM_Rand_Init()
{
    vector<int> seed;
    
    srand(time(NULL));
    int tmp = rand() % V0.size();
    seed.push_back(tmp);

    for(int i = 1; i < num_cluster; i ++)
    {
        double total_sum = 0;
        vector<double> dis;
        for(int j = 0; j < V0.size(); j ++)
        {
            double min_dis = 999999;
            for(int k = 0; k < i; k++)
            {
                double tmp_dis = cal_dis(V0[j], V0[seed[k]]);
                if(tmp_dis < min_dis) min_dis = tmp_dis;
            }
            total_sum += min_dis;
            dis.push_back(min_dis);
        }
        srand(time(NULL));
        double Rand_num =  ((rand() % 10 + 1) * total_sum) / 10;
        int seed_num = 0;
        while(Rand_num > 0)
        {
            Rand_num -= dis[seed_num];
            seed_num ++;
        }
        bool flag = true;
        while(flag)
        {
            vector<int>::iterator it = find(seed.begin(), seed.end(), seed_num);
            if(it != seed.end() && seed_num < V0.size() && seed_num >= 0) flag = false;
            else if(seed_num >= V0.size()) seed_num --;
            else seed_num ++;
        }
        
        seed.push_back(seed_num);
    }
    for(int i = 0; i < num_cluster; i ++)
    {
        C0[i].center_x = V0[seed[i]].x;
        C0[i].center_y = V0[seed[i]].y;
        C0[i].V.push_back(V0[seed[i]]);
        V0[seed[i]].flag = i;
    }

}

void KMeans()
{
    for(int i = 0; i < V0.size(); i ++)
    {
        int n_point = -1;
        double min_dis = 999999;
        for(int j = 0; j < C0.size(); j ++)
        {
            point center;
            center.x = C0[j].center_x;
            center.y = C0[j].center_y;
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
    KM_Rand_Init();
    for(int i = 0; i < iteration; i ++)
	{
		KMeans();
		if(i % 100 == 0) cout << i << endl;
	} 
    Output(output_file);
    return 0;
}


