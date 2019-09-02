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
#include<cmath>
#include<queue>
#include<sstream>

using namespace std;

struct point
{
    double x;
    double y;
    int id;
    int label;
    int flag; //-1未处理，0噪声，1核心，2边界
};

struct cluster
{
    vector<point> p_set;
};

double epsion, density;
vector<point> P;
vector<cluster> C;
string input_file, output_file;

bool check_core(point p)
{
    int cnt = 0;
    for(int i = 0; i < P.size(); i ++)
    {
        double tmp_dis = pow((p.x - P[i].x), 2) +  pow((p.y - P[i].y), 2);
        double tmp_ep = pow(epsion, 2);
        if(tmp_dis <= tmp_ep && tmp_dis > 0)
        {
            cnt ++;
        }
    }
    if(cnt >= density) return true;
    else return false;
}

void Init()
{
    input_file = "datasets\\Spiral_cluster=3.txt";
    output_file = "result\\DBSCAN_Spiral_cluster=3.txt";
    epsion  = 3.6;
    density = 4;
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
    //int test_num = 0;
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
        p.id = input_cnt;
        P.push_back(p);
        input_cnt ++;
    }
}

void DBSCAN()
{
    for(int i = 0 ; i < P.size(); i++)
    {
        cout << i << endl;
		if(P[i].flag == -1)
        {
            if(check_core(P[i]))
            {
                P[i].flag = 1;
                cluster tmp_c;
                tmp_c.p_set.push_back(P[i]);
                queue<point> Q;
                for(int j = 0; j < P.size(); j ++)
                {
                    double tmp_dis = pow((P[j].x - P[i].x), 2) +  pow((P[j].y - P[i].y), 2);
                    double tmp_ep = pow(epsion, 2);
                    if(tmp_dis <= tmp_ep && tmp_dis > 0)
                    {
                        Q.push(P[j]);
                        tmp_c.p_set.push_back(P[j]);
                        P[j].flag = 2;
                    }
                }
                while(! Q.empty())
                {
                    if(check_core(Q.front()))
                    {
                        P[Q.front().id].flag = 1;
                        for(int k = 0; k < P.size(); k ++)
                        {
                            if(P[k].flag == -1)
                            {
                                double tmp_dis = pow((P[k].x - Q.front().x), 2) +  pow((P[k].y - Q.front().y), 2);
                                double tmp_ep = pow(epsion, 2);
                                {
                                    if(tmp_dis <= tmp_ep && tmp_dis > 0)
                                    {
                                        Q.push(P[k]);
                                        tmp_c.p_set.push_back(P[k]);
                                        P[k].flag = 2;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        P[Q.front().id].flag = 2;
                    }
                    Q.pop();
                }
                C.push_back(tmp_c);
            }
            else P[i].flag = 2;
        }
    }
}

void Output(string file_name)
{
    for(int i = 0; i < C.size(); i ++)
    {
        for(int j = 0; j < C[i].p_set.size(); j ++)
        {
            P[C[i].p_set[j].id].label = i;
        }
    }
    ofstream fout(file_name);
    if(!fout)
	{
		cout<<"WRITE DATA ERROR!"<<endl;
		return;
	}
    for(int i = 0; i < P.size(); i ++)
	{
		if(P[i].flag != 2) fout << P[i].label << endl;
		else fout << -1 << endl;
	} 
    fout.close();
}

int main()
{
    Init();
    Input(input_file);
    cout <<  "Starting DBSCAN!" << endl;
    DBSCAN();
    Output(output_file);
    return 0;
}
