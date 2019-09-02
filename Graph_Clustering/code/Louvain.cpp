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

struct cluster
{
    vector<int> ele;
    int label;
    double in;
    double tot;
    double ki;
};

vector<int> Label;
vector<vector<int> > adj;
vector<int> degree;
int degree_sum;
vector<cluster> clu;
string input_file, output_file;
int num;
int aim_cluster;


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
        vector<int> tmp_vec;
		stringstream ss(buffer);
        string tmp;
        while(ss >> tmp)
        {
            if(tmp == "0") tmp_vec.push_back(0);
            else if(tmp == "1") tmp_vec.push_back(1);
            else
            {
                cout << "Input ERROR" << endl;
                return;
            }
        }
        adj.push_back(tmp_vec);
    }
}


void Init()
{
    aim_cluster = 5;
    input_file = "wisconsin_adj.txt";
    output_file = "Louvain_wisconsin_adj.txt";
    Input(input_file);
    num = adj.size();
    clu.resize(num);
    degree.resize(num,0);
    Label.resize(num);
    degree_sum = 0;
    for(int i = 0; i < num; i ++) Label[i] = i;
    for(int i = 0; i < num; i ++)
    {
        for(int j = 0; j < num; j ++)
        {
            degree[i] += adj[i][j];
            degree_sum += adj[i][j];
        }
    }
    degree_sum /= 2;
    for(int i = 0; i < num; i ++)
    {
        clu[i].ele.push_back(i);
        clu[i].label = i;
        clu[i].in = 0;
        clu[i].tot = degree[i];
        clu[i].ki = degree[i];
    }
}

double cal_delta_Q(int i, int c)
{
    double kiin = 0;
    for(int j = 0 ; j < clu[c].ele.size(); j ++)
    {
        kiin += adj[i][clu[c].ele[j]];
    }
    kiin *= 2;
    double delta_Q = 0;
    double tmp_0 = (clu[c].in + kiin) / degree_sum;
    double tmp_1 = pow(((clu[c].tot + clu[c].ki) / (2 * degree_sum)), 2);
    double tmp_2 = clu[c].in / (2 * degree_sum);
    double tmp_3 = pow((clu[c].tot / (2 * degree_sum)), 2);
    double tmp_4 = pow((clu[c].ki / (2 * degree_sum)), 2);
    delta_Q  = (tmp_0 - tmp_1) - (tmp_2 - tmp_3 - tmp_4);
    return delta_Q;
}

bool connected(int m, int n)
{
	bool conn = false;
	for(int i = 0; i < clu[m].ele.size(); i ++)
	{
		for(int j = 0; j < clu[n].ele.size(); j ++)
		{
			if(adj[clu[m].ele[i]][clu[n].ele[j]] == 1)
			{
				conn  = true;
				return conn;
			}
		}
	}
}

bool merge0()
{
    bool update = false;
    for(int i = 0; i < num; i ++)
    {
        double max_delta_Q = -1;
        double max_point;
        for(int j = 0 ; j < num; j ++)
        {
            if(i != j && connected(i, j))
            {
                double tmp_delta_Q = cal_delta_Q(i, j);
                if(tmp_delta_Q > max_delta_Q)
                {
                    max_delta_Q = tmp_delta_Q;
                    max_point = j;
                }
            }
        }
        if(max_delta_Q != -1)
        {
            update = true;
            for(int j = 0; j < clu[i].ele.size(); j ++)
            {
                Label[clu[i].ele[j]] = Label[clu[max_point].ele[0]];
                //cout<<"merge"<<i<<" "<<max_point<<endl;
            }
            
        }  
    }
    return update;
}

void Louvain()
{
    while(merge0())
    {
        cout<<"UPDATE"<<endl;
		clu.clear();
        vector<int> new_label;
        for(int i = 0; i < Label.size(); i ++)
        {
            int j = 0;
            while(j < new_label.size())
            {
            	if(new_label[j] == Label[i]) break;
            	j ++;
			}
            if(j == new_label.size())
            {
                new_label.push_back(Label[i]);
                cluster c;
                clu.push_back(c);
                clu.back().ele.push_back(i);
                clu.back().label = clu.size() - 1;
            }
            else
            {
				clu[j].ele.push_back(i);
            }
        }
        num = new_label.size();
        //cout<<num<<endl;
        for(int i = 0; i < clu.size(); i ++)
        {
            clu[i].in = 0;
            for(int j = 0 ; j < clu[i].ele.size(); j ++)
            {
                for(int k = j + 1 ; k < clu[i].ele.size(); k ++)
                {
                    clu[i].in += adj[clu[i].ele[j]][clu[i].ele[k]];
                }
            }
            clu[i].tot = 0;
            for(int j = 0; j < clu[i].ele.size(); j ++)
            {
                clu[i].tot += degree[clu[i].ele[j]];
            }
            clu[i].ki = clu[i].tot - 2 * clu[i].in;
            for(int j = 0; j < clu[i].ele.size(); j ++)
            {
                Label[clu[i].ele[j]] = i;
            }
        }
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
    for(int i = 0; i < Label.size(); i ++) fout << Label[i] + 1 << endl;
    fout.close();
}

int main()
{
	Init();
    Louvain();
    Output(output_file);
    return 0;
}


