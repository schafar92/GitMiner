#include "square.h"
#include <time.h>

int main()
{
    int v;
    cout<<"Please enter the order v = ";
    cin>>v;
    clock_t start = clock();
    square s(v,2);
    s.non_uniform(0.5);
    clock_t stop = clock();
    cout<<(double)(stop-start)/CLOCKS_PER_SEC<<"s"<<endl;

}
square::square(int _v, int _d):v(_v),d(_d)
{
    v_fac = 1;
    gA = 1;
    gB = 1;
    int temp = v+2;
    do
    {
        v_fac *= temp;
        temp--;
    }while(temp > 1);
}
void square::non_uniform(double delta_miu)
{
    generate_hopping();
    calculate_hop();
}
void square::generate_hopping()
{
    int comp_res = 0;
    int d_weight;
    int num_sites = 0;
    vector<int> hopArray, hopArray_new;
    vector<int> :: iterator it_vect;
    //int Count = 0;
	//static omp_lock_t lock;
	//omp_init_lock(&lock);
//#pragma omp parallel for
    for(int num = 0; num < (int)pow((float)4,v); num++)
    {
        hopArray = find_hop(num);
        /*for(it_vect = hopArray.begin(); it_vect != hopArray.end(); it_vect++)
            cout<<*it_vect;
        cout<<endl;*/
        hopArray_new = smid(hopArray);
        comp_res = comp_diagram(hopArray, hopArray_new);
        if(comp_res == 0)
        {
            /*for(it_vect = hopArray.begin(); it_vect != hopArray.end(); it_vect++)
                cout<<*it_vect;
            cout<<endl;*/
            hopArray_new = cfour(hopArray, d_weight);
            comp_res = comp_diagram(hopArray, hopArray_new);
            //cout<<comp_res<<endl;
            if(comp_res == 0)
            {
                num_sites = get_sitenum(hopArray);
                //cout<<num_sites<<endl;

                d_weight = 8 / d_weight;
                //Count++;
                //omp_set_lock(&lock);
                adiagram(hopArray, d_weight, num_sites);
                //omp_unset_lock(&lock);
            }
        }

    }
    //omp_destroy_lock(&lock);
    //cout<<Count<<endl;
}
vector<int> square::find_hop(int num)
{
    vector<int> hopArray(v,0);
    for(int i = 0; i < v; i++)
    {
        if(num < 2*d)
        {
            hopArray.at(i) = num;
            break;
        }
        else
        {
            hopArray.at(i) = num % (2*d);
            num = (num - hopArray.at(i)) / (2*d);
        }
    }
    return hopArray;
}
vector<int> square::smid(vector<int> hopArray_in)
{
    vector<int> hopArray_out;
    int ***arrow_sites = new int**[2*v+1];
    for(int i = 0; i < 2*v+1; i++)
    {
        arrow_sites[i] = new int*[2*v+1];
        for(int j = 0; j < 2*v+1; j++)
            arrow_sites[i][j] = new int[2*d];
    }
    for(int i = 0; i < 2*v+1; i++)
        for(int j = 0; j < 2*v+1; j++)
            for(int k = 0; k < 2*d; k++)
                arrow_sites[i][j][k] = 0;
    int x = v;
    int y = v;
    for(int i = 0; i < v; i++)
    {
        if(hopArray_in[i] == 0)
        {
            arrow_sites[x][y][0] += 1;
            x += 1;
        }
        else if(hopArray_in[i] == 1)
        {
            arrow_sites[x][y][1] += 1;
            y -= 1;
        }
        else if(hopArray_in[i] == 2)
        {
            arrow_sites[x][y][2] += 1;
            x -= 1;
        }
        else if(hopArray_in[i] == 3)
        {
            arrow_sites[x][y][3] += 1;
            y += 1;
        }
    }
    hopArray_out = find_small(arrow_sites);
    for(int i = 0; i < 2*v+1; i++)
    {
        for(int j = 0; j < 2*v+1; j++)
            delete [] arrow_sites[i][j];
        delete [] arrow_sites[i];
    }
    delete [] arrow_sites;
    return hopArray_out;
}
vector<int> square::find_small(int ***arrow_sites)
{
    bool jump_cycle;
    bool judge_end;
    vector<int> hopArray_out(v, 0);
    int l = 0;
    int x = v;
    int y = v;
    int n = 0;

    int **c = new int*[v];
	for(int i = 0; i < v; i++)
		c[i] = new int[2*d];
    for(int i = 0; i < v; i++)
        for(int j = 0; j < 2*d; j++)
            c[i][j] = -1;
/*
    for(int i = 0; i < 2*d; i++)
        cout<<arrow_sites[v][v][i];
    cout<<endl;
    system("pause");
*/
    jump_cycle = false;
    do{
        n += 1;
        if(n > v)
        {
            cout<<"wrong with the hopping loop!"<<endl;
            break;
        }
        if(c[n-1][0] == -1)
        {
            l = 0;
            for(int k = 0; k < 2*d; k++)
            {

                if(arrow_sites[x][y][k] > 0)
                {
                    c[n-1][l] = k;
                    l += 1;
                }
            }
        }
        //--------------------------------------
        arrow_sites[x][y][c[n-1][0]] -= 1;
        hopArray_out[n-1] = c[n-1][0];
        //cout<<a_out[n-1];
        //system("pause");
        //cout<<"("<<x<<","<<y<<")->";
        if(c[n-1][0] == 0)
            x += 1;
        else if(c[n-1][0] == 1)
            y -= 1;
        else if(c[n-1][0] == 2)
            x -= 1;
        else if(c[n-1][0] == 3)
            y += 1;
        //cout<<"("<<x<<","<<y<<")->";
        //system("pause");
        //---------------------------------------
        //---------------------------------------
        /*for(int i = 0; i < l; i++)
            cout<<c[i][n-1];
        cout<<endl;
        system("pause");*/
        //if(x == v && y == v)
        //{
            judge_end = true;
            for(int k = 0; k < 2*d; k++)
            {
                if(arrow_sites[x][y][k] > 0)
                {
                    judge_end = false;
                    break;
                }
            }
            if(judge_end)
            {
                if( n == v )
                    jump_cycle = true;
                else
                {

                    //while是先判断条件在进入循环体。而do-while则是先进入循环体，然后结束一次循环后进行判断！
                    while(c[n-1][1] == -1)
                    {
                        hopArray_out[n-1] = -1;
                        if(c[n-1][0] == 0)
                            x -= 1;
                        else if(c[n-1][0] == 1)
                            y += 1;
                        else if(c[n-1][0] == 2)
                            x += 1;
                        else if(c[n-1][0] == 3)
                            y -= 1;
                        arrow_sites[x][y][c[n-1][0]] += 1;
                        c[n-1][0] = -1;
                        n -= 1;
                    }
                    hopArray_out[n-1] = -1;
                    if(c[n-1][0] == 0)
                        x -= 1;
                    else if(c[n-1][0] == 1)
                        y += 1;
                    else if(c[n-1][0] == 2)
                        x += 1;
                    else if(c[n-1][0] == 3)
                        y -= 1;
                    arrow_sites[x][y][c[n-1][0]] += 1;
                    c[n-1][0] = c[n-1][1];
                    c[n-1][1] = c[n-1][2];
                    c[n-1][2] = c[n-1][3];
                    c[n-1][3] = -1;
                    n -= 1;

                }
        }
    }while(jump_cycle == false);
    //system("pause");
    for(int i = 0; i < v; i++)
    {
        delete [] c[i];
    }
    delete [] c;

    return hopArray_out;
}
int square::comp_diagram(vector<int> hopArray, vector<int> hopArray_new)
{
    int comp_res;
    for(int i = 0; i < v; i++)
    {
        if(hopArray[i] > hopArray_new[i])
        {
            comp_res = 1;
            break;
        }
        else if(hopArray[i] < hopArray_new[i])
        {
            comp_res = -1;
            break;
        }
        else
            comp_res = 0;
    }
    return comp_res;
}
vector<int> square::cfour(vector<int> hopArray, int &d_weight)
{
    vector<int> hopArray_mid(v,0);
    vector<int> hopArray_new(v,0);
    vector<int> hopArray_out(v,0);
    int comp_res1 = 0;
    int comp_res2 = 0;
    for(int i = 0; i < v; i++)
    {
        hopArray_out[i] = hopArray[i];
        //cout<<a_out[i];
    }
    //cout<<endl;
    //system("pause");
    d_weight = 1;


   //========================找对称性===========================
   //90°，180°，270°，x=0, y=0; y = x, y = -x 七种
   for(int ios = 1; ios <= 7; ios++)
   {
       hopArray_mid = sym_oper(hopArray, ios);
/*
       for(int i = 0; i < v; i++)
            cout<<a_mid[i];
       cout<<endl;
*/
       hopArray_new = smid(hopArray_mid);
       /*for(int i = 0; i < v; i++)
            cout<<a_new1[i];
       cout<<endl;
       system("pause");*/

       comp_res1 = comp_diagram(hopArray, hopArray_new);
       //cout<<comp_res1<<endl;
       if(comp_res1 == 0)
       {
           d_weight += 1;
       }
       else if(comp_res1 > 0)
       {
            comp_res2 = comp_diagram(hopArray_out, hopArray_new);

            //cout<<"comp_res2 = "<<comp_res2<<endl;
            //system("pause");

            if(comp_res2 == 1)
            {
                hopArray_out = hopArray_new;
            }
       }

   }
   return hopArray_out;

}
vector<int> square::sym_oper(vector<int> hopArray, int ios)
{
    vector<int> hopArray_out(hopArray);

    if(ios == 1)   //关于x对称
    {
        for(int i = 0; i < v; i++)
        {
            if(hopArray[i] == 1)
                hopArray_out[i] = 3;
            else if(hopArray[i] == 3)
                hopArray_out[i] = 1;
        }
    }
    else if(ios == 2)//关于y对称
    {
        for(int i = 0; i < v; i++)
        {
            if(hopArray[i] == 0)
                hopArray_out[i] = 2;
            else if(hopArray[i] == 2)
                hopArray_out[i] = 0;
        }
    }
    else if(ios == 3)//关于y=x对称
    {
        for(int i = 0; i < v; i++)
        {
            if(hopArray[i] == 0)
                hopArray_out[i] = 3;
            else if(hopArray[i] == 1)
                hopArray_out[i] = 2;
            else if(hopArray[i] == 2)
                hopArray_out[i] = 1;
            else if(hopArray[i] == 3)
                hopArray_out[i] = 0;
        }
    }
    else if(ios == 4)//关于y=-x对称
    {
        for(int i = 0; i < v; i++)
        {
            if(hopArray[i] == 0)
                hopArray_out[i] = 1;
            else if(hopArray[i] == 1)
                hopArray_out[i] = 0;
            else if(hopArray[i] == 2)
                hopArray_out[i] = 3;
            else if(hopArray[i] == 3)
                hopArray_out[i] = 2;
        }
    }
    else if(ios == 5)//90°旋转
    {
        for(int i = 0; i < v; i++)
        {
            if(hopArray[i] == 0)
                hopArray_out[i] = 3;
            else if(hopArray[i] == 1)
                hopArray_out[i] = 0;
            else if(hopArray[i] == 2)
                hopArray_out[i] = 1;
            else if(hopArray[i] == 3)
                hopArray_out[i] = 2;
        }
    }
    else if(ios == 6)//180°旋转
    {
        for(int i = 0; i < v; i++)
        {
            if(hopArray[i] == 0)
                hopArray_out[i] = 2;
            else if(hopArray[i] == 1)
                hopArray_out[i] = 3;
            else if(hopArray[i] == 2)
                hopArray_out[i] = 0;
            else if(hopArray[i] == 3)
                hopArray_out[i] = 1;
        }
    }
    else if(ios == 7)//270°旋转
    {
        for(int i = 0; i < v; i++)
        {
            if(hopArray[i] == 0)
                hopArray_out[i] = 1;
            else if(hopArray[i] == 1)
                hopArray_out[i] = 2;
            else if(hopArray[i] == 2)
                hopArray_out[i] = 3;
            else if(hopArray[i] == 3)
                hopArray_out[i] = 0;
        }
    }
    else
        cout<<"sym_oper: Error! The ios must between 1 and 7!\n";
    return hopArray_out;
}
int square::get_sitenum(vector<int> hopArray)
{
    int **site_coordinate = new int*[v+1];//记录跃迁从头到尾的格点坐标
	for(int i = 0; i < v+1; i++)
		site_coordinate[i] = new int[2];
    bool same_site;
    int num_sites_temp;
    int x = v;
    int y = v;
    //跃迁坐标
    site_coordinate[0][0] = x;
    site_coordinate[0][1] = y;
    num_sites_temp = 0;
    for(int i = 0; i < v; i++)
    {
        if(hopArray[i] == 0)
            x += 1;
        else if(hopArray[i] == 1)
            y -= 1;
        else if(hopArray[i] == 2)
            x -= 1;
        else if(hopArray[i] == 3)
            y += 1;

        same_site = false;
        for(int j = 0; j <= num_sites_temp; j++)
        {
            if(site_coordinate[j][0] == x && site_coordinate[j][1] == y )
            {
                same_site = true;
                break;
            }
        }
        //cout<<same_site<<endl;
        if(same_site == false)
        {
            num_sites_temp += 1;
            site_coordinate[num_sites_temp][0] = x;//储存的是不同的格点的位置信息,v+1个格点
            site_coordinate[num_sites_temp][1] = y;
        }
    }
    /*cout<<num_sites_temp<<endl;
    system("pause");*/

    return num_sites_temp;
}
void square::adiagram(vector<int> hopArray,int d_weight, int num_sites)
{
    diaList.push_back(hopArray);
    weight.push_back(d_weight);
    site.push_back(num_sites);
}
void square::calculate_hop()
{
//=======================================================
    //目的是把已经生成的kato.txt以及NumberOfKatoterm.txt数据读出
    char filename1[50], filename2[50];
    sprintf(filename1,"NumberOfkatoterms%d.txt",v+2);
    sprintf(filename2,"kato%d.txt",v+2);
    ifstream fin1,fin2;
    fin1.open(filename1);
    fin2.open(filename2);
    int num;
    fin1>>num;

    int *eng_weight = new int[num];//这个部分我就不改成stl了。感觉麻烦and没有必要
    int **eng_katolist = new int*[num];
    for(int i = 0; i < num; i++)
        eng_katolist[i] = new int[v+1];

    if(!fin2.is_open())
    {
        cerr<<"Something wrong during open kato.txt!"<<endl;
        exit(0);
    }
    int katoTerm;
    int countNum = 1;
    int j = 0;
    int i = 0;
    while(fin2>>katoTerm)
    {
        if(countNum % (v+2) != 0)
        {
            eng_katolist[i][j] = katoTerm;
            //cout<<eng_katolist[i][j];
            countNum++;
            j++;
        }
        else
        {
            eng_weight[i] = katoTerm;
            //cout<<" "<<eng_weight[j]<<endl;
            countNum++;
            j = 0;
            i++;
        }
    }
//=========================================================
    //计算跃迁过程部分
    vector<int> hopArray(v, 0);
    int site_num;
    int num_process = 0;
    int d_weight;
    int t_weight;
    list<vector<int> >   topo_list;
    list<vector<int> >::iterator it_list;
    vector<int>::iterator it_vect;
    vector<int>::iterator it_weight;
    vector<int>::iterator it_site;
    int **oper;
    oper = new int*[v];
    for(int i = 0; i < v; i++)
        oper[i] = new int[2*d];
    //arrow中存的是跃迁过程中格点的标记。000222过程为1234321
    int **arrow = new int*[v];
    for(int i = 0; i < v; i++)
        arrow[i] = new int[2];

    for(it_list = diaList.begin(), it_weight = weight.begin(), it_site = site.begin(); it_list != diaList.end(), it_weight != weight.end(), it_site != site.end(); \
        it_list++, it_weight++, it_site++)
    {
        hopArray = *it_list;
        d_weight = *it_weight;
        site_num = *it_site + 1;

        oper = mapping_oper(hopArray, oper);

        arrow = find_topo(oper, arrow, site_num);
/*
        for(int i = 0; i < v; i++)
            cout<<arrow[i][0];
        cout<<arrow[v-1][1];
        cout<<endl;
*/
        topo_list = get_topo_list(arrow, topo_list, d_weight);
    }
/*
    for(it_list = topo_list.begin(); it_list != topo_list.end(); it_list++)
    {
        for(it_vect = (*it_list).begin(); it_vect != (*it_list).end(); it_vect++)
            cout<<*it_vect<<" ";
        cout<<endl;
    }
    cout<<endl;
    system("pause");
*/
//---------------------------------------------------------
    double *energy_diagram = new double[2];
    double *eOrder_V = new double[2];
    for(int i = 0; i < 2; i++)
    {
        energy_diagram[i] = 0.0;
        eOrder_V[i] = 0.0;
    }
    int ***oper_collection;
    oper_collection = new int **[v+2];
    for(int i = 0; i < v+2; i++)
    {
        oper_collection[i] = new int *[2*d];
        for(int j = 0; j < 2*d; j++)
        {
            oper_collection[i][j] = new int [v_fac];
        }
    }
//---------------------------------------------------------
    ofstream energyA, energyB;
    char name1[50];
    char name2[50];
    double delta_miu = 0.5;
    sprintf(name1,"alphaA_sl_%d_deltamiu=%f.txt", v,delta_miu);
    sprintf(name2,"alphaB_sl_%d_deltamiu=%f.txt", v,delta_miu);
    energyA.open(name1);
    energyB.open(name2);
    list<vector<int> >   topo_list_temp;
    vector<list<vector<int> > >   topo_temp(49,topo_list);
    //static omp_lock_t lock;
    //omp_init_lock(&lock);

//#pragma omp parallel for private(arrow), private(oper), private(oper_collection), private(energy_diagram), private(eOrder_V)
    for(int count = 1; count < 50; count++){
//---------------------------------------------------------
    vector<int> hopArray(v, 0);
    int site_num;
    int num_process = 0;
    int d_weight;
    int t_weight;
    list<vector<int> >   topo_list;
    list<vector<int> >::iterator it_list;
    vector<int>::iterator it_vect;
    vector<int>::iterator it_weight;
    vector<int>::iterator it_site;

    double miu = 0.0;
    miu = double(count)/100;
        topo_list_temp = topo_temp[count-1];
    for(it_list = diaList.begin(), it_weight = weight.begin(), it_site = site.begin(); it_list != diaList.end(), it_weight != weight.end(), it_site != site.end(); \
        it_list++, it_weight++, it_site++)
    {
        hopArray = *it_list;
        d_weight = *it_weight;
        site_num = *it_site + 1;

        oper = mapping_oper(hopArray, oper);

        arrow = find_topo(oper, arrow, site_num);

        topo_temp[count-1] = map_topo_list(arrow, topo_temp[count-1], t_weight);
        //cout<<t_weight<<endl;
        if(t_weight != 0)
        {
            //cout<<t_weight<<endl;
            oper_collection = platoon(oper_collection, oper, num_process);
            //cout<<num_process<<endl;
            //miu = 0.10;
            energy_diagram = cal_eng_dia(oper_collection, num, eng_weight, eng_katolist, num_process, energy_diagram, miu, delta_miu);
            //cout<<energy_diagram[0]<<" "<<energy_diagram[1]<<endl;
            energy_diagram[0] *= t_weight;
            energy_diagram[1] *= t_weight;
            //cout<<miu<<" "<<energy_diagram<<endl;//system("pause");
            eOrder_V[0] += energy_diagram[0];
            eOrder_V[1] += energy_diagram[1];
        }
    }
        //cout<<fixed<<setprecision(6)<<eOrder_V[0]<<endl;system("pause");
        //omp_set_lock(&lock);
        energyA<<miu<<" "<<fixed<<setprecision(6)<<eOrder_V[0]<<endl;
        energyB<<miu<<" "<<fixed<<setprecision(6)<<eOrder_V[1]<<endl;
        //omp_unset_lock(&lock);
        topo_temp[count-1] = topo_list_temp;
        eOrder_V[0] = 0.0;
        eOrder_V[1] = 0.0;

//------------------------------------------------------------------------------
        delete [] energy_diagram;
        delete [] eOrder_V;

        for(int i = 0; i < v; i++)
            delete [] oper[i];
        delete [] oper;

        for(int i = 0; i < v; i++)
        {
            for(int j = 0; j < 2*d; j++)
            {
                delete [] oper_collection[i][j];
            }
            delete [] oper_collection[i];
        }
        delete [] oper_collection;

        for(int i = 0; i < v; i++)
        {
            delete [] arrow[i];
        }
        delete [] arrow;
//------------------------------------------------------------------------------
    }
    //omp_destroy_lock(&lock);
//===============================================================
    for(int i = 0; i < v; i++)
        delete [] oper[i];
    delete [] oper;

    for(int i = 0; i < v; i++)
    {
        delete [] arrow[i];
    }
    delete [] arrow;
//=================================================================
}
int** square::mapping_oper(vector<int> hopArray, int **oper_out)
{
    for(int i = 0; i < v; i++)
        for(int j = 0; j < 2*d; j++)
            oper_out[i][j] = 0;

    for(int i = 0; i < v; i++)
    {
        //2d维表示的是从某点(00)到某点(01)的坐标变化。例，0102表示点从01->02
        if(hopArray[i] == 0)
        {
            oper_out[i][2] = oper_out[i][0] + 1;
            oper_out[i][3] = oper_out[i][1];
        }
        else if(hopArray[i] == 1)
        {
            oper_out[i][2] = oper_out[i][0];
            oper_out[i][3] = oper_out[i][1] - 1;
        }
        else if(hopArray[i] == 2)
        {
            oper_out[i][2] = oper_out[i][0] - 1;
            oper_out[i][3] = oper_out[i][1];
        }
        else if(hopArray[i] == 3)
        {
            oper_out[i][2] = oper_out[i][0];
            oper_out[i][3] = oper_out[i][1] + 1;
        }
        if(i < v-1)
        {
            oper_out[i+1][0] = oper_out[i][2];
            oper_out[i+1][1] = oper_out[i][3];
        }
    }
    return oper_out;
}
int **square::find_topo(int **oper, int **arrow, int site_num)
{
    bool judge_same;
    int numm;
    int **site_mapping = new int*[site_num];
    for(int i = 0; i < site_num; i++)
        site_mapping[i] = new int[2];
    arrow[0][0] = 1;
    arrow[0][1] = 2;
    numm = 2;
    site_mapping[0][0] = oper[0][0];
    site_mapping[0][1] = oper[0][1];
    site_mapping[1][0] = oper[0][2];
    site_mapping[1][1] = oper[0][3];
    //cout<<site_num<<endl;
    //cout<<oper[2][1]<<oper[3][1]<<endl;
    //system("pause");
    for(int i = 1; i < v; i++)
    {
        arrow[i][0] = arrow[i-1][1];
        judge_same = false;
        for(int j = 0; j < numm; j++)
        {
            if(oper[i][2] == site_mapping[j][0] && oper[i][3] == site_mapping[j][1])
            {
                arrow[i][1] = j+1;
                judge_same = true;
                break;
            }
        }
        if(judge_same == false)
        {
            numm += 1;
            arrow[i][1] = numm;
            site_mapping[numm-1][0] = oper[i][2];
            site_mapping[numm-1][1] = oper[i][3];
        }
    }
    //cout<<numm<<endl;
    //------------------------------------------------------------------------------------
/*
    for(int i = 0; i < v; i++)
        cout<<arrow[i][0]<<" ";
    cout<<arrow[v-1][1];
    cout<<endl;
    system("pause");
*/
    //---------------------
    for(int i = 0; i < site_num; i++)
    {
        delete [] site_mapping[i];
    }
    delete [] site_mapping;
    //----------------------
    return arrow;
}
list<vector<int> > square::get_topo_list(int **arrow, list<vector<int> > topo_list, int d_weight)
{
    list<vector<int> >::iterator it_list;
    vector<int>::iterator it_vec1;
    vector<int>::iterator it_vec2;
    vector<int> vec;
    bool judge_same1;
    judge_same1 = true;
    for(int i = 0; i < v; i++)
    {
        vec.push_back(arrow[i][0]);
    }
    vec.push_back(arrow[v-1][1]);
/*
    for(int i = 0; i < v+1; i++)
        cout<<vec.at(i);
    cout<<endl;
*/
    if(topo_list.size() == 0)
    {
        vec.push_back(d_weight);//weight
        vec.push_back(1);//是否被使用过的标记
        topo_list.push_back(vec);
    }
    else
    {
        for(it_list = topo_list.begin(); it_list != topo_list.end(); it_list++)
        {
            for(it_vec1 = (*it_list).begin(), it_vec2 = vec.begin(); it_vec1 != (*it_list).end(), it_vec2 != vec.end(); it_vec1++, it_vec2++)
            {
                if(*it_vec1 != *it_vec2)
                {
                    judge_same1 = false;
                    break;
                }
                judge_same1 = true;
            }
            if(judge_same1 == true)
            {
                it_vec1 = (*it_list).end() - 2;//获得d_weight的位置
                *it_vec1 += d_weight;
                break;
            }
        }
        if(judge_same1 == false)
        {
            vec.push_back(d_weight);
            vec.push_back(1);
            topo_list.push_back(vec);
        }
    }
    return topo_list;
}
list<vector<int> > square::map_topo_list(int **arrow, list<vector<int> > topo_list, int &t_weight)
{
    list<vector<int> >::iterator it_list;
    vector<int>::iterator it_vec1;
    vector<int>::iterator it_vec2;
    vector<int> vec;
    for(int i = 0; i < v; i++)
    {
        vec.push_back(arrow[i][0]);
    }
    vec.push_back(arrow[v-1][1]);
    /*for(it_vec1 = vec.begin(); it_vec1 != vec.end(); it_vec1++)
        cout<<*it_vec1;
    cout<<endl;*/
    bool judge_same;
    for(it_list = topo_list.begin(); it_list != topo_list.end(); it_list++)
    {
        for(it_vec1 = (*it_list).begin(), it_vec2 = vec.begin(); it_vec1 != (*it_list).end()-2, it_vec2 != vec.end(); it_vec1++, it_vec2++)
        {
            if(*it_vec2 != *it_vec1)
            {
                judge_same = false;
                break;
            }
            judge_same = true;
        }
        if(judge_same)
        {
            /*for(it_vec1 = (*it_list).begin(); it_vec1 != (*it_list).end(); it_vec1++)
                cout<<*it_vec1;
            cout<<endl;*/
            it_vec1 = (*it_list).end() - 1;

            if(*it_vec1 == 1)
            {
                t_weight = *(it_vec1-1);
                *it_vec1 = *it_vec1 - 1;
                break;
            }
            else
            {
                t_weight = 0;
                break;
            }

        }
    }
    /*cout<<t_weight<<endl;
    //system("pause");
    for(it_list = topo_list.begin(); it_list != topo_list.end(); it_list++)
    {
        for(it_vec1 = (*it_list).begin(); it_vec1 != (*it_list).end(); it_vec1++)
        {
            cout<<*it_vec1;
        }
        cout<<endl;
    }
    cout<<endl;
    system("pause");*/
    return topo_list;
}
int*** square::platoon(int ***oper_collection, int **oper, int &num_process)
{
   //cout<<v_fac<<endl;
   //system("pause");
   //如果是用有效场的算法，需要加入外源项，外源项当然是在全排列算法内。
   //我这里(1000,1000)->(某点)表示为外源流入，(某点)->(1000,1000)表示外源流出
    oper_collection[0][0][0] = 1000;
    oper_collection[0][1][0] = 1000;
    oper_collection[0][2][0] = oper[0][0];
    oper_collection[0][3][0] = oper[0][1];
    oper_collection[v+1][0][0] = oper[v-1][2];
    oper_collection[v+1][1][0] = oper[v-1][3];
    oper_collection[v+1][2][0] = 1000;
    oper_collection[v+1][3][0] = 1000;
   for(int i = 1; i < v+1; i++)
   {
       for(int j = 0; j < 2*d; j++)
       {
           for(int k = 0; k < v_fac; k++)
           {
                if(k == 0)
                    oper_collection[i][j][k] = oper[i-1][j];
                else
                    oper_collection[i][j][k] = 0;
           }
       }
   }
/*
   for(int i = 0; i < 2*d; i++)
   {
       for(int j = 0; j < v; j++)
       {
            cout<<oper_collection[i][j][0]<<" ";
       }
       cout<<endl;
   }
   cout<<endl;
   system("pause");
*/
   int *m = new int[v+2];
   int **oper_copy = new  int*[v+2];
   for(int i = 0; i < v+2; i++)
        oper_copy[i] = new int[2*d];
   bool judge_swap = true;
   for(int i =0 ; i < v+2; i++)
        m[i] = 0;
   m[0] = 0;
   //一个全排列算法。没懂。
   for(int i = 0; i < v+1; i++)
   {
       if(i < v+1)
            m[i+1] = m[i];
       for(int k = 0; k <= m[i]; k++)
       {
           for(int row = 0; row < v+2; row++)
           {
               for(int col = 0; col < 2*d; col++)
               {
                   oper_copy[row][col] = oper_collection[row][col][k];
                   //cout<<oper_collection[row][col][0]<<" ";
               }
               //cout<<endl;
           }
           //cout<<endl;
           for(int j = i+1; j < v+2; j++)
           {

               judge_swap = jud_swap(oper_copy, i, j, judge_swap);
               if(judge_swap)
               {
                  oper_copy = Swap(oper_copy, i, j);
                  m[i+1] += 1;
                  //cout<<m[i+1]<<"\t";
                  for(int row = 0; row < v+2; row++)
                      for(int col = 0; col < 2*d; col++)
                          oper_collection[row][col][m[i+1]] = oper_copy[row][col];
               }
           }
       }
   }
   num_process = m[v+1];
   //cout<<"num_process: "<<num_process<<endl;
   //system("pause");
   /*
    for(int k = 0; k < num_process+1; k++)
    {
        for(int i = 0; i < 2*d; i++)
        {
            for(int j = 0; j < v; j++)
            {
                cout<<oper_collection[i][j][k]<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
    }
    system("pause");
    */
   //------------------------------------
   delete [] m;
   for(int i = 0; i < v+2; i++)
        delete [] oper_copy[i];
   delete oper_copy;
   //------------------------------------

   return oper_collection;
}
bool square::jud_swap(int **oper, int m, int i, bool judge_swap)
{
/*    for(int i = 0; i < 2*d; i++)
    {
        for(int j = 0; j < v; j++)
        {
            cout<<oper[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    //cin.get();*/
    bool kkkk;
    judge_swap = true;
    for(int j = m; j <= i - 1; j++)
    {
        kkkk = true;
        //当前跃迁如果与之前的跃迁不同则交换
        for(int k = 0; k < 2*d; k++)
        {
            if(oper[i][k] != oper[j][k])
            {
                kkkk = false;
                break;
            }
        }
        if(kkkk)
            judge_swap = false;

    }
    //cout<<judge_swap<<endl;
    return judge_swap;
}
int** square::Swap(int **oper, int i, int j)
{
    int temp;
    for(int col = 0; col < 2*d; col++)
    {
        temp = oper[i][col];
        oper[i][col] = oper[j][col];
        oper[j][col] = temp;

    }
    return oper;
}
double* square::cal_eng_dia(int ***oper_collection, int m, int *eng_weight, int **eng_katolist, int num_process, double *energy_diagram, double miu, double delta_miu)
{
    int sites_num;
    bool same_site, judge_ground_A, judge_ground_B, judge_good_kato_A, judge_good_kato_B;
    int *ground_excited_A = new int[v+1];
    int *ground_excited_B = new int[v+1];
    double *wA = new double[v+1];//以A为出发点
    double *wB = new double[v+1];//以B为出发点
    double *h = new double[v+2];
    int **site = new int*[v+1];
    for(int i = 0; i < v+1; i++)//v+1个格点
        site[i] = new int[2];
    int **particle_num_A = new int*[v+3];//格点最多是v+1个，过程的话加上外源的流入流出和起始态总共应该是v+3
    for(int i = 0; i < v+3; i++)
        particle_num_A[i] = new int[v+1];
    int **particle_num_B = new int*[v+3];
    for(int i = 0; i < v+3; i++)
        particle_num_B[i] = new int[v+1];
    double energy_kato_A;
    double energy_kato_B;
    int num_zero1A, num_zero1B, num_zero2;
    energy_diagram[0] = 0.0;
    energy_diagram[1] = 0.0;
    for(int k = 0; k <= num_process; k++)
    {
        //site[0][0] = oper_collection[0][0][k];         //site是点信息，由数组的两位组成坐标，第二维则为第几次跃迁
        //site[0][1] = oper_collection[0][1][k];
        sites_num = 0;                                    //sites_num从0开始
        for(int i = 0; i < v+2; i++)
        {
            if(oper_collection[i][0][k] != 1000)
            {

                if(sites_num == 0)
                {
                    for(int j = 0; j < 2; j++)
                    {
                        site[sites_num][j] = oper_collection[i][j][k];
                    }
                    sites_num += 1;
                }
                else
                {
                    same_site = false;
                    for(int j = 0; j < i; j++)
                    {
                        //当到i位置时，遍历i前面的所有点，是否有跟i位置一样，有一样的话就是相同的site，不记录在site[][]中
                        //这里注意到，对oper_collection只针对前第一个点判断，这是因为第二个点将是下一维的第一个点，所以。。。
                        if(oper_collection[i][0][k] == oper_collection[j][0][k] && oper_collection[i][1][k] == oper_collection[j][1][k])
                        {
                            same_site = true;
                            break;// 这个break很重要!!!
                        }
                    }
                    if(same_site == false)
                    {
                        for(int j = 0; j < 2; j++)
                        {
                            site[sites_num][j] = oper_collection[i][j][k];
                        }
                        sites_num += 1;
                    }
                }
            }
        }
/*
    if(sites_num == 3)
    {
        cout<<sites_num<<endl;
        for(int i = 0; i < sites_num; i++)
        {
            cout<<site[i][0]<<site[i][1]<<" ";
        }
        cout<<endl;
        system("pause");
        for(int i = 0; i < v+2; i++)
        {
            cout<<"("<<oper_collection[i][0][k]<<","<<oper_collection[i][1][k]<<")"<<"->"<<"("<<oper_collection[i][2][k]<<","<<oper_collection[i][3][k]<<")"<<" ";
        }
        cout<<endl;
    }
*/

        for(int i = 0; i < v+3; i++)
        {
            for(int j = 0; j < sites_num; j++)
            {

                if(i == 0)//第一维为初始填充数
                {
                    if(abs(site[j][0] + site[j][1])%2 == 0)
                    {
                        particle_num_A[i][j] = gA;
                        particle_num_B[i][j] = gB;
                    }

                    else if(abs(site[j][0] + site[j][1])%2 == 1)
                    {
                        particle_num_A[i][j] = gB;
                        particle_num_B[i][j] = gA;
                    }
                }
                else
                {
                    particle_num_A[i][j] = 0;
                    particle_num_B[i][j] = 0;
                }

            }
        }
        //对于一个确定的i和k，oper_collection[][i][k]的前两个是跃迁起始点，后两位是跃迁终止点
        //site中第一维的两位是坐标，第二维是点的数目。保存着所有跃迁过程的点，在做sites_num循环
        //时，如果遍历到跃迁起始点，则particle_num-1，如果是跃迁终止点，则+1。这样，particle_num
        //数组中就存着各个过程某点的粒子数了。在起始位置遇到真空态(1000,1000)，第一个格点粒子加1
        //在结束位置则在末格点减1
        for(int i = 0; i < v+2; i++)
        {
            for(int j = 0; j < v+1; j++)
            {
                particle_num_A[i+1][j] = particle_num_A[i][j];
                particle_num_B[i+1][j] = particle_num_B[i][j];
            }

            //这里的sites_num是确切的格点数目，没有小1
            //cout<<oper_collection[i][0][k]<<endl;
            for(int j = 0; j < sites_num; j++)
            {
                if(site[j][0] == oper_collection[i][0][k] && site[j][1] == oper_collection[i][1][k])//这里需要注意的是总共只有site_num个格点，所以判断该格点在某个跃迁过程是在
                {                                                                                   //oper_collection的前两个后两个就可以了。site中的顺序其实无所谓，site主要就是记住格点
                    particle_num_A[i+1][j] = particle_num_A[i+1][j] - 1;
                    particle_num_B[i+1][j] = particle_num_B[i+1][j] - 1;
                }
                else if(site[j][0] == oper_collection[i][2][k] && site[j][1] == oper_collection[i][3][k])
                {
                    particle_num_A[i+1][j] = particle_num_A[i+1][j] + 1;
                    particle_num_B[i+1][j] = particle_num_B[i+1][j] + 1;
                }
            }
        }
        //如果particle_num为初始填充数，则该态处在基态，标记为0，否则标记为1
        for(int i = 0; i < v+1; i++)
        {
            judge_ground_A = true;
            for(int j = 0; j < sites_num; j++)
            {
                if(abs(site[j][0] + site[j][1])%2 == 0){
                    if(particle_num_A[i+1][j] != gA)
                    {
                        judge_ground_A = false;
                        break;
                    }
                }
                if(abs(site[j][0] + site[j][1])%2 == 1){
                    if(particle_num_A[i+1][j] != gB)
                    {
                        judge_ground_A = false;
                        break;
                    }
                }

            }
            if(judge_ground_A)
                ground_excited_A[i] = 0;
            else
                ground_excited_A[i] = 1;

        }
        for(int i = 0; i < v+1; i++)
        {
            judge_ground_B = true;
            for(int j = 0; j < sites_num; j++)
            {
                if(abs(site[j][0] + site[j][1])%2 == 0){
                    if(particle_num_B[i+1][j] != gB)
                    {
                        judge_ground_B = false;
                        break;
                    }
                }
                if(abs(site[j][0] + site[j][1])%2 == 1){
                    if(particle_num_B[i+1][j] != gA)
                    {
                        judge_ground_B = false;
                        break;
                    }
                }

            }
            if(judge_ground_B)
                ground_excited_B[i] = 0;
            else
                ground_excited_B[i] = 1;

        }
/*
        if(sites_num == 3)
        {

            for(int i = 0; i < v+3; i++)
            {
               for(int j = 0; j < sites_num; j++)
               {
                    cout<<particle_num_A[i][j]<<" ";
               }
               cout<<"   ";
               for(int l = 0; l < sites_num; l++)
               {
                    cout<<particle_num_B[i][l]<<" ";
               }
                cout<<endl;
            }
            cout<<endl;
            for(int i = 0; i < v+1; i++)
                cout<<ground_excited_A[i];
            cout<<endl;
            for(int i = 0; i < v+1; i++)
                cout<<ground_excited_A[i];
            cout<<endl;
            system("pause");
        }
*/
        num_zero1A = 0;
        for(int j = 0; j < v+1; j++)
        {
            if(ground_excited_A[j] == 0)
                num_zero1A += 1;
        }
        num_zero1B = 0;
        for(int j = 0; j < v+1; j++)
        {
            if(ground_excited_B[j] == 0)
                num_zero1B += 1;
        }
        //一行一行的判断，所以顺序不能变
        //生成的跃迁过程和katolist一一比对，如果符合则就保留然后计算。
        //对比的内容是基态的数目
        //如果是对于计算有效场的时候，v次跃迁要对应于v+2阶的katoterm
        for(int i = 0; i < m; i++)
        {
            judge_good_kato_A = true;
            judge_good_kato_B = true;
            num_zero2 = 0;
            for(int j = 0; j < v+1; j++)
                if(eng_katolist[i][j] == 0)
                    num_zero2 += 1;
            if(num_zero1A != num_zero2)
                judge_good_kato_A = false;
            else
            {
                for(int j = 0; j < v+1; j++)
                {
                    if(eng_katolist[i][j] == 0)
                        if(ground_excited_A[j] != eng_katolist[i][j])
                            judge_good_kato_A = false;
                }
            }
            if(num_zero1B != num_zero2)
                judge_good_kato_B = false;
            else
            {
                for(int j = 0; j < v+1; j++)
                {
                    if(eng_katolist[i][j] == 0)
                        if(ground_excited_B[j] != eng_katolist[i][j])
                            judge_good_kato_B = false;
                }
            }
            if(judge_good_kato_A)
            {
                //w和h就是计算的式子，一个是分子一个是分母
                for(int j = 0; j < v+1; j++)
                {
                    wA[j] = 0.0;
                }
                for(int j = 0; j < v+1; j++)
                {

                    //cout<<i<<" "<<j<<" "<<eng_katolist[i][j]<<endl;
                    if(eng_katolist[i][j] == 0)
                    {
                        wA[j] = -1.0;
                    }
                    else
                    {
                        for(int l = 0; l < sites_num; l++)
                        {
                            if(abs(site[l][0] + site[l][1])%2 == 0){

                                wA[j] = wA[j]  + 0.5*(gA*(gA-1)) - (miu+delta_miu)*gA - 0.5*((particle_num_A[j+1][l]-1)*particle_num_A[j+1][l]) + (miu+delta_miu)*particle_num_A[j+1][l];//Em-Ei,m是基态

                            }
                            else if(abs(site[l][0] + site[l][1])%2 == 1){

                                wA[j] = wA[j]  + 0.5*(gB*(gB-1)) - miu*gB - 0.5*((particle_num_A[j+1][l]-1)*particle_num_A[j+1][l]) + miu*particle_num_A[j+1][l];

                            }
                            //cout<<particle_num_A[j+1][l];
                        }
                        //cout<<wA[j]<<endl;
                        /*if(sites_num == 3)
                        {
                            cout<<"A: "<<wA[j]<<endl;
                            system("pause");
                        }*/
                        wA[j] = pow(wA[j], -1.0*eng_katolist[i][j]);
                        //system("pause");

                    }
                }
                for(int j = 0; j < v+2; j++)
                    h[j] = 1.0;
                for(int j = 0; j < v+2; j++)
                {
                    for(int l = 0; l < sites_num; l++)
                    {
                        if(site[l][0] == oper_collection[j][0][k] && site[l][1] == oper_collection[j][1][k])
                        {
                            if(particle_num_A[j][l] > 0)
                                h[j] = h[j] * sqrt(particle_num_A[j][l]*1.0);
                            else
                                h[j] = 0.0;
                        }
                        if(site[l][0] == oper_collection[j][2][k] && site[l][1] == oper_collection[j][3][k])
                        {
                             if(particle_num_A[j][l]+1 > 0)
                                h[j] = h[j] * sqrt((particle_num_A[j][l]+1)*1.0);
                             else
                                h[j] = 0.0;
                        }
                    }
                }
                energy_kato_A = 1.0;
                for(int j = 0; j < v+1; j++)
                {
                    energy_kato_A = energy_kato_A * wA[j];
                }
                for(int j = 0; j < v+2; j++)
                {
                    energy_kato_A = energy_kato_A * h[j];
                }
                energy_kato_A = energy_kato_A * eng_weight[i];
                energy_diagram[0] = energy_diagram[0] + energy_kato_A;
            }
            if(judge_good_kato_B)
            {
                //w和h就是计算的式子，一个是分子一个是分母
                for(int j = 0; j < v+1; j++)
                {
                    wB[j] = 0.0;
                }
                for(int j = 0; j < v+1; j++)
                {
                    //cout<<i<<" "<<j<<" "<<eng_katolist[i][j]<<endl;
                    if(eng_katolist[i][j] == 0)
                    {
                        wB[j] = -1.0;
                    }
                    else
                    {
                        for(int l = 0; l < sites_num; l++)
                        {
                            if(abs(site[l][0] + site[l][1])%2 == 0){
                                wB[j] = wB[j]  + 0.5*(gB*(gB-1)) - miu*gB - 0.5*((particle_num_B[j+1][l]-1)*particle_num_B[j+1][l]) + miu*particle_num_B[j+1][l];//Em-Ei,m是基态

                                //cout<<"B: "<<wB[j]<<endl;

                            }
                            else if(abs(site[l][0] + site[l][1])%2 == 1){

                                wB[j] = wB[j]  + 0.5*(gA*(gA-1)) - (miu+delta_miu)*gA - 0.5*((particle_num_B[j+1][l]-1)*particle_num_B[j+1][l]) + (miu+delta_miu)*particle_num_B[j+1][l];

                                //cout<<"A: "<<wB[j]<<endl;

                            }
                            //cout<<particle_num_B[j+1][l]<<endl;
                        }
                        wB[j] = pow(wB[j], -1.0*eng_katolist[i][j]);


                    }
                }
                for(int j = 0; j < v+2; j++)
                    h[j] = 1.0;
                for(int j = 0; j < v+2; j++)
                {
                    for(int l = 0; l < sites_num; l++)
                    {
                        if(site[l][0] == oper_collection[j][0][k] && site[l][1] == oper_collection[j][1][k])
                        {
                            if(particle_num_B[j][l] > 0)
                                h[j] = h[j] * sqrt(particle_num_B[j][l]*1.0);
                            else
                                h[j] = 0.0;
                        }
                        if(site[l][0] == oper_collection[j][2][k] && site[l][1] == oper_collection[j][3][k])
                        {
                            if(particle_num_B[j][l]+1 > 0)
                                h[j] = h[j] * sqrt((particle_num_B[j][l]+1)*1.0);
                            else
                                h[j] = 0.0;
                        }
                    }
                }
                energy_kato_B = 1.0;
                for(int j = 0; j < v+1; j++)
                {
                    energy_kato_B = energy_kato_B * wB[j];
                }

                for(int j = 0; j < v+2; j++)
                {
                    energy_kato_B = energy_kato_B * h[j];
                }
                energy_kato_B = energy_kato_B * eng_weight[i];
                energy_diagram[1] = energy_diagram[1] + energy_kato_B;
            }
        }
    }

    //-------------------------------------
    //cout<<energy_diagram[0]<<" "<<energy_diagram[1]<<endl;system("pause");
    delete [] ground_excited_A;
    delete [] ground_excited_B;
    delete [] wA;
    delete [] wB;
    delete [] h;
    for(int i = 0; i < v+1; i++)
    {
        delete [] site[i];
        delete [] particle_num_A[i];
        delete [] particle_num_B[i];
    }
    delete [] site;
    delete [] particle_num_A;
    delete [] particle_num_B;
    //-------------------------------------
    return energy_diagram;
}
