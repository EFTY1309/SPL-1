#include<bits/stdc++.h>
using namespace std;
#define N 7

void transpose( double W[N][N],double transpose[N][N], int row,int col)
{
    for (int i = 0; i < row; i++)
    {
      for (int j = 0; j < col; j++)
        {
            transpose[j][i] = W[i][j];
        }
    }
}


void multiplyMatrices( double firstMatrix[N][N],  double secondMatrix[N][N], double mult[N][N], int rowFirst, int columnFirst, int rowSecond, int columnSecond)
{
	int i, j, k;
	for(i = 0; i < rowFirst; ++i)
	{
		for(j = 0; j < columnSecond; ++j)
		{
			mult[i][j] = 0;
		}
	}

	for(i = 0; i < rowFirst; ++i)
	{
		for(j = 0; j < columnSecond; ++j)
		{
			for(k=0; k<columnFirst; ++k)
			{
				mult[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
			}
		}
	}
}

double frobeniusnorm( double y[N][N],int row,int col)
{
    double sum=0,norm;

    for(int i=0;i<row;i++)
        for(int j=0;j<col;j++)
    {
        sum+=pow(y[i][j],2);
    }
    norm=sqrt(sum);
    return norm;
}

void display(double a[N][N],int row,int col)
{
    for(int i=0;i<row;i++){
      for(int j=0;j<col;j++)
      cout<<a[i][j]<<"  ";
      cout<<endl;
    }
}

void ElementwiseMultiplication(double a[N][N],double b[N][N],int row,int col)
{
    for(int i=0;i<row;i++)
        for(int j=0;j<col;j++)
    {
        a[i][j]=(a[i][j])*(b[i][j]);
    }


}

void FindD(double X[N][N],double ab_D[N][N],double D[N][N],int row,int col)
{
    double a1[N];
    for(int i=0;i<row;i++)
        for(int j=0;j<col;j++)
    {
        D[i][j]=0;
    }

    for(int j=0;j<col;j++)
    {
        double power=0;
        for(int i=0;i<row;i++)
        {
            power+=pow((X[i][j]-ab_D[i][j] ),2);
        }
        a1[j]=sqrt(power);
    }


    for(int i=0;i<row;i++)
    {
       for(int j=0;j<col;j++)
        {
            if(i==j) D[i][j]=1/a1[j];
        }
    }

}

int main(void)
{
  double x[N][N];
  int row,col;
  cout<<"Enter the number of row and column for X matrix:"<<endl;
  cin>>row>>col;
  cout<<"Enter the element of the matrices:"<<endl;

  for(int i=0;i<row;i++)
    for(int j=0;j<col;j++)
  {
      cin>>x[i][j];
  }
  double a[N][N];
  int row1,col1;

  cout<<"Enter the number of row and column for A matrix:"<<endl;
  cin>>row1>>col1;
  cout<<"Enter the element of the matrices:"<<endl;

  for(int i=0;i<row1;i++)
    for(int j=0;j<col1;j++)
  {
      cin>>a[i][j];
  }

  double b[N][N];
  int row2,col2;

  cout<<"Enter the number of row and column for B matrix:"<<endl;
  cin>>row2>>col2;
  cout<<"Enter the element of the matrices:"<<endl;

  for(int i=0;i<row2;i++)
    for(int j=0;j<col2;j++)
  {
      cin>>b[i][j];
  }

  int t=200;
  while(t--)
  {
  double a_t[N][N],b_t[N][N],xb_t[N][N],ab[N][N],abb_t[N][N],a_tx[N][N],a_ta[N][N] ,a_tab[N][N],ab1[N][N],L[N][N],D[N][N],ab_d[N][N];
  double a1[N][N],b1[N][N],abD[N][N],abDb_t[N][N],a_txD[N][N],a_tabD[N][N];
  double xD[N][N],xDb_t[N][N];
  transpose(a,a_t,row1,col1);
  transpose(b,b_t,row2,col2);

  multiplyMatrices(a,b,ab_d,row1,col1,row2,col2);

  cout<<"abd matrix is:"<<endl;
  display(ab_d,row,col);

  FindD(x,ab_d,D,row,col);
  cout<<"D matrix is:"<<endl;
  display(D,row,col);

  multiplyMatrices(x,D,xD,row,col,row,col);
  multiplyMatrices(xD,b_t,xDb_t,row,col,col2,row2);
  multiplyMatrices(a,b,ab,row1,col1,row2,col2);
  multiplyMatrices(ab,D,abD,row1,col2,row,col);
  multiplyMatrices(abD,b_t,abDb_t,row1,col,col2,row2);



  for(int i=0;i<row1;i++){
    for(int j=0;j<col1;j++){
        a[i][j]=(a[i][j])*(xDb_t[i][j])/(abDb_t[i][j]);
    }
  }
 // cout<<endl<<"Updated a1 matrix is:"<<endl;
  //display(a1,row1,col1);
  //cout<<endl;

 // ElementwiseMultiplication(a,a1,row1,col1);
   transpose(a,a_t,row1,col1);

  multiplyMatrices(a_t,x,a_tx,col1,row1,row,col);
  multiplyMatrices(a_tx,D,a_txD,col1,col,row,col);

  multiplyMatrices(a_t,a,a_ta,col1,row1,row1,col1);
  multiplyMatrices(a_ta,b,a_tab,col1,col1,row2,col2);
  multiplyMatrices(a_tab,D,a_tabD,col1,col2,row,col);



  cout<<"updated A matrix is:"<<endl;
  display(a,row1,col1);
  cout<<endl;

  for(int i=0;i<row2;i++){
    for(int j=0;j<col2;j++){
        b[i][j]=(b[i][j])*(a_txD[i][j])/(a_tabD[i][j]);
    }
    cout<<endl;
  }
  //cout<<endl<<"Updated b1 matrix is:"<<endl;
  //display(b1,row2,col2);
  //cout<<endl;


  //ElementwiseMultiplication(b,b1,row2,col2);

  cout<<"updated B matrix is:"<<endl;

   display(b,row2,col2);

  multiplyMatrices(a,b,ab1,row1,col1,row2,col2);
  cout<<endl;
  cout<<"Iteration:"<<t<<endl;
  //cout<<endl<<"updated ab matrix is:"<<endl;
  //display(ab1,row1,col2);
  //cout<<endl;

  for(int i=0;i<row;i++)
     for(int j=0;j<col;j++)
  {
      L[i][j]=x[i][j]-ab1[i][j];
  }

  //cout<<endl<<"updated L matrix is:"<<endl;
  //display(L,row,col);
  //cout<<endl;
 cout<<endl<<"The frobenius norm  is:"<<endl;
 double frobenius=frobeniusnorm(L,row,col);
 cout<<endl<<frobenius<<endl;
 cout<<endl;


  }

}

/*5 5
1 3 4 5 6
1 2 4 5 1
1 2 3 5 6
6 7 8 9 10
11 12 13 14 15
5 3
2 3 4
3 4 5
6 7 8
9 10 11
2 4 5
3 5
3 4 5 6 7
1 2 3 4 5
2 7 8 9 10*/


